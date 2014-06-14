#ifndef Denoiser_hpp
#define Denoiser_hpp

/*
Class to manage denoising.

Set fDenoiseAPDs or fDenoiseUWires to decide which things to denoise.  (The default is APDs only.)

Then call InitDenoiser with the first event to denoise; this will set up a channel list, create
the lightmap and noise matrices, possibly other stuff.

Then...

*/

#include "Denoising/WaveformCache.hpp"

class Denoiser
{


 public:

  Denoiser(EXOSetDenoisedScintModule& SetDenoisedScintModule)
  : fDenoiseAPDs(true),
    fDenoiseUWires(false),
    fNumMulsAccumulate(500),
    fReclusterDistForScint(1.*CLHEP::meter),
    fReclusterDistForCharge(1.*CLHEP::centimeter),
    fReclusterTimeForCharge(1.*CLHEP::microsecond),
    fSetDenoisedScintMod(SetDenoisedScintMod),
  { }

  // Parameters to set *before* calling InitDenoiser.
  bool fDenoiseAPDs;
  bool fDenoiseUWires;
  std::string fLightmapFilename;
  std::string fNoiseCorrFilename;
  size_t fNumMulsAccumulate;
  double fReclusterDistForScint;
  double fReclusterDistForCharge;
  double fReclusterTimeForCharge;

  // Initialize the denoiser with the first event.
  // Return the channel map which we should use.
  void InitDenoiser(const EXOEventData* event) {
    assert(event);

    // Get the channel map for this event.
    const EXOChannelMap& ChannelMap = GetChanMapForHeader(ED->fEventHeader);

    // Pick channels to use.
    // Don't use bad channels; we assume this list does not change event-by-event.
    // Also, here we decide whether to denoise u-wires.
    assert(fDenoiseAPDs or fDenoiseUWires);
    assert(fChannelMap.MaxIndex() == 0);
    for(size_t chan = 0; chan < NUMBER_READOUT_CHANNELS; chan++) {
      if(ChannelMap.channel_suppressed_by_daq(i) or not ChannelMap.good_channel(i)) continue;

      EXOMiscUtil::EChannelType channelType = EXOMiscUtil::TypeOfChannel(i);
      if(not(fDenoiseAPDs and channelType == EXOMiscUtil::kAPDGang or
             fDenoiseUWires and channelType == EXOMiscUtil::kUWire)) continue;

      fChannelMap.InsertKey(chan);
    }

    // If we're denoising APDs (as we generally are), read in the lightmap.
    if(fDenoiseAPDs) {
      fAPDGainSnapshot = LightMapIO::ReadLightmapAtRun(fLightmapFilename,
                                                       event->fRunNumber,
                                                       fAPDPosFunc);
    }

    // Initialize the noise multiplier.
    // fNoiseMultiplier.fNumVectorsMax = 1000; // Set here to allocate bigger buffers.
    fNoiseMultiplier.FillNoiseCorrelations(fNoiseCorrFilename, fChannelMap);

    // Let the waveform cacher know which channels are worth caching.
    fWFCache.SetChannelIndex(fChannelMap);
  }

  void StepDenoiser(const EXOEventData* event) {
    assert(event);

    // Go ahead and create an EventContainer.
    // True, we don't know whether this event is denoisable, but there are only so many events total.
    // Note that shared_ptr should be created *exactly* this way to be exception-safe.
    boost::shared_ptr<EventContainer> evtContainer(new EventContainer);

    // Form a list of the scintillation clusters we'll denoise, and how we'll recluster them.
    if(fDenoiseAPDs) {
      for(size_t iscint = 0; iscint < event->GetNumScintillationClusters(); iscint++) {
        EXOScintillationCluster* scint = event->GetScintillationCluster(iscint);
        std::vector<std::set<std::pair<EXOChargeCluster*, size_t> > > charge_reclustered;
        for(size_t icharge = 0; icharge < scint->GetNumChargeClusters(); icharge++) {
          EXOChargeCluster* charge = scint->GetChargeClusterAt(icharge);

          // We can't make use of a charge cluster without well-defined position,
          // because we wouldn't know where to evaluate the lightmap.
          if(std::fabs(charge->fX) > 200*CLHEP::millimeter or
             std::fabs(charge->fY) > 200*CLHEP::millimeter or
             std::fabs(charge->fZ) > 200*CLHEP::millimeter) continue;

          // If the charge cluster is estimated to have no energy, then we can't use it either;
          // however, we can defer to the energy estimate without a purity correction for now
          // knowing that the run will be marked as missing a calibration.
          // If neither of those is positive, though, then give up on this cluster.
          if(charge->fPurityCorrectedEnergy < 1 and charge->fCorrectedEnergy < 1) continue;

          // If we've made it to this point, we want to include this cluster as a location for scintillation.
          // Add it to charge_reclustered, either to an existing reclustering or a new one.
          Recluster(charge_reclustered, charge, icharge,
                    "3DPos", fReclusterDistForScint);
        }

        for(size_t i = 0; i < charge_reclustered.size(); i++) {
          ScintDepositSummary summary;
          summary.fScintCluster = iscint;
          for(std::set<std::pair<EXOChargeCluster*, size_t> >::iterator it = charge_reclustered[i].begin();
              it != charge_reclustered[i].end();
              it++) {
            summary.fIncludedChargeClusters.push_back(it->second);
          }
          evtContainer->fEventSummary.fScintDeposits.push_back(summary);
        }
      } // End loop over scintillation clusters.
    } // End if(fDenoiseAPDs).

    // Form a list of the charge clusters we'll denoise, and how we'll recluster them.
    if(fDenoiseUWires) {
      std::vector<std::set<std::pair<EXOChargeCluster*, size_t> > > charge_reclustered;
      for(size_t icharge = 0; icharge < event->GetNumChargeClusters(); icharge++) {
        EXOChargeCluster* charge = event->GetChargeCluster(icharge);

        // Since we currently only use u-wires, we only need a well-defined
        // u-position and collection time.
        if(std::fabs(charge->fU) > 200*CLHEP::millimeter or
           charge->fCollectionTime < 0 or
           charge->fCollectionTime > 2048*CLHEP::microsecond) continue;

        // If the charge cluster is estimated to have no energy, then call it suspicious;
        // however, we can defer to the energy estimate without a purity correction for now
        // knowing that the run will be marked as missing a calibration.
        // If neither of those is positive, though, then give up on this cluster.
        if(charge->fPurityCorrectedEnergy < 1 and charge->fCorrectedEnergy < 1) continue;

        // If we've made it to this point, we want to include this cluster as a location for ionization.
        // Add it to charge_reclustered, either to an existing reclustering or a new one.
        Recluster(charge_reclustered, charge, icharge,
                  "UandT", fReclusterDistForCharge, fReclusterTimeForCharge);
      }

      for(size_t i = 0; i < charge_reclustered.size(); i++) {
        ChargeDepositSummary summary;
        for(std::set<std::pair<EXOChargeCluster*, size_t> >::iterator it = charge_reclustered[i].begin();
            it != charge_reclustered[i].end();
            it++) {
          summary.fIncludedChargeClusters.push_back(it->second);
        }
        evtContainer->fEventSummary.fChargeDeposits.push_back(summary);
      }
    } // End if(fDenoiseUWires).

    // If we have no signals to denoise, then we're done.
    if(evtContainer->EventSummary.fScintDeposits.size() +
       evtContainer->EventSummary.fChargeDeposits.size() == 0) {
      fSetDenoisedScintMod.InsertResult(event->fRunNumber, event->fEventNumber, evtContainer->EventSummary);
      return;
    }

    // Now build the model waveforms for the deposits we want to denoise.














 private:

  LightMap::PositionFunc fAPDPosFunc;
  LightMap::GainSnapshot fAPDGainSnapshot;

  NoiseMultiplier fNoiseMultiplier;
  WaveformCache fWFCache;
  EXOSetDenoisedScintModule& fSetDenoisedScintMod;

  MapIndexHandler<unsigned char> fChannelMap; // Channels which we are denoising.

  // Convenience function to handle the logic of reclustering.
  // We follow the definition which leads to a unique reclustering -- for a distance D,
  // if a pair of clusters are within a distance D then they are placed within the same clustering.
  // That means that a reclustering can be quite large, but it is uniquely defined.
  static void Recluster(std::vector<std::set<EXOChargeCluster*> >& existing_clusters,
                        const EXOChargeCluster* new_charge,
                        size_t PosOfNewCharge,
                        std::string TypeOfDistance,
                        double Distance,
                        double Time = -1) {
    assert(TypeOfDistance == "3DPos" or TypeOfDistance == "UandT");
    if(TypeOfDistance == "UandT") assert(Time > 0);
    // See if we should recluster this charge with other charge from the same scint.
    bool is_merged = false;
    for(size_t i = 0; i < existing_clusters.size(); i++) {
      for(std::set<std::pair<EXOChargeCluster*, size_t> >::iterator it = existing_clusters[i].begin();
          it != existing_clusters[i].end();
          it++) {
        bool mergeable = true;
        if(TypeOfDistance == "3DPos") {
          double dx = new_charge->fX - it->first->fX;
          double dy = new_charge->fY - it->first->fY;
          double dz = new_charge->fZ - it->first->fZ;
          if(dx*dx + dy*dy + dz*dz >= Distance*Distance) mergeable = false;
        }
        else {
          assert(TypeOfDistance == "UandT" and Time > 0) {
          double du = new_charge->fU - it->first->fU;
          if(dU*dU >= Distance*Distance) mergeable = false;
          double dT = new_charge->fCollectionTime - it->first->fCollectionTime;
          if(dT*dT >= Time*Time) mergeable = false;
        }
        if(mergeable) {
          // Recluster charge with this group.
          existing_clusters[i].insert(std::make_pair(new_charge, PosOfNewCharge));
          is_merged = true;
          // Also verify that no other charge reclusterings get merged with this one too.
          for(size_t j = existing_clusters.size()-1; j > i; j--) {
            for(std::set<std::pair<EXOChargeCluster*, size_t> >::iterator it2 = existing_clusters[j].begin();
                it2 != existing_clusters[j].end();
                it2++) {
              mergeable = true;
              if(TypeOfDistance == "3DPos") {
                double dx = new_charge->fX - it2->first->fX;
                double dy = new_charge->fY - it2->first->fY;
                double dz = new_charge->fZ - it2->first->fZ;
                if(dx*dx + dy*dy + dz*dz >= Distance*Distance) mergeable = false;
              }
              else {
                assert(TypeOfDistance == "UandT" and Time > 0) {
                double du = new_charge->fU - it2->first->fU;
                if(dU*dU >= Distance*Distance) mergeable = false;
                double dT = new_charge->fCollectionTime - it2->first->fCollectionTime;
                if(dT*dT >= Time*Time) mergeable = false;
              }
              if(mergeable) {
                // Merge these two sets together.
                existing_clusters[i].insert(existing_clusters[j].begin(),
                                            existing_clusters[j].end());
                // Erase entry j, which has been merged into entry i.
                if(j < existing_clusters.size()-1) {
                  existing_clusters[j] = existing_clusters.back();
                }
                existing_clusters.pop_back();
                break; // Step back and examine another set.
              }
            }
          } // End verifying that there are no other charge clusterings to merge into i.
          break;
        } // End if(dx*dx + dy*dy + dz*dz < Distance*Distance).
      } // End loop through existing_clusters[i].
      if(is_merged) break;
    } // End for(size_t i = 0; i < charge_reclustered.size(); i++).

    // If we didn't merge this charge cluster with any others, then it forms its own new grouping.
    if(not is_merged) {
      existing_clusters.resize(existing_clusters.size()+1);
      existing_clusters.back().insert(std::make_pair(new_charge, PosOfNewCharge));
    }
  }



};
#endif
