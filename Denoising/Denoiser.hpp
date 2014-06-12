#ifndef Denoiser_hpp
#define Denoiser_hpp

/*
Class to manage denoising.

Set fDenoiseAPDs or fDenoiseUWires to decide which things to denoise.  (The default is APDs only.)

Then call InitDenoiser with the first event to denoise; this will set up a channel list, create
the lightmap and noise matrices, possibly other stuff.  The channel list is also returned.

Then...

*/


class Denoiser
{


 public:

  Denoiser()
  : fDenoiseAPDs(true),
    fDenoiseUWires(false),
    fNumMulsAccumulate(500),
  { }

  // Parameters to set *before* calling InitDenoiser.
  bool fDenoiseAPDs;
  bool fDenoiseUWires;
  std::string fLightmapFilename;
  std::string fNoiseCorrFilename;
  size_t fNumMulsAccumulate;

  // Initialize the denoiser with the first event.
  // Return the channel map which we should use.
  MapIndexHandler<unsigned char> InitDenoiser(const EXOEventData* event) {
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

    // Read in the noise correlations.
    fNoiseCorr.SetChannelIndex(fChannelMap);
    NoiseCorrelationsIO::ReadNoiseCorrelations(fNoiseCorrFilename, fNoiseCorr);

    return fChannelMap;
  }


 private:

  LightMap::PositionFunc fAPDPosFunc;
  LightMap::GainSnapshot fAPDGainSnapshot;  
  NoiseCorrelations fNoiseCorr;
  MapIndexHandler<unsigned char> fChannelMap; // Channels which we are denoising.
};
#endif
