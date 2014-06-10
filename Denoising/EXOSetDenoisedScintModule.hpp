#ifndef EXOSetDenoisedScintModule_hpp
#define EXOSetDenoisedScintModule_hpp

/*
A module which accepts denoised results, stores them, and later installs them into events during processing.

Currently it does not handle charge clusters, but in mainline code we don't try to denoised charge
so it should be fine for now.
*/

#include "EXOAnalysisManager/EXOAnalysisModule.hh"
#include "EXOAnalysisManager/EXOAnalysisModuleFactory.hh"

class EXOSetDenoisedScintModule : public EXOAnalysisModule
{
 public:

  // Result for an event (with possibly many or no scintillation clusters).
  struct DenoisingResultSet {
    int fResultSetCode;
    std::vector<double> fResults;
  };

  EventStatus ProcessEvent(EXOEventData *ED) {
    assert(ED);
    if(ED->GetNumScintillationClusters() == 0) {
      // Nothing to set, so return.
      return kOk;
    }

    DenoisingResultsT::iterator it =
      fDenoisingResults.find(std::make_pair<int, int>(ED->fRunNumber, ED->fEventNumber));
    assert(it != fDenoisingResults.end()); // We should ALWAYS have some kind of result.

    if(it->fResultSetCode != 0) {
      // We failed to even attempt denoising.
      // Install the error code for all scintillation clusters,
      // and set the scintillation energy to 0.
      for(size_t i = 0; i < ED->GetNumScintillationClusters(); i++) {
        EXOScintillationCluster* scint = ED->GetScintillationCluster(i);
        scint->fDenoisedEnergy = scint->fRawEnergy = scint->fEnergy = 0;
        scint->fDenoisedEnergyError = 0;
        scint->fDenoisingInternalCode = it->fResultSetCode;
      }
    }
    else {
      // Denoising was successful.
      assert(ED->GetNumScintillationClusters() == it->fResults.size());
      for(size_t i = 0; i < ED->GetNumScintillationClusters(); i++) {
        EXOScintillationCluster* scint = ED->GetScintillationCluster(i);
        scint->fDenoisedEnergy = scint->fRawEnergy = scint->fEnergy = it->fResults[i];
        scint->fDenoisedEnergyError = 0;
        scint->fDenoisingInternalCode = 0;
      }
    }

    return kOk;
  }


 private:

  typedef std::map<std::pair<int, int>,
          DenoisingResultSet> DenoisingResultsT;

  // Map from run/event number to a list of scintillation energies.
  // At the end of denoising, absence from this list means that there was nothing to denoise.
  // If denoising fails, we should add an entry indicating why.
  DenoisingResultsT fDenoisingResults;

  DEFINE_EXO_ANALYSIS_MODULE(EXOSetDenoisedScintModule)
};

// No need for IMPLEMENT_EXO_ANALYSIS_MODULE, since we won't build this through the factory.




#endif
