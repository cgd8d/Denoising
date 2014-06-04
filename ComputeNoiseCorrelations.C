/*
Program which takes in a list of low-background or noise files,
computes the noise correlations observed in the solicited triggers,
and writes those correlations out to a file.

Although there are some cuts on which solicited triggers get used,
generally you should only pass in files in which you expect the solicited trigger
events to be clean.

The output file is in hdf5 format, so it should be portable between very different
machines.  It is also intended to be self-describing, though no guarantee is made
that the file format will remain portable.  The current version is version 1,
so if you change the format, please also increment the version before creating
any files.

This program does the expensive thing and stores noise correlations for all channels;
this makes the files usable for basically all possible purposes, even though most
people will only use the APD noise.  Note that even non-existent channels will
typically be written out.

The program should be called as:
./ComputeNoiseCorrelations <outfilename> <num_events_per_file> <infile1> <infile2 etc>
where <num_events_per_file> can be set to -1 to use all available events.
*/

#include "NoiseCorrelations/NoiseCorrelations.hpp"
#include "NoiseCorrelations/NoiseCorrelationsIO.hpp"
#include "EXOUtilities/EXODimensions.hh"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOEventHeader.hh"
#include "EXOUtilities/EXOCoincidences.hh"
#include "EXOUtilities/EXOTemplWaveform.hh"
#include "EXOUtilities/EXOWaveformFT.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXOFastFourierTransformFFTW.hh"
#include "EXOUtilities/EXORunInfoManager.hh"
#include <vector>
#include <string>
#include <cstdlib>

bool IsEventAcceptable(const EXOEventData* event, const EXOCoincidences& coinc)
{
  // Check that this event is an acceptable example of pure electronic noise.
  const EXOEventHeader& header = event->fEventHeader;

  // Only use solicited triggers.
  if(header.fIndividualTriggerRequest or header.fSumTriggerRequest) return false;

  // Require full-length events.
  if(header.fSampleCount != 2047) return false; // 0 = 1 sample; 2047 = 2048 samples.

  // Don't use events with abnormal tags.
  if(event->fHasSaturatedChannel or header.fTaggedAsNoise or header.fSirenActiveInCR) return false;

  // Reconstruction has the lowest threshold available -- if it found a signal, skip the event.
  // Ignore v-wire signals.
  if(event->GetNumAPDSignals() or event->GetNumUWireSignals()) return false;

  // Skip events occurring during a bad-environment time.
  // This should be sufficient for both noise and LB runs;
  // if we cut events at the beginning and end of the run, that kills a lot of noise-run statistics.
  if(coinc.IsVetoed_BadEnvironment(*event)) return false;

  // OK :)
  return true;
}

int main(int argc, char** argv)
{
  assert(argc > 3);
  const int EntriesPerRun = std::atoi(argv[2]);

  std::vector<int> Runs;
  for(int i = 3; i < argc; i++) Runs.push_back(std::atoi(argv[i]));

  // Build temporary structures to enable rapid calculation of FT products.
  std::map<std::pair<unsigned char, unsigned char>, std::vector<double> > ChanChanProductRR;
  std::map<std::pair<unsigned char, unsigned char>, std::vector<double> > ChanChanProductRI;
  std::map<std::pair<unsigned char, unsigned char>, std::vector<double> > ChanChanProductIR;
  std::map<std::pair<unsigned char, unsigned char>, std::vector<double> > ChanChanProductII;
  for(unsigned char i = 0; i < (unsigned char)NUMBER_READOUT_CHANNELS; i++) {
    for(unsigned char j = i; j < (unsigned char)NUMBER_READOUT_CHANNELS; j++) {
      ChanChanProductRR[std::make_pair(i, j)].resize(1024);
      ChanChanProductRI[std::make_pair(i, j)].resize(1024);
      ChanChanProductIR[std::make_pair(i, j)].resize(1024);
      ChanChanProductII[std::make_pair(i, j)].resize(1024);
    }
  }

  size_t NumEntriesAccepted = 0;

  for(size_t runIndex = 0; runIndex < Runs.size(); runIndex++) {
    std::cout<<"Starting on run "<<Runs[runIndex]<<" ("<<runIndex<<"/"<<Runs.size()<<")"<<std::endl;

    const EXORunInfo::RunList& rawRunList =
      EXORunInfoManager::GetRunInfo(Runs[runIndex], "Data/Raw/root").GetRunFiles();
    const EXORunInfo::RunList& procRunList =
      EXORunInfoManager::GetRunInfo(Runs[runIndex], "Data/Processed/masked").GetRunFiles();
    if(rawRunList.size() != procRunList.size()) {
      std::cout<<"Run "<<Runs[runIndex]<<" has a mismatch between raw and processed file count."<<std::endl;
      continue;
    }

    TChain rawChain("tree");
    TChain procChain("tree");
    for(EXORunInfo::RunList::const_iterator it = rawRunList.begin(); it != rawRunList.end(); it++) {
      rawChain.AddFile(it->GetFileLocation().c_str());
    }
    for(EXORunInfo::RunList::const_iterator it = procRunList.begin(); it != procRunList.end(); it++) {
      procChain.AddFile(it->GetFileLocation().c_str());
    }

    EXOEventData* RawEvent = NULL;
    EXOEventData* ProcEvent = NULL;
    rawChain.SetBranchAddress("EventBranch", &RawEvent);
    procChain.SetBranchAddress("EventBranch", &ProcEvent);
    assert(rawChain.BuildIndex("fRunNumber", "fEventNumber") >= 0);

    EXOCoincidences coinc;
    coinc.Load(procChain); // Clears any previously-loaded chain.

    int NumAcceptedFromThisRun = 0;
    Long64_t EntryNum = 0;
    while((NumAcceptedFromThisRun < EntriesPerRun or EntriesPerRun == -1) and
          EntryNum < procChain.GetEntries()) {
      if(EntryNum % 1000 == 0) std::cout<<"\tTrying entry "<<EntryNum<<std::endl;
      procChain.GetEntry(EntryNum);
      EntryNum++;
      if(not IsEventAcceptable(ProcEvent, coinc)) continue;

      // Get the raw entry by index, since entry numbers won't generally match (due to masking).
      assert(rawChain.GetEntryWithIndex(ProcEvent->fRunNumber, ProcEvent->fEventNumber));

      // Convert EXOIntWaveforms to EXODoubleWaveforms.
      // Also establish the channel ordering.
      std::vector<EXOWaveformFT> FourierWaveforms;
      for(Int_t channel = 0; channel < NUMBER_READOUT_CHANNELS; channel++) {
        EXOWaveform* rawWF = RawEvent->GetWaveformData()->GetWaveformWithChannelToEdit(channel);
        if(rawWF) {
          rawWF->Decompress();
          EXODoubleWaveform dblWF = *rawWF;
          EXOWaveformFT ftWF;
          assert(dblWF.GetLength() == 2048);
          EXOFastFourierTransformFFTW::GetFFT(2048).PerformFFT(dblWF, ftWF);
          assert(ftWF.GetLength() == 1025); // Just to make sure I'm not confused.
          FourierWaveforms.push_back(ftWF);
        }
        else {
          // This channel doesn't exist in this run -- but don't leave an empty space,
          // since that introduces the nuisance of including a mapping from index to software channel.
          // Instead, just fill with zeros.
          EXOWaveformFT tempFT;
          tempFT.SetLength(1025);
          tempFT.Zero();
          FourierWaveforms.push_back(tempFT);
        }
      }
      assert(FourierWaveforms.size() == NUMBER_READOUT_CHANNELS);

      // Take products of the vectors, and add them to the noise correlations container.
      // Note that we include all channels, but not the f=0 component.
      for(unsigned char i = 0; i < (unsigned char)NUMBER_READOUT_CHANNELS; i++) {
        for(unsigned char j = i; j < (unsigned char)NUMBER_READOUT_CHANNELS; j++) {
          const EXOWaveformFT& ft1 = FourierWaveforms[i];
          const EXOWaveformFT& ft2 = FourierWaveforms[j];
          std::vector<double>& outRR = ChanChanProductRR.find(std::make_pair(i, j))->second;
          std::vector<double>& outRI = ChanChanProductRI.find(std::make_pair(i, j))->second;
          std::vector<double>& outIR = ChanChanProductIR.find(std::make_pair(i, j))->second;
          std::vector<double>& outII = ChanChanProductII.find(std::make_pair(i, j))->second;

          for(size_t f = 0; f < 1024; f++) {
            outRR[f] += ft1[f+1].real()*ft2[f+1].real();
            outRI[f] += ft1[f+1].real()*ft2[f+1].imag();
            outIR[f] += ft1[f+1].imag()*ft2[f+1].real();
            outII[f] += ft1[f+1].imag()*ft2[f+1].imag();
          }
        }
      }

      NumAcceptedFromThisRun++;
      NumEntriesAccepted++;
    } // Finish loop over entries in a particular run.
  } // Finish loop over runs.

  // Convert ChanChanProduct to a NoiseCorrelations object.
  MapIndexHandler<unsigned char> ChannelIndex;
  for(unsigned char i = 0; i < (unsigned char)NUMBER_READOUT_CHANNELS; i++) ChannelIndex.InsertKey(i);
  assert(ChannelIndex.MaxIndex() == NUMBER_READOUT_CHANNELS);
  NoiseCorrelations NoiseCorr(ChannelIndex);
  assert(NoiseCorr.GetFrequencyIndex().MaxIndex() == 1024);

  for(size_t i = 0; i < NUMBER_READOUT_CHANNELS; i++) {
    for(size_t j = i; j < NUMBER_READOUT_CHANNELS; j++) {
      std::pair<unsigned char, unsigned char> tmpChanPair(i, j);
      std::vector<double>& oldFormRR = ChanChanProductRR.find(tmpChanPair)->second;
      std::vector<double>& oldFormRI = ChanChanProductRI.find(tmpChanPair)->second;
      std::vector<double>& oldFormIR = ChanChanProductIR.find(tmpChanPair)->second;
      std::vector<double>& oldFormII = ChanChanProductII.find(tmpChanPair)->second;

      for(size_t f = 0; f < 1024; f++) {

        // Normalize for number of entries here.
        oldFormRR[f] /= NumEntriesAccepted;
        oldFormRI[f] /= NumEntriesAccepted;
        oldFormIR[f] /= NumEntriesAccepted;
        oldFormII[f] /= NumEntriesAccepted;

        NoiseMatrix& mat = NoiseCorr.GetMatrixForIndex(f);

        // real-real
        mat.GetCorrByIndex(i, j) = oldFormRR[f];
        mat.GetCorrByIndex(j, i) = oldFormRR[f];

        // imag-imag
        mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j+NUMBER_READOUT_CHANNELS) = oldFormII[f];
        mat.GetCorrByIndex(j+NUMBER_READOUT_CHANNELS, i+NUMBER_READOUT_CHANNELS) = oldFormII[f];

        // real-imag and imag-real
        // if i == j, oldForm.IR == oldForm.RI.
        mat.GetCorrByIndex(i, j+NUMBER_READOUT_CHANNELS) = oldFormRI[f];
        mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j) = oldFormIR[f];
        mat.GetCorrByIndex(j, i+NUMBER_READOUT_CHANNELS) = oldFormIR[f];
        mat.GetCorrByIndex(j+NUMBER_READOUT_CHANNELS, i) = oldFormRI[f];
      }
    }
  }

  // Apply symmetries and identities as appropriate, to ensure the best use of redundant information.
  // Note that we've already ensured i <-> j symmetry, but we must be careful to retain it.
  for(size_t f = 0; f < 1023; f++) { // symmetries which don't apply to the last frequency.
    NoiseMatrix& mat = NoiseCorr.GetMatrixForIndex(f);
    for(size_t i = 0; i < NUMBER_READOUT_CHANNELS; i++) {
      for(size_t j = i; j < NUMBER_READOUT_CHANNELS; j++) {
        double tmp;

        // real*real = imag*imag.
        tmp = mat.GetCorrByIndex(i, j) +
              mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j+NUMBER_READOUT_CHANNELS);
        tmp /= 2;
        mat.GetCorrByIndex(i, j) = tmp;
        mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j+NUMBER_READOUT_CHANNELS) = tmp;
        mat.GetCorrByIndex(j, i) = tmp;
        mat.GetCorrByIndex(j+NUMBER_READOUT_CHANNELS, i+NUMBER_READOUT_CHANNELS) = tmp;

        // real*imag = -imag*real.
        tmp = mat.GetCorrByIndex(i, j+NUMBER_READOUT_CHANNELS) -
              mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j);
        tmp /= 2;
        if(i == j) tmp = 0; // Anti-symmetry means entries are zero on the diagonal.
        mat.GetCorrByIndex(i, j+NUMBER_READOUT_CHANNELS) = tmp;
        mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j) = -tmp;
        mat.GetCorrByIndex(j, i+NUMBER_READOUT_CHANNELS) = -tmp;
        mat.GetCorrByIndex(j+NUMBER_READOUT_CHANNELS, i) = tmp;
      }
    }
  }
  // Next handle the last f, f=1023.
  // All imaginary terms are perfectly zero.
  // Real-real term symmetry is already guaranteed, so nothing we can do about them.
  {
    NoiseMatrix& mat = NoiseCorr.GetMatrixForIndex(1023);
    for(size_t i = 0; i < NUMBER_READOUT_CHANNELS; i++) {
      for(size_t j = i; j < NUMBER_READOUT_CHANNELS; j++) {
        mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j) = 0;
        mat.GetCorrByIndex(i, j+NUMBER_READOUT_CHANNELS) = 0;
        mat.GetCorrByIndex(i+NUMBER_READOUT_CHANNELS, j+NUMBER_READOUT_CHANNELS) = 0;
        mat.GetCorrByIndex(j+NUMBER_READOUT_CHANNELS, i) = 0;
        mat.GetCorrByIndex(j, i+NUMBER_READOUT_CHANNELS) = 0;
        mat.GetCorrByIndex(j+NUMBER_READOUT_CHANNELS, i+NUMBER_READOUT_CHANNELS) = 0;
      }
    }
  }

  // Finally, write to file.
  std::cout<<"Writing to file."<<std::endl;
  NoiseCorrelationIO::WriteNoiseCorrelations(argv[1], NoiseCorr);
  std::cout<<"Done writing file."<<std::endl;

  std::cout<<NumEntriesAccepted<<" entries were used in creating this noise correlation file."<<std::endl;
}
