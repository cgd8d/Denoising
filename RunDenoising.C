/*
IN PROGRESS

Main program to read in binary data files, reconstruct it, and then do denoising.

The first argument should be a filename for an exo instruction file.  The toutput module
should be used, but the output filename should not be specified.

Invoke with:

RunDenoising <exo_filename> <temp_output_filename> <final_output_filename>

The exo file should include use of the toutput module, but should *not* specify
the filename of the temporary output file; we take care of that.
It also should not invoke begin or exit.

*/

#include "Denoising/EXOSetDenoisedScintModule.hpp"

#include "EXOAnalysisManager/EXOAnalysisModule.hh"
#include "EXOAnalysisManager/EXOAnalysisManager.hh"
#include "EXOUtilities/EXOTalkToManager.hh"
#include "EXOUtilities/EXOEventData.hh"

#include <cassert>
#include <iostream>


int main(int argc, char** argv)
{
  assert(argc == 3);

  // Create the object which will retain, and later set, denoised results.
  EXOSetDenoisedScintModule SetDenoisedScintMod;

  {
    // Start up an EXOAnalysis session.  Initialize it with the settings from the exo.
    EXOTalkToManager talktoManager;
    talktoManager.SetFilename(argv[1]);
    EXOAnalysisManager analysisManager(&talktoManager);
    talktoManager.InterpretStream();

    // Set the output filename for a temporary output file.
    // NOTE: If we migrate to root 5.34, more efficient would be to use a TMemFile.
    talktoManager.InterpretCommand(std::string("/toutput/file ") + argv[2]);

    // Initialize analysis.
    analysisManager.InitAnalysis();

    // Run processing loop.
    EXOEventData* event = analysisManager.StepAnalysis();
    //if(event) InitDenoiser(event); // Get appropriate lightmap, etc.
    while(event) {
      //bool EventIsDenoisable = StepDenoiser(event);
      //if(EventIsDenoisable) CacheWaveforms(event);
      event = analysisManager.StepAnalysis();
    }

    // Finish EXOAnalysis processing.
    if(analysisManager.EndOfRunSegment(eventData) < 0 or
       analysisManager.EndOfRun(eventData) < 0) {
      throw EXOAnalysisManager::FailedProcess("Failed finishing up processing");
    }
    analysisManager.ShutDown();

    // Print running statistics
    std::cout << "EXOAnalysis module statistics (excluding denoising):" << std::endl;
    analysisManager.print_statistics();
  }

  // Reset shared objects so that we can do another loop of processing.
  EXOAnalysisModule::ClearSharedObjects();

  // Denoising is done; just use tinput/toutput to fill the scintillation variables.
  {
    // Start up an EXOAnalysis session.  Initialize it with the settings from the exo.
    EXOTalkToManager talktoManager;
    EXOAnalysisManager analysisManager(&talktoManager);

    // Specify modules to use.
    analysisManager.UseModule("tinput");
    analysisManager.UseModule(&SetDenoisedScintMod);
    analysisManager.UseModule("toutput");

    // Pass commands to tinput and toutput.
    talktoManager.InterpretCommand(std::string("/input/file ") + argv[2]);
    talktoManager.InterpretCommand(std::string("/toutput/file") + argv[3]);

    // Run.
    talktoManager.InterpretCommand("begin");

    // Print running statistics
    std::cout << "EXOAnalysis module statistics, part two (excluding denoising):" << std::endl;
    analysisManager.print_statistics();
  }
}
