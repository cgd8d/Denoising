/*
IN PROGRESS

Main program to read in binary data files, reconstruct it, and then do denoising.

The first argument should be a filename for an exo instruction file.


Thoughts for self:

One option is to break out the inner loop of EXOAnalysisManager::RunAnalysis into a new
function, EXOAnalysisManager::StepAnalysis, which only does one event *and returns it*.
This would take EXOAnalysisManager out of the loop entirely.

Another option would be to optionally let myself install a callback in EXOAnalysisManager which gets invoked
at the end of the event loop.  It would need to have a fixed known signature -- for instance, I would need
to decide whether the callback function is in a class or not.

A slightly more flexible option would be to create a callback module which accepts as arguments
one or more callback functions.  This gives me the flexibility to do callback at any point in the
processing sequence, but I probably don't need that flexibility.

Let's go with the simple callback option.

*/




int main(int argc, char** argv)
{







}
