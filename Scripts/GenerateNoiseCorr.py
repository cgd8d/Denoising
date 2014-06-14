'''
Generate noise correlation files for each of the periods we've identified as having constant noise.
This will result in the submission of a number of batch jobs.
Note that exactly where the boundary between time periods occurs is not always 100% certain,
but where possible we have tried to identify a cause to help with the precise placement.

In the future, one way to improve the placement could be to try shifting them slightly (with all else
held constant) and see if the resolution changes.  Another cool project could be to design
a criterion where a run's noise is or is not consistent with the average of its window.  Note one
difficulty that we currently extract noise only from physics run and (where available) noise runs; but
there is noise information in source and charge injection runs, so we could try to use that to really
pinpoint when the noise changes.  Also, we deal in units of one run, but the noise is
temperature-dependent, and no one has tried to identify that.

Finally, there is an issue that we have a type of noise which is correlated with the trigger time,
referred as "glitch noise".  It has changed magnitude over time, and we don't understand its origin.
Samuel Homiller developed a short script which first identifies the average glitch noise of a run by
taking the average of many noise traces, then subtracts it off before computing the noise correlations
of that channel.  It has not been tested for denoising, but it would be nice to understand if it makes a
difference.
'''

import subprocess
import os
import shutil
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libEXOUtilities")

EntriesPerRun = -1
Queue = 'xxl'

# Inclusive run ranges.
# This set of windows is meant to only include Run2ab for now.
# Boundaries which are explicitly commented are confirmed exactly;
# others are not independently confirmed by me with environmental correlations.
# (I still trust them, though.)
RunWindows = [(2424, 2690), # FEC voltage adjustment
              (2691, 2852), # Ebox 1 fan installed
              (2853, 2891), # Ebox 2 fan installed
              (2892, 3117), # Power outage here.
              (3118, 3326), # APD board swap
              (3327, 3700), # Utility power swap
              (3701, 3949),
              (3950, 4140), # Ralph's diode box installed
              (4141, 4579),
              (4580, 4779),
              (4780, 5197), # LC filters removed from FECs
              (5198, 5590), # Toward end of 5590 CR temps are elevated; there was a lasting effect
              (5591, 5892)] # Run2c ends.

# Make a directory for log files.
LogFileDir = 'Tmp/NoiseCorrLogFiles'
try: shutil.rmtree(LogFileDir)
except OSError: pass
os.makedirs(LogFileDir)

for runWindow in RunWindows:
    ds = ROOT.EXORunInfoManager.GetDataSet('Data/Processed/masked', '%i<=run&&run<=%i&&quality==\"GOLDEN\"&&runType==\"Data-Physics\"' % runWindow)
    runList = [str(ri.GetRunNumber()) for ri in ds]
    subprocess.call(['bsub', '-q', Queue, '-R', 'rhel60',
                     '-o', os.path.join(LogFileDir, '%i_to_%i.log' % runWindow),
                     './ComputeNoiseCorrelations',
                     'Tmp/NoiseCorr_%i_to_%i.hdf5' % runWindow, str(EntriesPerRun)] + runList)

# We also generate a noise window for runs 2401-2423 (09-28-11 APD biases).
# But, lacking physics data, we use a noise run there.
# Always use the full run because it's all we've got.
subprocess.call(['bsub', '-q', Queue, '-R', 'rhel60',
                 '-o', os.path.join(LogFileDir, '2401_to_2423.log'),
                 './ComputeNoiseCorrelations',
                 'Tmp/NoiseCorr_2401_to_2423.hdf5', '-1', '2401'])

