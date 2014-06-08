# Run from the base directory with:
# bsub -q xlong -R rhel60 -o RunComputeRotationAngle.log python Scripts/RunComputeRotationAngle.py

import subprocess
import shutil
import os
import time
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libEXOUtilities")

# Grab all source datasets from the data catalog.
# Note we use a few suspect runs, so we have to run this script over those too.
# Note: is there a filter, or could there be one, for picking out Th runs specifically?
ds = ROOT.EXORunInfoManager.GetDataSet("Data/Processed/masked", "quality!=\"BAD\"&&run>=2401&&run<=5892&&runType==\"Data-Source calibration\"")

# Create the needed directories, and ensure they are empty.
OutDir = 'Tmp/ComputeRotationAngle_oldversion'
try: shutil.rmtree(OutDir)
except OSError: pass
os.makedirs(OutDir)

# For each run in parallel, extract the single-site full-energy Th events into a database.
proc = []
for runInfo in ds:
    runNo = runInfo.GetRunNumber()
    # Try to clear out processes; simultaneously limit the number of running processes.
    while len(proc) > 30:
        map(lambda x: x[0].poll(), proc)
        finishedProcs = filter(lambda x: x[0].returncode != None, proc)
        for p in finishedProcs:
            proc.remove(p)
            if p[0].returncode != 0: print "Run %08i failed." % p[1]
        if len(proc) > 30: time.sleep(10) # Just to avoid burning cycles for no reason.

    if runNo == 4435:
        # This is a Cs run misidentified as Th.  (Need to get Tony to patch this.)
        continue
    # Submit the job.
    proc.append((subprocess.Popen(['bsub', '-q', 'xlong', '-R', 'rhel60', '-K',
                                   '-o', os.path.join(OutDir, 'RotationAngle_%08i.log' % runNo),
                                   'python',
                                   'Scripts/ComputeRotationAngle.py', str(runNo)]),
                 runNo))

# Wait for all of the processes to be done; none should fail.
# Many won't produce output databases, for good reasons; we'll just omit them and interpolate the gain.
for p in proc:
    p[0].wait()
    if p[0].returncode != 0: print "Run %08i failed." % p[1]
print "Done."
