"""
Call this module like: wrap_python.sh GrabFullEnergySSTh_SingleRun.py <run_number>
It uses the results from ComputeRotationAngle.py, and other criteria, to identify
full-energy SS Th events.
It then produces a sqlite database of such events from this run.
"""

import Common
import os
import sqlite3
import re
import math
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libEXOUtilities")
ROOT.gSystem.Load("libEXOCalibUtilities")

# Create an empty sqlite3 database to hold event information.
conn = sqlite3.connect('Tmp/ThoriumLightmapEvents.db')
conn.execute('CREATE TABLE events (runNo integer, xpos real, ypos real, zpos real, ' +
                                   ', '.join(['apd_%03i_magnitude real' % apd for apd in Common.APDs]) + ')')
InsertStmt = 'INSERT INTO events VALUES (%s)' % ','.join(['?']*(4 + len(Common.APDs)))

# Make a list of ComputeRotationAngle log files.
LogFileDir = 'Tmp/ComputeRotationAngle_oldversion'
ListOfLogFiles = os.listdir(LogFileDir)

# Make a regular expression which verifies that the log file name matches the expected pattern
# and extracts the run number when it does.
file_regexp = re.compile('^RotationAngle_(\d{8})\.log$')

# Also precompile regular expressions for the terms we'll search for within the log files.
FloatRegExp = '[\-\+]?\d+(?:\.\d+)?(?:e[\-\+]?\d+)?' # Regular expression for a floating-point number (!).
theta_regexp = re.compile('\'Theta_ss\': \((' + FloatRegExp + ')')
peakpos_regexp = re.compile('\'PeakPos_ss\': \((' + FloatRegExp + ')')
resol_regexp = re.compile('\'Resolution_ss\': \((' + FloatRegExp + ')')

# Get the calib manager.
calibManager = ROOT.EXOCalibManager.GetCalibManager()

# Keep some statistics on cuts applied.
NumRunsTotal = 0
NumRunsCut_NoCalib = 0
NumRunsCut_RotAngleFailed = 0
NumRunsCut_Res = 0
NumRunsCut_PeakPos = 0
NumRunsUsed = 0
NumEntriesTotal = 0
NumEntriesCut_Noise = 0
NumScintTotal = 0
NumScintCut_MS = 0
NumScintCut_NumWires = 0
NumScintCut_Energy = 0
NumScintCut_AmbiguousCharge = 0
NumScintCut_Separation = 0
NumScintCut_Suspicious = 0
NumScintCut_Diag = 0
NumScintCut_NearEdge = 0
NumScintUsed = 0

# Now iterate through the log files and add entries to the database as available.
for logfilename in ListOfLogFiles:
    filename_match = file_regexp.match(logfilename)
    if filename_match == None: continue
    runNo = int(filename_match.group(1))

    # Get the denoised data chain for this run, and verify it's OK.
    chain = ROOT.TChain("tree")
    num_files = chain.Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/DN_Source_LJPurity_Jan2014/%i/denoised%08i-*.root" % (runNo, runNo))
    if num_files == 0:
        print "Skipping run %i because no files could be found." % runNo
        continue
    if chain.GetEntries() == 0:
        print "Skipping run %i because it is empty." % runNo
        continue
    chain.LoadTree(0)
    beginRec = chain.GetTree().GetUserInfo().At(1).GetNextRecord('EXOBeginRecord')()
    if not isinstance(beginRec, ROOT.EXOBeginSourceCalibrationRunRecord):
        print "Skipping run %i because it has no begin-run record." % runNo
        continue
    if (beginRec.GetSourceType() != ROOT.EXOBeginSourceCalibrationRunRecord.kThWeak and
        beginRec.GetSourceType() != ROOT.EXOBeginSourceCalibrationRunRecord.kThStrong):
        print "Skipping run %i because it is not a thorium run." % runNo
        continue

    # Runs which make it to this point are, in principle, permissible.
    # Collect statistics on accepted/cut runs starting here.
    NumRunsTotal += 1

    if not Common.VerifyCalibrations(chain):
        print "Skipping run %i because some calibrations weren't available." % runNo
        NumRunsCut_NoCalib += 1
        continue
    event = ROOT.EXOEventData()
    chain.SetBranchAddress("EventBranch", event)

    # Extract relevant information from the log file.
    logfilecontents = open(logfilename).read()
    theta_match = theta_regexp.search(logfilecontents)
    peakpos_match = peakpos_regexp.search(logfilecontents)
    resol_match = resol_regexp.search(logfilecontents)
    if theta_match == None or peakpos_match == None or resol_match == None:
        print "Skipping run %i because its ComputeRotationFile script didn't finish properly." % runNo
        NumRunsCut_RotAngleFailed += 1
        continue
    theta = float(theta_match.group(1))
    peakpos = float(peakpos_match.group(1))
    resol = float(resol_match.group(1))

    # Apply some cuts on whether the run is usable.
    if resol < 0.01 or resol > 0.03:
        # This is a pretty effective cut on bad peak fits.
        print "Skipping run %i because it's resolution was %f." % (runNo, resol)
        NumRunsCut_Res += 1
        continue
    if peakpos < 3000:
        # This cut has helped to reject Cs falsely identified as Th.
        print "Skipping run %i because peakpos was %f." % (runNo, peakpos)
        NumRunsCut_PeakPos += 1
        continue

    # If we've made it to this point, the run will be used.
    print "Using run %i." % runNo
    NumRunsUsed += 1

    # Get a cursor for the database.
    cursor = conn.cursor()

    for i in xrange(chain.GetEntries()):
        chain.GetEntry(i)
        NumEntriesTotal += 1

        if event.fEventHeader.fTaggedAsNoise:
            NumEntriesCut_Noise += 1
            continue

        for iscint in range(event.GetNumScintillationClusters()):
            scint = event.GetScintillationCluster(iscint)
            NumScintClustersTotal += 1

            # If it's not SS, throw it away.
            if scint.GetNumChargeClusters() != 1:
                NumScintCut_MS += 1
                continue
            charge = scint.GetChargeClusterAt(0)

            # If it has u-wire signals on more than two distinct channels, it's not SS.
            if len(set(charge.GetUWireSignalAt(j).fChannel for j in range(charge.GetNumUWireSignals()))) > 2:
                NumScintCut_NumWires += 1
                continue

            # Compute anticorrelated energy.
            if charge.fZ < 0:
                corr_factor = 0.938 + 0.6892*pow(abs(charge.fZ)/1000, 1.716)
            else:
                corr_factor = 0.9355 + 1.289*pow(abs(charge.fZ)/1000, 2.004)
            scintE = scint.fDenoisedEnergy/corr_factor
            AntiCorrE = charge.fPurityCorrectedEnergy*math.cos(theta) + scintE*math.sin(theta)

            # Accept events within 1.5sigma of the peak.
            Sigmas = (AntiCorrE - peakpos)/(resol*peakpos)
            if abs(Sigmas) > 1.5:
                NumScintCut_Energy += 1
                continue

            # Verify that all charge can be assigned unambiguously.
            # In other words, there are no charge clusters which *could* come from this scint
            # and this charge cluster can't belong to any other scint clusters.
            charge_OK = True
            for jscint in range(event.GetNumScintillationClusters()):
                if jscint == iscint: continue
                drift_time = charge.fCollectionTime - event.GetScintillationCluster(jscint).fTime
                if drift_time > -5*ROOT.CLHEP.microsecond and drift_time < 120*ROOT.CLHEP.microsecond:
                    charge_OK = False
                    break
            if not charge_OK:
                NumScintCut_AmbiguousCharge += 1
                continue
            for icharge in range(event.GetNumChargeClusters()):
                if event.GetChargeCluster(icharge) == charge: continue
                chargeTime = event.GetChargeCluster(icharge).fCollectionTime
                drift_time = chargeTime - scint.fTime
                if drift_time > -5*ROOT.CLHEP.microsecond and drift_time < 120*ROOT.CLHEP.microsecond:
                    charge_OK = False
                    break
            if not charge_OK:
                NumScintCut_AmbiguousCharge += 1
                continue

            # Verify that it is well-separated from other scintillation.
            # This is to ensure a good fit for the APD signals (though maybe not really needed?).
            well_separated = True
            for jscint in range(event.GetNumScintillationClusters()):
                if iscint == jscint: continue
                diff_time = scint.fTime - event.GetScintillationCluster(jscint).fTime
                if abs(diff_time) < 150*ROOT.CLHEP.microsecond:
                    well_separated = False
                    break
            if not well_separated:
                 NumScintCut_Separation += 1
                 continue

            # Verify basic good properties of the charge.
            # (Should add a check on the light range as well.)
            if (charge.fPurityCorrectedEnergy < 1 or
                charge.fPurityCorrectedEnergy > 5000 or
                abs(charge.fX) > 500 or abs(charge.fY) > 500 or abs(charge.fZ) > 500):
                NumScintCut_Suspicious += 1
                continue

            # Passes diagonal cut.
            diagCut = calibManager.getCalib("diagonal-cut", "2013-0nu-denoised", event.fEventHeader)
            if not diagCut.SurvivesSingleSiteCut(scint.fEnergy, charge.fPurityCorrectedEnergy): continue
                NumScintCut_Diag += 1
                continue

            # Check that the time of the scintillation is not near the edge of the waveform.
            if (scint.fTime < 150*ROOT.CLHEP.microsecond or
                scint.fTime > (event.fEventHeader.fSampleCount-150)*ROOT.CLHEP.microsecond):
                NumScintCut_NearEdge += 1
                continue

            # OK, we will use this scint cluster.
            NumScintUsed += 1

            # Prepare an sql row, and insert it.
            # Note that we use the fit magnitude from reconstruction, which is undenoised;
            # denoised individual-APD magnitudes could improve the quality of the lightmap.
            RowToInsert = [runNo, cluster.fX, cluster.fY, cluster.fZ]
            scaling_factor = ROOT.APD_ADC_FULL_SCALE_ELECTRONS/(ROOT.ADC_BITS*ROOT.APD_GAIN)
            for apd in Common.APDs:
                APDsignal = scint.GetAPDSignal(ROOT.EXOAPDSignal.kGangFit, apd)
                if APDsignal == None: RowToInsert.append(0.)
                else: RowToInsert.append(APDsignal.fRawCounts/scaling_factor)
            cursor.execute(InsertStmt, RowToInsert)

    # At the end of each run, commit the transaction.
    conn.commit()

# Build indices.
print "Building indices."
conn.execute('CREATE INDEX x_index ON events (xpos)')
conn.execute('CREATE INDEX y_index ON events (ypos)')
conn.execute('CREATE INDEX z_index ON events (zpos)')
conn.execute('CREATE INDEX run_index ON events (runNo)')
conn.execute('ANALYZE events') # Helps sqlite pick which index to use intelligently.
# Note: we could be more intelligent if we used a version of sqlite3 compiled
# with the SQLITE_ENABLE_STAT4 macro.  But this is awkward to do from the python interface.
print "Done building indices."

conn.close()

# Print run statistics.
print
print "STATISTICS:"
print "There were %i runs which were potentially usable Th runs." % NumRunsTotal
print "\t%i cut due to missing calibrations." % NumRunsCut_NoCalib
print "\t%i cut because the rotation angle script did not finish." % NumRunsCut_RotAngleFailed
print "\t%i cut because the resolution was unreasonable." % NumRunsCut_Res
print "\t%i cut because the peak position was unreasonable." % NumRunsCut_PeakPos
print "%i runs were actually used." % NumRunsUsed
print
print "These runs contained % events (total)." % NumEntriesTotal
print "\t%i were cut due to being tagged as noise." % NumEntriesCut_Noise
print "The acceptable events contained %i scint clusters (before additional cuts)." % NumScintTotal
print "\t%i were cut as multi-site (multiple charge clusters)." % NumScintCut_MS
print "\t%i were cut for hitting too many u-wires." % NumScintCut_NumWires
print "\t%i were cut for not being near enough to 2615 keV." % NumScintCut_Energy
print "\t%i were cut because charge/scint clustering was ambiguous." % NumScintCut_AmbiguousCharge
print "\t%i were cut because another scint cluster was too close." % NumScintCut_Separation
print "\t%i were cut because the amount of charge looked suspicious." % NumScintCut_Suspicious
print "\t%i were cut because they failed the diagonal cut." % NumScintCut_Diag
print "\t%i were cut because they were too close to the waveform edge." % NumScintCut_NearEdge
print "%i scintillation clusters were deemed usable for the lightmap." % NumScintUsed

