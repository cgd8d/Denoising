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
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libEXOUtilities")
ROOT.gSystem.Load("libEXOCalibUtilities")

# modules below this line may not actually be needed.
import sys
import math

# Create an empty sqlite3 database to hold event information.
conn = sqlite3.connect('Tmp/ThoriumLightmapEvents.db')
conn.execute('CREATE TABLE events (runRange integer, xpos real, ypos real, zpos real, ' +
                                   ', '.join(['apd_%03i_magnitude real' % apd for apd in Common.APDs]) + ')')

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
    NumRunsUsed += 1



####################
# BELOW THIS POINT, STILL IN PROGRESS.
####################


















for i in xrange(chain.GetEntries()):
    if i % 1000 == 0: print "Getting entry %i." % i
    chain.GetEntry(i)

    if event.fEventHeader.fTaggedAsNoise: continue

    for iscint in range(event.GetNumScintillationClusters()):
        scint = event.GetScintillationCluster(iscint)
        TotalScint += 1

        # If it's not SS, throw it away.
        if scint.GetNumChargeClusters() != 1: continue
        NumChargeClusters += 1
        charge = scint.GetChargeClusterAt(0)
        if len(set(charge.GetUWireSignalAt(j).fChannel for j in range(charge.GetNumUWireSignals()))) > 2:
            continue # If it has more than 2 distinct u-wire signals, it's not SS.
        NumUWires += 1

        # Use the best-available energy estimate to select full-energy events.
        if charge.fZ < 0:
            corr_factor = 0.938 + 0.6892*pow(abs(charge.fZ)/1000, 1.716)
        else:
            corr_factor = 0.9355 + 1.289*pow(abs(charge.fZ)/1000, 2.004)
        if UseDenoised:
            scintE = scint.fDenoisedEnergy/corr_factor
        else:
            scintE = scint.fRawEnergy
        AntiCorrE = charge.fPurityCorrectedEnergy*math.cos(Theta) + scintE*math.sin(Theta)
        Sigmas = (AntiCorrE - PeakPos)/(Res*PeakPos)
        if abs(Sigmas) > 1: continue
        EnergyCut += 1

        # Verify that all charge can be assigned unambiguously, at least with reference to this scintillation.
        charge_OK = True
        for icharge in range(event.GetNumChargeClusters()):
            chargeTime = event.GetChargeCluster(icharge).fCollectionTime
            if Utilities.CanPair(scint.fTime, chargeTime):
                # OK, this charge could be associated to this scint.
                # Can it be associated to any other scint?
                for jscint in range(event.GetNumScintillationClusters()):
                    if iscint == jscint: continue
                    if Utilities.CanPair(event.GetScintillationCluster(jscint).fTime, chargeTime):
                        charge_OK = False
                        break
            if not charge_OK: break
        if not charge_OK: continue
        UnambiguousCharge += 1

        # Verify that it is well-separated from other scintillation.
        # Here we are overly aggressive, but we want enough space to have a good baseline.
        well_separated = True
        for jscint in range(event.GetNumScintillationClusters()):
            if iscint == jscint: continue
            if scint.fTime > event.GetScintillationCluster(jscint).fTime:
                # Ensure that pretrace is protected.
                tdiff = scint.fTime - event.GetScintillationCluster(jscint).fTime
                if tdiff < (Utilities.ModelPretrace + Utilities.DiffTime)*ROOT.CLHEP.microsecond:
                    well_separated = False
                    break
            else:
                # Ensure that posttrace is protected.
                tdiff = event.GetScintillationCluster(jscint).fTime - scint.fTime
                if tdiff < Utilities.ModelPosttrace*ROOT.CLHEP.microsecond:
                    well_separated = False
                    break
        if not well_separated: continue
        ScintSeparation += 1

        # Basic good properties of the charge and scint clusters.
        if charge.fPurityCorrectedEnergy < 1: continue
        if charge.fPurityCorrectedEnergy > 5000: continue
        if abs(charge.fX) > 500 or abs(charge.fY) > 500: continue
        if not any(map(lambda z_range: charge.fZ >= z_range[0] and charge.fZ < z_range[1], FiducialZ_Ranges)):
            continue
        if UseDenoised:
            if scint.fEnergy < 1 or scint.fEnergy > 15000: continue
        else:
            if scint.fRawEnergy < 1 or scint.fRawEnergy > 15000: continue
        BasicChargeProperties += 1

        # Passes diagonal cut.
        diagCut = calibManager.getCalib("diagonal-cut", "vanilla", event.fEventHeader)
        if UseDenoised:
            if not diagCut.SurvivesSingleSiteCut(scint.fEnergy, charge.fPurityCorrectedEnergy): continue
        else:
            if not diagCut.SurvivesSingleSiteCut(scint.fRawEnergy, charge.fPurityCorrectedEnergy): continue
        DiagonalCut += 1

        # Compute the dot-product and integral of the raw waveforms with our model.
        SignalIndex = int(scint.fTime/ROOT.CLHEP.microsecond + 0.5)
        FirstIndex = SignalIndex - Utilities.ModelPretrace
        LastIndex = SignalIndex + Utilities.ModelPosttrace
        if FirstIndex < 0 or LastIndex > event.fEventHeader.fSampleCount + 1:
            NearEdge += 1
            continue
        WaveformMagnitudes = []
        for apd in Utilities.APDs:
            APDsignal = scint.GetAPDSignal(ROOT.EXOAPDSignal.kGangFit, apd)
            if APDsignal == None:
                amplitude = 0.
            else:
                scaling_factor = ROOT.APD_ADC_FULL_SCALE_ELECTRONS/(ROOT.ADC_BITS*ROOT.APD_GAIN)
                amplitude = APDsignal.fRawCounts / scaling_factor
            WaveformMagnitudes.append(amplitude)

        # Add it to our list; we want to update the database in one big operation.
        Entries.append( (Utilities.IndexOfRun(event.fRunNumber),
                         Utilities.GetBinAtPoint(charge.fX, charge.fY, charge.fZ)) +
                        tuple(WaveformMagnitudes))

# Print statistics for why we lost events.
print "TotalScint =", TotalScint
print "NumChargeClusters =", NumChargeClusters
print "NumUWires =", NumUWires
print "EnergyCut =", EnergyCut
print "UnambiguousCharge =", UnambiguousCharge
print "ScintSeparation =", ScintSeparation
print "NearEdge =", NearEdge
print "BasicChargeProperties =", BasicChargeProperties
print "DiagonalCut =", DiagonalCut

if len(Entries) < 200:
    print "Not enough accepted events in the run (%i)." % len(Entries)
    sys.exit("Not enough events accepted.")

dbname = os.path.join('PassingEvents', 'PassingEvents%s_%04i.db' % (segment_str, runNo))
conn = sqlite3.connect(dbname)
c = conn.cursor()
TableList = ('runRange integer, posBin integer,' +
             ', '.join(['apd_' + str(apd) + '_signal real' for apd in Utilities.APDs]))
c.execute('CREATE TABLE events (' + TableList + ')')
c.executemany('INSERT INTO events VALUES (' +
              ','.join(['?']*(2 + len(Utilities.APDs))) +
              ')', Entries)
conn.commit()
conn.close()

print "Done."
