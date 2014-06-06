
"""
For a given theta, find the resolution; then optimize resolution by varying theta.
This script is adapted by Clayton from a script by Manuel,
at http://java.freehep.org/svn/repos/exo/list/EXOCalibration/trunk/?revision=6540.
It is only intended to run on data, not MC.  (Different thresholds are not taken into account.)
It has its own rudimentary peak-finder, so it's intended to be a little robust against changes to the detector.
(That's not to say it can't fail.  Please report any unsound output it may produce.
 Note that in addition to the output values, it produces two pdf files with fits and so forth.
 This can help to pinpoint a problem.)

IMPORTANT ToDO: It would be very nice to unify this script with the one run during processing;
this would reduce the cost of generating a lightmap, reduce code duplication and the
difficulty of code maintenance, and possible let us incorporate improvements to the algorithm.
However, the standard script doesn't handle events with multiple scintillation clusters,
and we need to handle those to make use of strong source runs.  If someone can make the standard script
work with strong source runs, please get rid of this script and redirect the lightmap generation
code to use the standard pipeline output.
"""

import numpy
import ROOT
import math
import sys
ROOT.gROOT.SetBatch()

MaxRotatedEnergy = 10000.
UseDenoised = True

canvas = None
FileNames = { }

RunNumber = 0
FiducialZ_Ranges = [(-190, -5), (5, 190)]

class PeakFinderError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def IsFiducial(cluster):
    if not any(map(lambda z_range: cluster.fZ >= z_range[0]*ROOT.CLHEP.millimeter and
                                   cluster.fZ < z_range[1]*ROOT.CLHEP.millimeter, FiducialZ_Ranges)):
        return False # Not in any of the acceptable z-ranges.
    if pow(cluster.fX, 2) + pow(cluster.fY, 2) > pow(171*ROOT.CLHEP.millimeter, 2): return False
    return True

def IsSingleSite(scint):
    if scint.GetNumChargeClusters() != 1: return False
    cluster = scint.GetChargeClusterAt(0)
    UniqueU = set([cluster.GetUWireSignalAt(i).fChannel for i in range(cluster.GetNumUWireSignals())])
    if len(UniqueU) > 2: return False
    else: return True

def DivideWithErrors(x, y):
    # If x and y have uncorrelated errors, return their quotient (with errors).
    return (x[0]/y[0],
            math.sqrt(pow(x[1]/y[0], 2) +
                      pow(x[0]*y[1] / pow(y[0], 2), 2)))

def FindPeak(dataSet):
    # Given an unbinned RooDataSet, locate the fit range for the photopeak.
    # OK -- have to find peaks near the end.  First bin it coursely.
    hist = ROOT.TH1F('PeakFinder', 'PeakFinder', int(MaxRotatedEnergy/200), 0., MaxRotatedEnergy)
    if dataSet.fillHistogram(hist, ROOT.RooArgList(dataSet.get())) == None:
        raise PeakFinderError('Failed to bin data for peak-finding.')
    for i in range(hist.GetNbinsX()-1, 0, -1): # Step from right to left.
        # We want the chance of an accidental peak to be very small.
        # Here we demand it's more than 4 sigmas.  (3 would probably be sufficient, if 4 is too restrictive.)
        if hist.GetBinContent(i+1) <= 0: continue
        sigma = math.sqrt(hist.GetBinContent(i) + hist.GetBinContent(i+1))
        if (hist.GetBinContent(i+1) - hist.GetBinContent(i))/sigma >= 4:
            Peak = hist.GetBinCenter(i+1)
            break
    try:
        if i == 1: raise PeakFinderError('Failed to find a peak.')
        if Peak < 300: raise PeakFinderError('Peak found, but the energy is only ' + str(Peak) + '.')
        return Peak
    except PeakFinderError:
        canvas = ROOT.TCanvas()
        hist.Draw()
        canvas.Print('ComputeRotationAngle_oldversion/FailedToFindPeak%08i.pdf' % RunNumber)
        raise

def DoFit(dataSet, PeakPosGuess, SigmaGuess, ErfcFracGuess, ID = None):
    # Using input guesses, pick a fit window and refine those guesses.
    # If ID is given, use it to save a plot.
    RotatedEnergy = dataSet.get()['RotatedEnergy']

    # Note, the bounds used here imply resolution can never be worse than 30% at the peak.
    mean = ROOT.RooRealVar('Mean', 'Mean', PeakPosGuess, PeakPosGuess - SigmaGuess, PeakPosGuess + SigmaGuess)
    sigma = ROOT.RooRealVar('Sigma', 'Sigma', SigmaGuess, 0., 0.3*PeakPosGuess)
    erfc_coeff = ROOT.RooRealVar('coeff', 'Erfc coeff', ErfcFracGuess, 0., 1.)
    GaussianComponent = ROOT.RooGaussian('GaussianComponent', 'GaussianComponent',
                                         RotatedEnergy, mean, sigma)
    ErfcComponent = ROOT.RooGenericPdf('ErfcComponent', 'ErfcComponent',
                                       'TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))',
                                       ROOT.RooArgList(RotatedEnergy, mean, sigma))
    FitFunction = ROOT.RooAddPdf('FitFunction', 'FitFunction',
                                 GaussianComponent, ErfcComponent, erfc_coeff)
    FitFunction.fitTo(dataSet,
                      ROOT.RooFit.Range(PeakPosGuess - 4*SigmaGuess, PeakPosGuess + 3*SigmaGuess),
                      ROOT.RooFit.PrintLevel(-1))

    if ID:
        # ID[0] should identify the fit as single-site or multi-site.
        # ID[1] should be a formatted identifier of theta.

        # The plot should be wider than the fit window, so we can see what wasn't fit.
        frame = RotatedEnergy.frame(ROOT.RooFit.Bins(200),
                                    ROOT.RooFit.Range(PeakPosGuess - 6*SigmaGuess,
                                                      PeakPosGuess + 6*SigmaGuess))
        dataSet.plotOn(frame,
                       ROOT.RooFit.DrawOption('ZP'),
                       ROOT.RooFit.MarkerStyle(20),
                       ROOT.RooFit.MarkerSize(0.5))
        FitFunction.plotOn(frame)
        FitFunction.plotOn(frame,
                           ROOT.RooFit.Components("ErfcComponent"),
                           ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2))
        FitFunction.plotOn(frame,
                           ROOT.RooFit.Components("GaussianComponent"),
                           ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2))
        FitFunction.paramOn(frame,
                            ROOT.RooFit.Layout(0.6))
        frame.SetTitle(ID[0] + ' rotated spectrum (Theta = ' + ID[1] + ' rad)')
        frame.SetXTitle('Rotated Energy (uncalibrated)')
        frame.Draw()
        canvas.Print(FileNames[ID[0]])

    return ((mean.getVal(), mean.getError()),
            (sigma.getVal(), sigma.getError()),
            (erfc_coeff.getVal(), sigma.getError()))

def GetResolutionForTheta(ID, EnergyList2D, Theta):
    # Build rotated energy dataset for this theta
    RotatedEnergy = ROOT.RooRealVar('RotatedEnergy', 'RotatedEnergy', 0., MaxRotatedEnergy)
    ArgSet = ROOT.RooArgSet(RotatedEnergy)
    dataSet = ROOT.RooDataSet('dataSet', 'dataset', ArgSet)
    for Energy2D in EnergyList2D:
        RotatedEnergy.setVal(Energy2D[0]*math.cos(Theta) + Energy2D[1]*math.sin(Theta))
        dataSet.add(ArgSet)

    # Find the peak.
    InitPeakPos = FindPeak(dataSet)

    # Initialize guesses.  (Uncertainties are bogus.)
    PeakPos = (InitPeakPos, 0.1*InitPeakPos)
    Sigma = (0.02*InitPeakPos, 0.01*InitPeakPos)
    ErfcFrac = (0.1, 0.01)

    # Do a series of intermediate fits, to improve guesses.
    # This turns out to be necessary to get smooth results in some cases.
    for NumFit in range(3):
        PeakPos, Sigma, ErfcFrac = DoFit(dataSet, PeakPos[0], Sigma[0], ErfcFrac[0])

    # Final fit.  Save the results.
    PeakPos, Sigma, ErfcFrac = DoFit(dataSet,
                                     PeakPos[0],
                                     Sigma[0],
                                     ErfcFrac[0],
                                     ID = (ID, '%.5f' % Theta))

    Resolution = DivideWithErrors(Sigma, PeakPos)
    return Resolution, PeakPos

def Run(prefix, **kwargs):

    # Only process thorium source runs.
    beginRecord = kwargs['ControlRecordList'].GetNextRecord('EXOBeginRecord')()
    if not isinstance(beginRecord, ROOT.EXOBeginSourceCalibrationRunRecord):
        print "This is not a source run; skipping."
        return
    if 'Th' not in beginRecord.GetSourceTypeString():
        print "We only handle thorium runs here; skipping."
        return

    # Set up plotting 
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetOptFit(11)
    global canvas
    global FileNames
    FileNames['SS'] = prefix + '_ss.pdf'
    canvas = ROOT.TCanvas()
    for name in FileNames.itervalues(): canvas.Print(name + '[')

    # We'll need the drift velocity.
    ROOT.gSystem.Load("libEXOCalibUtilities")
    CalibManager = ROOT.EXOCalibManager.GetCalibManager()
    print "There are " + str(kwargs['EventTree'].GetEntries()) + " entries."
    EnergyList2D_ss = []
    EnergyList2D_ms = []
    # List to hold 2D energy points.
    event = ROOT.EXOEventData()
    kwargs['EventTree'].SetBranchAddress("EventBranch", event)

    # Do this once, and get the charge and light values after all cuts.
    print "Starting to read in data...",
    for i_event in xrange(kwargs['EventTree'].GetEntries()):
        kwargs['EventTree'].GetEntry(i_event)

        # Verify that the electron lifetime is sufficient -- say, more than 2ms.
        lifecalib = ROOT.EXOCalibManager.GetCalibManager().getCalib(
            "a-el-life", "TPC1_2013_0nu", event.fEventHeader.fTriggerSeconds)
        if lifecalib.lifetime(float(event.fEventHeader.fTriggerSeconds)) < 2000*ROOT.CLHEP.microsecond:
            print "Killing run on account of inadequate purity."
            sys.exit() # Kill the whole run if purity is insufficient at any point.


        numScint = int(event.GetNumScintillationClusters())
        if numScint == 0: continue
        ScintList = map(lambda i: (i, event.GetScintillationCluster(i)), range(numScint))
        ScintList.sort(key = lambda scint: scint[1].fTime)

        # Remove any scintillation clusters within a drift time of each other, or the end of the trace.
        DriftVelocity = CalibManager.getCalib("drift", "vanilla", event.fEventHeader)
        MaxDriftTime = ROOT.CATHODE_ANODE_x_DISTANCE*min(DriftVelocity.get_drift_velocity_TPC1(),
                                                         DriftVelocity.get_drift_velocity_TPC2())
        TempScintList = []
        for i in range(len(ScintList)):
            if ScintList[i][1].fTime < MaxDriftTime: continue
            if (event.fEventHeader.fSampleCount+1)*ROOT.SAMPLE_TIME - ScintList[i][1].fTime < MaxDriftTime: continue
            if i != 0 and ScintList[i][1].fTime - ScintList[i-1][1].fTime < MaxDriftTime: continue
            if i+1 < len(ScintList) and ScintList[i+1][1].fTime - ScintList[i][1].fTime < MaxDriftTime: continue
            TempScintList.append(ScintList[i])
        ScintList = TempScintList
        # ScintList has now been filtered -- it only contains scintillation clusters with unambiguous charge assignment.

        for scint in ScintList:
            if scint[1].GetNumChargeClusters() == 0: continue
            # If any charge cluster is non-fiducial, skip the whole scintillation cluster.
            if any(not IsFiducial(scint[1].GetChargeClusterAt(i))
                   for i in range(scint[1].GetNumChargeClusters())): continue
            # OK, if we got to this point we're keeping the scintillation cluster.
            ChargeEnergy = 0.
            for i_charge in range(scint[1].GetNumChargeClusters()):
                ChargeEnergy += scint[1].GetChargeClusterAt(i_charge).fPurityCorrectedEnergy
            if IsSingleSite(scint[1]):
                charge_z = scint[1].GetChargeClusterAt(0).fZ
                if charge_z < 0:
                    corr_factor = 0.938 + 0.6892*pow(abs(charge_z)/1000, 1.716)
                else:
                    corr_factor = 0.9355 + 1.289*pow(abs(charge_z)/1000, 2.004)
                if UseDenoised:
                    scintE = scint[1].fDenoisedEnergy/corr_factor
                else:
                    scintE = scint[1].fRawEnergy
                EnergyList2D_ss.append( (ChargeEnergy, scintE) )
            else:
                EnergyList2D_ms.append( (ChargeEnergy, scint[1].fRawEnergy) )
    print "Done."
    print len(EnergyList2D_ss),'passed the single-site cuts.'
    print len(EnergyList2D_ms),'passed the multi-site cuts.'
    if len(EnergyList2D_ss) < 500:
        print "Skipping run for too few events."
        return
    if len(EnergyList2D_ss) < 1000: print "ALERT: low-statistics run." # For future review.

    ThetaToTry = [0.3 + 0.01*i for i in range(40)]
    graph_ss = ROOT.TGraphErrors()

    BestTheta_ss = None
    BestRes_ss = None

    try:
        for TestTheta in ThetaToTry:
            ssResolution, _ = GetResolutionForTheta('SS', EnergyList2D_ss, TestTheta)

            graph_ss.Set(graph_ss.GetN()+1)
            graph_ss.SetPoint(graph_ss.GetN()-1, TestTheta, ssResolution[0])
            graph_ss.SetPointError(graph_ss.GetN()-1, 0., ssResolution[1])

            if BestRes_ss == None or ssResolution[0] < BestRes_ss:
                BestRes_ss = ssResolution[0]
                BestTheta_ss = TestTheta

        # ToDO:  Consider trying to pick a fit window more intelligently.
        ssPolynomialFunc = ROOT.TF1("ssFit", "[0] + [1]*(x-[2])**2", BestTheta_ss - 0.05, BestTheta_ss + 0.05)
        ssPolynomialFunc.SetParameters(BestRes_ss, 4.0, BestTheta_ss)
        graph_ss.Fit(ssPolynomialFunc, "ER")

        ROOT.gStyle.SetOptFit(11)
        graph_ss.SetTitle('SS Resolution vs Theta')
        graph_ss.GetXaxis().SetTitle('Theta (radians)')
        graph_ss.GetYaxis().SetTitle('Resolution (at thorium peak)')
        graph_ss.Draw("A*")
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs(FileNames['SS'])

        ssResolution, ssMean = GetResolutionForTheta('SS', EnergyList2D_ss, ssPolynomialFunc.GetParameter(2))

        for name in FileNames.itervalues(): canvas.Print(name + ']')

        return { 'Theta_ss' : (ssPolynomialFunc.GetParameter(2), ssPolynomialFunc.GetParError(2)),
                 'PeakPos_ss' : ssMean,
                 'Resolution_ss' : ssResolution }

    except PeakFinderError, err:
        # Clean up and terminate -- hopefully a person will look at this, and fix it.
        for name in FileNames.itervalues(): canvas.Print(name + ']')
        print 'Peak-finding failed: ', err
        print 'Terminating this script!!!  (Someone please figure out what went wrong.)'
        return

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'We need the run number to be passed in as an argument:'
        print 'python ComputeRotationAngle.py <run_number>'
        sys.exit()
    RunNumber = int(sys.argv[1])
    import ROOT
    ROOT.TH1.AddDirectory(False)
    ROOT.gROOT.SetBatch()
    if ROOT.gSystem.Load("libEXOUtilities") < 0: sys.exit('Failed to load EXOUtilities.')
    import glob
    import time
    import random
    globname = '/nfs/slac/g/exo_data3/exo_data/data/WIPP/DN_Source_LJPurity_Jan2014/' + str(int(sys.argv[1])) + '/denoised*.root'
    EventTree = ROOT.TChain('tree')
    numFilesAddedToTree = EventTree.Add(globname)
    AllMaskedFiles = glob.glob(globname)
    while len(AllMaskedFiles) != numFilesAddedToTree:
        print "Failed to get the number of masked files expected; try again."
        time.sleep(120*random.random())
        EventTree.Reset()
        numFilesAddedToTree = EventTree.Add(globname)
        AllMaskedFiles = glob.glob(globname)
    if len(AllMaskedFiles) == 0:
        # In case, for instance, we chose not to denoise a particular run.
        print "No files exist for this run; exiting."
        sys.exit()

    AllMaskedFiles.sort(reverse = True)
    LastFile = ROOT.TFile(AllMaskedFiles[0])
    LastTree = LastFile.Get('tree')
    ControlRecordList = LastTree.GetUserInfo().At(1)

    result = Run('Tmp/ComputeRotationAngle_oldversion/RotationAngle_%08i' % int(sys.argv[1]),
                 EventTree=EventTree, ControlRecordList=ControlRecordList)
    print result

