'''
Contain information which is common to multiple scripts.
'''

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libEXOUtilities")

# List of apd channels to use.
# Channels 178, 191, and 205 were never good.
# Some others have been disconnected during running, but all of these apds
# have been good at some point during physics running.
# For purposes of the lightmap, recently bad channels will be handled by
# forcing their gain to zero.
APDs = [ channel for channel in range(152, 226) if channel not in [178, 191, 205] ]

# Verify that all trees in the chain, for all processing stages,
# had all calibration database quantities available.
def VerifyCalibrations(chain):
    event = ROOT.EXOEventData()
    chain.SetBranchAddress("EventBranch", event)
    NextEntryToLoad = 0
    for i in range(chain.GetNtrees()):
        chain.LoadTree(NextEntryToLoad)
        tree = chain.GetTree()
        ProcessingInfo = tree.GetUserInfo().At(0)
        while ProcessingInfo:
            if not ProcessingInfo.GetCalibrationsFromDatabase(): return False
            ProcessingInfo = ProcessingInfo.GetPrevProcInfo()
        if tree.GetEntries() <= 0: sys.exit('A tree has non-positive entries; how?')
        NextEntryToLoad += tree.GetEntries()
    return True
