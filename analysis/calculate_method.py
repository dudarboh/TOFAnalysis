import ROOT
import sys
ROOT.gInterpreter.Declare('#include "extract.hpp"')

ch = ROOT.TChain("TOFAnalysis")

# # 2f_Z_hadronic data
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")


def calculate_method(momentum="pTrackAtCalo", length="lengthTrackCalo", tof="Closest", smearing="0."):
    '''
    Return root file with momentum, beta, PDG for every PFO in the barrel.
    '''
    df = ROOT.RDataFrame(ch)#.Filter('if(rdfentry_ % 1000000 == 0) {cout<<"Event: "<<rdfentry_<< endl;} return true;')
    df = df.Filter("nECALHits > 0 && abs(xyzCluster.Z()) < 2200.")

    df = df.Define("mom", "{}.R()".format(momentum))
    df = df.Define("tofHit", "tofHit(tECALHit, {})".format(smearing))
    if tof == "Closest":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("tof", "tofClosest(tofHit, dToImpact)")
    elif tof == "Fastest":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("tof", "tofFastest(tofHit, dToImpact)")
    elif tof == "Frank":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
               .Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
               .Define("tof", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 0)")
    elif tof == "Cyl":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
               .Define("sel_cyl", "selectHits(dToLine, layerECALHit, false, 10, 5.)")\
               .Define("tof", "fitFunc(tofHit[sel_cyl], dToImpact[sel_cyl], 0)")
    elif tof == "All":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("tof", "fitFunc(tofHit, dToImpact, 0)")
    elif tof == "FrankAvg":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
               .Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
               .Define("tof", "tofAvg(tofHit[sel_frank], dToImpact[sel_frank])")

    df = df.Define("beta", "{}/(tof * SPEED_OF_LIGHT)".format(length))

    method = "_".join([momentum, length, tof, smearing])
    print("Started run for", method)
    df.Snapshot(method, "./{}.root".format(method), ["mom", "beta", "PDG"])

# Example of arguments that can be passed to the comman line
# methods = [["pTrackAtIP", "lengthTrackIP", "FrankAvg", "10.0"], ["pTrackAtIP", "lengthTrackIP", "FrankAvg", "50.0"], ["pTrackAtIP", "lengthTrackIP", "FrankAvg", "100.0"],
#            ["pTrackAtCalo", "lengthTrackCalo", "Frank", "10.0"], ["pTrackAtCalo", "lengthTrackCalo", "Frank", "50.0"], ["pTrackAtCalo", "lengthTrackCalo", "Frank", "100.0"],]


if __name__ == "__main__":
    calculate_method(*sys.argv[1:])
