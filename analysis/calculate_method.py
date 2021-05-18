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
    df = ROOT.RDataFrame(ch)
    df = df.Filter("nECALHits > 0")

    df = df.Define("mom", "{}.R()".format(momentum))\
           .Define("pt", "{}.Rho()".format(momentum))
    df = df.Define("tofNoSmearingHit", "tofHit(tECALHit, 0.0)").Define("tofHit", "tofHit(tECALHit, {})".format(smearing))

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
               .Define("tof", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 0, 0)")
               # .Define("slope0", "fitFunc(tofNoSmearingHit[sel_frank], dToImpact[sel_frank], 1, 0)")\
               # .Define("tof100", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 0, 0)")\
               # .Define("slope100", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 1, 0)")\
               # .Define("tof100_lim", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 0, 1)")\
               # .Define("slope100_lim", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 1, 1)")
    elif tof == "Cyl":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
               .Define("sel_cyl", "selectHits(dToLine, layerECALHit, false, 10, 5.)")\
               # .Define("tof", "fitFunc(tofHit[sel_cyl], dToImpact[sel_cyl], tofNoSmearingHit, 0)")\
               # .Define("slope", "1./fitFunc(tofHit[sel_cyl], dToImpact[sel_cyl], tofNoSmearingHit, 1)")
               # .Define("tof", "fitFunc(tofHit[sel_cyl], dToImpact[sel_cyl], 0)")\
               # .Define("slope", "1./fitFunc(tofHit[sel_cyl], dToImpact[sel_cyl], 1)")
    elif tof == "All":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               # .Define("tof", "fitFunc(tofHit, dToImpact, tofNoSmearingHit, 0)")\
               # .Define("slope", "1./fitFunc(tofHit, dToImpact, tofNoSmearingHit, 1)")
               #
               # .Define("tof", "fitFunc(tofHit, dToImpact, 0)")\
               # .Define("slope", "1./fitFunc(tofHit, dToImpact, 1)")
    elif tof == "AllAvg":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("tof", "tofAvg(tofHit, dToImpact)")
    elif tof == "FrankAvg":
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
               .Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
               .Define("tof", "tofAvg(tofHit[sel_frank], dToImpact[sel_frank])")

    df = df.Define("beta", "{}/(tof * SPEED_OF_LIGHT)".format(length))
    # df = df.Filter("slope0 < 20.").Range(30).Define("kek_test", "fit_analysis(tofHit[sel_frank], dToImpact[sel_frank], tofNoSmearingHit[sel_frank], PDG, pTrackAtIP, nECALHits)")

    method = "_".join([momentum, length, tof, smearing])
    print("Started run for", method)
    # df.Snapshot(method, "./{}.root".format(method), ["mom", "beta", "PDG", "tof0", "slope0", "tof100", "slope100", "tof100_lim", "slope100_lim", "{}".format(length)])
    df.Snapshot(method, "./{}.root".format(method), ["mom", "beta", "PDG", "tof"])


# Example of arguments that can be passed to the comman line
# methods = [["pTrackAtIP", "lengthTrackIP", "FrankAvg", "10.0"], ["pTrackAtIP", "lengthTrackIP", "FrankAvg", "50.0"], ["pTrackAtIP", "lengthTrackIP", "FrankAvg", "100.0"],
#            ["pTrackAtCalo", "lengthTrackCalo", "Frank", "10.0"], ["pTrackAtCalo", "lengthTrackCalo", "Frank", "50.0"], ["pTrackAtCalo", "lengthTrackCalo", "Frank", "100.0"],]


if __name__ == "__main__":
    calculate_method(*sys.argv[1:])
