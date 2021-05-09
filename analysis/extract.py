import ROOT
import numpy as np

ROOT.gStyle.SetPalette(1)

# # 2f_Z_hadronic data
ch = ROOT.TChain("TOFAnalysis")
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")
df = ROOT.RDataFrame(ch)
ROOT.gInterpreter.Declare('#include "extract.hpp"')
ROOT.gInterpreter.Declare('int eventCounter = 0;')

def get_photons(df):
    # Check only good photons in barrel
    df = df.Filter("PDG == 22 && nECALHits > 0 && abs(xyzCluster.Z()) < 2200. && xyzVtxMC.R() < 0.5").Range(1000000).Filter('eventCounter++; if(eventCounter % 10000 == 0) {cout << eventCounter << endl;} return true; ')

    # All dependency calculations
    df = df.Define("ECALPlane", "getECALPlane( pMC.Phi() )")\
           .Define("rImpact", "intersection(xyzVtxMC, pMC, ECALPlane)")\
           .Define("dToImpact", "dToImpact(xyzECALHit, rImpact)")\
           .Define("dToLine", "dToLine(xyzECALHit, rImpact, pMC)")\
           .Define("mom", "pMC.R()")
    # smearings of ECAL hits times
    df = df.Define("tofHit0", "tofHit(tECALHit, 0)")\
           .Define("tofHit10", "tofHit(tECALHit, 10)")
           # .Define("tofHit30", "tofHit(tECALHit, 30)")\
           # .Define("tofHit50", "tofHit(tECALHit, 50)")\
           # .Define("tofHit100", "tofHit(tECALHit, 100)")\
           # .Define("tofHit200", "tofHit(tECALHit, 200)")\
           # .Define("tofHit300", "tofHit(tECALHit, 300)")

    #selection of hits
    df = df.Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
           .Define("sel_cyl", "selectHits(dToLine, layerECALHit, false, 10, 5.)")

    #tofs
    df = df.Define("tofTrue", "(rImpact - xyzVtxMC).R() / SPEED_OF_LIGHT")\
           .Define("tofClosest0", "tofClosest(tofHit0, dToImpact)")\
           .Define("tofFastest0", "tofFastest(tofHit0, dToImpact)")\
           .Define("tofFrank0", "fitFunc(tofHit0[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl0", "fitFunc(tofHit0[sel_cyl], dToImpact[sel_cyl], 0)")\
           .Define("tofAvg0", "tofAvg(tofHit0[sel_frank], dToImpact[sel_frank])")\

    df = df.Define("tofClosest10", "tofClosest(tofHit10, dToImpact)")\
           .Define("tofFastest10", "tofFastest(tofHit10, dToImpact)")\
           .Define("tofFrank10", "fitFunc(tofHit10[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl10", "fitFunc(tofHit10[sel_cyl], dToImpact[sel_cyl], 0)")\
           .Define("tofAvg10", "tofAvg(tofHit10[sel_frank], dToImpact[sel_frank])")

    # df = df.Define("tofClosest30", "tofClosest(tofHit30, dToImpact)")\
    #        .Define("tofFastest30", "tofFastest(tofHit30, dToImpact)")\
    #        .Define("tofFrank30", "fitFunc(tofHit30[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl30", "fitFunc(tofHit30[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest50", "tofClosest(tofHit50, dToImpact)")\
    #        .Define("tofFastest50", "tofFastest(tofHit50, dToImpact)")\
    #        .Define("tofFrank50", "fitFunc(tofHit50[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl50", "fitFunc(tofHit50[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest100", "tofClosest(tofHit100, dToImpact)")\
    #        .Define("tofFastest100", "tofFastest(tofHit100, dToImpact)")\
    #        .Define("tofFrank100", "fitFunc(tofHit100[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl100", "fitFunc(tofHit100[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest200", "tofClosest(tofHit200, dToImpact)")\
    #        .Define("tofFastest200", "tofFastest(tofHit200, dToImpact)")\
    #        .Define("tofFrank200", "fitFunc(tofHit200[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl200", "fitFunc(tofHit200[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest300", "tofClosest(tofHit300, dToImpact)")\
    #        .Define("tofFastest300", "tofFastest(tofHit300, dToImpact)")\
    #        .Define("tofFrank300", "fitFunc(tofHit300[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl300", "fitFunc(tofHit300[sel_cyl], dToImpact[sel_cyl], 0)")\

    df.Snapshot("photons", "/nfs/dust/ilc/user/dudarboh/analysis/photons.root", ["mom", "nECALHits", "tofTrue", "tofClosest0" ,"tofFastest0" ,"tofFrank0", "tofAvg0", "tofClosest10" ,"tofFastest10" ,"tofFrank10", "tofAvg10"])


def get_pions(df):
    # Check only good photons in barrel
    df = df.Filter("abs(PDG) == 211 && nECALHits > 0 && abs(xyzCluster.Z()) < 2200. && xyzVtxMC.R() < 0.5").Range(1000000)

    # All dependency calculations
    df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
           .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
           .Define("mom", "pTrackAtCalo.R()")\
    # smearings of ECAL hits times
    df = df.Define("tofHit0", "tofHit(tECALHit, 0)")\
           .Define("tofHit10", "tofHit(tECALHit, 10)")\
           .Define("tofHit30", "tofHit(tECALHit, 30)")\
           .Define("tofHit50", "tofHit(tECALHit, 50)")\
           .Define("tofHit100", "tofHit(tECALHit, 100)")\
           .Define("tofHit200", "tofHit(tECALHit, 200)")\
           .Define("tofHit300", "tofHit(tECALHit, 300)")

    #selection of hits
    df = df.Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
           .Define("sel_cyl", "selectHits(dToLine, layerECALHit, false, 10, 5.)")

    #tofs
    df = df.Define("tofClosest0", "tofClosest(tofHit0, dToImpact)")\
           .Define("tofFastest0", "tofFastest(tofHit0, dToImpact)")\
           .Define("tofFrank0", "fitFunc(tofHit0[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl0", "fitFunc(tofHit0[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("tofClosest10", "tofClosest(tofHit10, dToImpact)")\
           .Define("tofFastest10", "tofFastest(tofHit10, dToImpact)")\
           .Define("tofFrank10", "fitFunc(tofHit10[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl10", "fitFunc(tofHit10[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("tofClosest30", "tofClosest(tofHit30, dToImpact)")\
           .Define("tofFastest30", "tofFastest(tofHit30, dToImpact)")\
           .Define("tofFrank30", "fitFunc(tofHit30[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl30", "fitFunc(tofHit30[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("tofClosest50", "tofClosest(tofHit50, dToImpact)")\
           .Define("tofFastest50", "tofFastest(tofHit50, dToImpact)")\
           .Define("tofFrank50", "fitFunc(tofHit50[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl50", "fitFunc(tofHit50[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("tofClosest100", "tofClosest(tofHit100, dToImpact)")\
           .Define("tofFastest100", "tofFastest(tofHit100, dToImpact)")\
           .Define("tofFrank100", "fitFunc(tofHit100[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl100", "fitFunc(tofHit100[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("tofClosest200", "tofClosest(tofHit200, dToImpact)")\
           .Define("tofFastest200", "tofFastest(tofHit200, dToImpact)")\
           .Define("tofFrank200", "fitFunc(tofHit200[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl200", "fitFunc(tofHit200[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("tofClosest300", "tofClosest(tofHit300, dToImpact)")\
           .Define("tofFastest300", "tofFastest(tofHit300, dToImpact)")\
           .Define("tofFrank300", "fitFunc(tofHit300[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl300", "fitFunc(tofHit300[sel_cyl], dToImpact[sel_cyl], 0)")\

    df.Snapshot("pions", "/nfs/dust/ilc/user/dudarboh/analysis/pions.root")


def beta_vs_p(df):

    df = df.Filter("cal.nHits > 0").Range(10000000)

    # All dependency calculations
    df = df.Define("nHitsCluster", "int(xyzECALHit.size())")\
           .Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
           .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
           .Define("mom", "pTrackAtCalo.R()")\
           .Define("lengthIP", "abs((piFit.phiIP-piFit.phiCalo)/piFit.omegaIP)*sqrt(1. + piFit.tanLIP*piFit.tanLIP)")\
           .Define("lengthCalo", "abs((piFit.phiIP-piFit.phiCalo)/piFit.omegaCalo)*sqrt(1. + piFit.tanLCalo*piFit.tanLCalo)")\
    # smearings of ECAL hits times
    df = df.Define("tofHit0", "tofHit(xyzECALHit, 0)")\
           # .Define("tofHit10", "tofHit(xyzECALHit, 10)")\
           # .Define("tofHit30", "tofHit(xyzECALHit, 30)")\
           # .Define("tofHit50", "tofHit(xyzECALHit, 50)")\
           # .Define("tofHit100", "tofHit(xyzECALHit, 100)")\
           # .Define("tofHit200", "tofHit(xyzECALHit, 200)")\
           # .Define("tofHit300", "tofHit(xyzECALHit, 300)")

    #selection of hits
    df = df.Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
           .Define("sel_cyl", "selectHits(dToLine, layerECALHit, false, 10, 5.)")

    #tofs
    df = df.Define("tofClosest0", "tofClosest(tofHit0, dToImpact)")\
           .Define("tofFastest0", "tofFastest(tofHit0, dToImpact)")\
           .Define("tofFrank0", "fitFunc(tofHit0[sel_frank], dToImpact[sel_frank], 0)")\
           .Define("tofCyl0", "fitFunc(tofHit0[sel_cyl], dToImpact[sel_cyl], 0)")\

    # df = df.Define("tofClosest10", "tofClosest(tofHit10, dToImpact)")\
    #        .Define("tofFastest10", "tofFastest(tofHit10, dToImpact)")\
    #        .Define("tofFrank10", "fitFunc(tofHit10[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl10", "fitFunc(tofHit10[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest30", "tofClosest(tofHit30, dToImpact)")\
    #        .Define("tofFastest30", "tofFastest(tofHit30, dToImpact)")\
    #        .Define("tofFrank30", "fitFunc(tofHit30[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl30", "fitFunc(tofHit30[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest50", "tofClosest(tofHit50, dToImpact)")\
    #        .Define("tofFastest50", "tofFastest(tofHit50, dToImpact)")\
    #        .Define("tofFrank50", "fitFunc(tofHit50[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl50", "fitFunc(tofHit50[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest100", "tofClosest(tofHit100, dToImpact)")\
    #        .Define("tofFastest100", "tofFastest(tofHit100, dToImpact)")\
    #        .Define("tofFrank100", "fitFunc(tofHit100[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl100", "fitFunc(tofHit100[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest200", "tofClosest(tofHit200, dToImpact)")\
    #        .Define("tofFastest200", "tofFastest(tofHit200, dToImpact)")\
    #        .Define("tofFrank200", "fitFunc(tofHit200[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl200", "fitFunc(tofHit200[sel_cyl], dToImpact[sel_cyl], 0)")\
    #
    # df = df.Define("tofClosest300", "tofClosest(tofHit300, dToImpact)")\
    #        .Define("tofFastest300", "tofFastest(tofHit300, dToImpact)")\
    #        .Define("tofFrank300", "fitFunc(tofHit300[sel_frank], dToImpact[sel_frank], 0)")\
    #        .Define("tofCyl300", "fitFunc(tofHit300[sel_cyl], dToImpact[sel_cyl], 0)")\

    df = df.Define("beta", "lengthCalo/(tofCyl0 * SPEED_OF_LIGHT)")

    nBinsX = 100
    nBinsY = 200
    # log10(min/max momentum / GeV)
    minBinX = 0
    maxBinX = 1
    # for TOF beta = v/c values from 0 to 1
    minBinY = 0.7
    maxBinY = 1.05
    histbinsX = []
    histbinsY = []
    for i in range(nBinsX +1):
        histbinsX.append( pow( 10, minBinX + (maxBinX-minBinX)*i/nBinsX ) )
    for i in range(nBinsY +1):
        histbinsY.append( minBinY+ (maxBinY-minBinY)*i/nBinsY )


    h = df.Histo2D(("name", "title", nBinsX, np.array(histbinsX), nBinsY, np.array(histbinsY) ), "mom", "beta")
    h.Draw("colz")
    input("wait")


get_photons(df)
# get_pions(df)
# beta_vs_p(df)
