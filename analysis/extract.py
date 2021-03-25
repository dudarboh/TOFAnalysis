import ROOT
import numpy as np

ROOT.gStyle.SetPalette(1)

# # 2f_Z_hadronic data
ch = ROOT.TChain("PFO")
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch1 = ROOT.TChain("Cluster")
ch1.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch2 = ROOT.TChain("PionTrack")
ch2.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch3 = ROOT.TChain("TrackerHits")
ch3.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch4 = ROOT.TChain("ECALHits")
ch4.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch.AddFriend(ch1, "cluster")
ch.AddFriend(ch2, "piFit")
ch.AddFriend(ch3, "tr")
ch.AddFriend(ch4, "cal")

df = ROOT.RDataFrame(ch)
ROOT.gInterpreter.Declare('#include "extract.hpp"')

def get_photons(df):
    # Check only good photons in barrel
    df = df.Filter("nMC == 1 && PDG[0] == 22 && cal.nHits > 0 \
            && isBackscatter[0] == 0 && isDecayedInTracker[0] == 0\
            && abs(cluster.posCluster.Z()) < 2200. && vtxMC[0].R() < 0.5").Range(1000000)

    # All dependency calculations
    df = df.Define("nHitsCluster", "int(cal.posECALHit.size())")\
           .Define("ECALPlane", "getECALPlane( pMC[0].Phi() )")\
           .Define("rImpact", "intersection(vtxMC[0].Vect(), pMC[0].Vect(), ECALPlane)")\
           .Define("dToImpact", "dToImpact(cal.posECALHit, rImpact)")\
           .Define("dToLine", "dToLine(cal.posECALHit, rImpact, pMC[0].Vect())")\
           .Define("mom", "pMC[0].Vect().R()")
    # smearings of ECAL hits times
    df = df.Define("tofHit0", "tofHit(cal.posECALHit, 0)")\
           .Define("tofHit10", "tofHit(cal.posECALHit, 10)")\
           .Define("tofHit30", "tofHit(cal.posECALHit, 30)")\
           .Define("tofHit50", "tofHit(cal.posECALHit, 50)")\
           .Define("tofHit100", "tofHit(cal.posECALHit, 100)")\
           .Define("tofHit200", "tofHit(cal.posECALHit, 200)")\
           .Define("tofHit300", "tofHit(cal.posECALHit, 300)")

    #selection of hits
    df = df.Define("sel_frank", "selectHits(dToLine, cal.layer, true, 10, 9999.)")\
           .Define("sel_cyl", "selectHits(dToLine, cal.layer, false, 10, 5.)")

    #tofs
    df = df.Define("tofTrue", "(rImpact - vtxMC[0].Vect()).R() / SPEED_OF_LIGHT")\
           .Define("tofClosest0", "tofClosest(tofHit0, dToImpact)")\
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

    df.Snapshot("photons", "/nfs/dust/ilc/user/dudarboh/analysis/photons.root")


def get_pions(df):
    # Check only good photons in barrel
    df = df.Filter("nMC == 1 && abs(PDG[0]) == 211 && cal.nHits > 0 \
            && isBackscatter[0] == 0\
            && abs(cluster.posCluster.Z()) < 2200. && vtxMC[0].R() < 0.5").Range(1000000)

    # All dependency calculations
    df = df.Define("nHitsCluster", "int(cal.posECALHit.size())")\
           .Define("dToImpact", "dToImpact(cal.posECALHit, XYZVector(piFit.refCalo))")\
           .Define("dToLine", "dToLine(cal.posECALHit, XYZVector(piFit.refCalo), piFit.pCalo)")\
           .Define("mom", "piFit.pCalo.R()")\
           .Define("lengthIP", "abs((piFit.phiIP-piFit.phiCalo)/piFit.omegaIP)*sqrt(1. + piFit.tanLIP*piFit.tanLIP)")\
           .Define("lengthCalo", "abs((piFit.phiIP-piFit.phiCalo)/piFit.omegaCalo)*sqrt(1. + piFit.tanLCalo*piFit.tanLCalo)")\
    # smearings of ECAL hits times
    df = df.Define("tofHit0", "tofHit(cal.posECALHit, 0)")\
           .Define("tofHit10", "tofHit(cal.posECALHit, 10)")\
           .Define("tofHit30", "tofHit(cal.posECALHit, 30)")\
           .Define("tofHit50", "tofHit(cal.posECALHit, 50)")\
           .Define("tofHit100", "tofHit(cal.posECALHit, 100)")\
           .Define("tofHit200", "tofHit(cal.posECALHit, 200)")\
           .Define("tofHit300", "tofHit(cal.posECALHit, 300)")

    #selection of hits
    df = df.Define("sel_frank", "selectHits(dToLine, cal.layer, true, 10, 9999.)")\
           .Define("sel_cyl", "selectHits(dToLine, cal.layer, false, 10, 5.)")

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

    # .Define("betaCalo", "lengthCalo/(tof*SPEED_OF_LIGHT)")\
    # .Define("betaAvg", "lengthAvg/(tof*SPEED_OF_LIGHT)")\
    # .Define("beta75", "length75/(tof*SPEED_OF_LIGHT)")\
    # .Define("massCalo", "pReco.Vect().R() / betaCalo * sqrt(1. - betaCalo * betaCalo) * 1000")\
    # .Define("massAvg", "pReco.Vect().R() / betaAvg * sqrt(1. - betaAvg * betaAvg) * 1000")\
    # .Define("mass75", "pReco.Vect().R() / beta75 * sqrt(1. - beta75 * beta75) * 1000")\


def beta_vs_p(df):
    # df = df.Filter("nMC == 1 && cal.nHits > 0 \
    #         && isBackscatter[0] == 0\
    #         && abs(cluster.posCluster.Z()) < 2200. && vtxMC[0].R() < 0.5").Range(1000000)

    df = df.Filter("cal.nHits > 0").Range(10000000)

    # All dependency calculations
    df = df.Define("nHitsCluster", "int(cal.posECALHit.size())")\
           .Define("dToImpact", "dToImpact(cal.posECALHit, XYZVector(piFit.refCalo))")\
           .Define("dToLine", "dToLine(cal.posECALHit, XYZVector(piFit.refCalo), piFit.pCalo)")\
           .Define("mom", "piFit.pCalo.R()")\
           .Define("lengthIP", "abs((piFit.phiIP-piFit.phiCalo)/piFit.omegaIP)*sqrt(1. + piFit.tanLIP*piFit.tanLIP)")\
           .Define("lengthCalo", "abs((piFit.phiIP-piFit.phiCalo)/piFit.omegaCalo)*sqrt(1. + piFit.tanLCalo*piFit.tanLCalo)")\
    # smearings of ECAL hits times
    df = df.Define("tofHit0", "tofHit(cal.posECALHit, 0)")\
           # .Define("tofHit10", "tofHit(cal.posECALHit, 10)")\
           # .Define("tofHit30", "tofHit(cal.posECALHit, 30)")\
           # .Define("tofHit50", "tofHit(cal.posECALHit, 50)")\
           # .Define("tofHit100", "tofHit(cal.posECALHit, 100)")\
           # .Define("tofHit200", "tofHit(cal.posECALHit, 200)")\
           # .Define("tofHit300", "tofHit(cal.posECALHit, 300)")

    #selection of hits
    df = df.Define("sel_frank", "selectHits(dToLine, cal.layer, true, 10, 9999.)")\
           .Define("sel_cyl", "selectHits(dToLine, cal.layer, false, 10, 5.)")

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


# get_photons(df)
# get_pions(df)
beta_vs_p(df)
