import ROOT
import time
import numpy as np

ROOT.gStyle.SetPalette(1)

# 2f_Z_hadronic data
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

canvas = ROOT.TCanvas()
df = ROOT.RDataFrame(ch)
ROOT.gInterpreter.Declare('#include "analysis.hpp"')

# Check only good photons in barrel
df = df.Filter("nMC == 1 && PDG[0] == 22 && cal.nHits > 0 \
                         && isBackscatter[0] == 0 && isDecayedInTracker[0] == 0\
                         && abs(cluster.posCluster.Z()) < 2200. && vtxMC[0].R() < 0.5")\
        .Range(800000)

# Find intersection of track with ECAL plane
df = df.Define("ECALPlane", "getECALPlane( pMC[0].Phi() )")\
       .Define("rImpact", "intersection(vtxMC[0].Vect(), pMC[0].Vect(), ECALPlane)")

# Find true tof based on track length
df = df.Define("tofTrue", "(rImpact - vtxMC[0].Vect()).R() / SPEED_OF_LIGHT")\
       .Define("dToImpact", "dToImpact(cal.posECALHit, rImpact)")\
       .Define("tofClosest", "tofClosest(cal.posECALHit, dToImpact)")\
       .Define("tofFastest", "tofFastest(cal.posECALHit, dToImpact)")

df = df.Define("dToLine", "dToLine(cal.posECALHit, vtxMC[0].Vect(), pMC[0].Vect())")\
       .Define("sel", "selectHits(dToLine, cal.layer)")\
       .Define("tofAvg", "tofAvg(cal.posECALHit[sel], dToImpact[sel])")\
       .Define("tofFit", "fitFunc(cal.posECALHit[sel], dToImpact[sel], 0)")\
       .Define("slopeFit", "fitFunc(cal.posECALHit[sel], dToImpact[sel], 1)")\
       .Define("chi2Fit", "fitFunc(cal.posECALHit[sel], dToImpact[sel], 2)")\
       .Define("ndfFit", "fitFunc(cal.posECALHit[sel], dToImpact[sel], 3)")\
       .Filter("ndfFit > 0")\
       .Define("chi2ndf", "chi2Fit/ndfFit")\
       .Define("residual", "residual(cal.posECALHit[sel], dToImpact[sel], tofFit, slopeFit)")\
       .Define("nHitsFit", "residual.size()")

# hFastest = df.Define("dtFastest", "tofFastest - tofTrue").Histo1D(("hFastest","TOF fastest; TOF, [ns]; N_{PFO}", 1000, -.1, .1), "dtFastest")
# hClosest = df.Define("dtClosest", "tofClosest - tofTrue").Histo1D(("hClosest","TOF closest; TOF, [ns]; N_{PFO}", 1000, -.1, .1), "dtClosest")
# hAvg = df.Define("dtAvg", "tofAvg - tofTrue").Histo1D(("hAvg","TOF avg; TOF, [ns]; N_{PFO}", 1000, -.1, .1), "dtAvg")
# hFit = df.Define("dtFit", "tofFit - tofTrue").Histo1D(("hFit","TOF fit; TOF, [ns]; N_{PFO}", 1000, -.1, .1), "dtFit")
#
# # 1 - fastest, 4 - closest, 6 - fit, 8 - avg
# hFastest.Draw()
# hClosest.Draw("sames")
# hClosest.SetLineColor(4)
# hFit.Draw("sames")
# hFit.SetLineColor(6)
# hAvg.Draw("sames")
# hAvg.SetLineColor(8)


# h_good = df.Filter("tofFit - tofTrue >= -0.008").Profile1D(("h_good", "Good PFOs; residual, [ns]; avg(E_{hit}), [GeV]", 500, -.1, .1), "residual", "energyFit")
# h_bad = df.Filter("tofFit - tofTrue < -0.008").Profile1D(("h_bad", "Bad PFOs; residual, [ns]; avg(E_{hit}), [GeV]", 500, -.1, .1), "residual", "energyFit")
# h_good = df.Filter("tofFit - tofTrue >= -0.008").Histo2D(("h_good", "Good PFOs; residual, [ns]; slope, [ns/mm]", 500, -.04, .04, 500, 0.002, .007), "residual", "slopeFit")
# h_bad = df.Filter("tofFit - tofTrue < -0.008").Histo2D(("h_bad", "Bad PFOs; residual, [ns];  slope, [ns/mm]", 500, -.04, .04, 500, 0.002, .007), "residual", "slopeFit")
h_good = df.Filter("tofFit - tofTrue >= -0.008").Histo2D(("h_good", "Good PFOs; nHits used for fit; slope, [ns/mm]", 13, 0, 13, 500, 0.002, .007), "nHitsFit", "slopeFit")
h_bad = df.Filter("tofFit - tofTrue < -0.008").Histo2D(("h_bad", "Bad PFOs; nHits used for fit;  slope, [ns/mm]", 13, 0, 13, 500, 0.002, .007), "nHitsFit", "slopeFit")


h_good.Draw("colz")
canvas.Update()
input("wait")
h_bad.Draw("colz")
# h_bad.SetLineColor(2)
# canvas.BuildLegend()
# canvas.SetGridx()
# canvas.SetGridy()
# hFastest.SetTitle("Histo title")
canvas.Update()
input("wait")



def analysis():
    df = selection(df_initial).Range(100000)




# analysis()
