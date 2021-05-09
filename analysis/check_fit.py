import ROOT
import sys
ROOT.gInterpreter.Declare('#include "extract.hpp"')

ch = ROOT.TChain("TOFAnalysis")

# # 2f_Z_hadronic data
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

canvas = ROOT.TCanvas("")
canvas.Divide(2, 1)

for i, event in enumerate(ch):
    if not (event.nECALHits > 0 and abs(event.xyzCluster.Z() < 2200.) and abs(event.xyzVtxMC.R() < .5)):
        continue
    canvas.cd(1)
    x0 = event.xyzTrackAtCalo.X()
    y0 = event.xyzTrackAtCalo.Y()
    px = event.pTrackAtCalo.X()
    py = event.pTrackAtCalo.Y()
    gr_line = ROOT.TF1("f1", "{0} + {1}/{2}*(x - {3})".format(y0, py, px, x0), -3000., 3000.)
    gr_point = ROOT.TGraph()
    gr_point.SetPoint(0, x0, y0)
    gr_point.SetMarkerStyle(29)
    gr_point.SetMarkerColor(4)
    gr_point.SetMarkerSize(2.5)

    gr_shower = ROOT.TGraph()
    for j, hit in enumerate(event.xyzECALHit):
        gr_shower.SetPoint(j, hit.X(), hit.Y())

    gr_shower.Draw("AP")
    gr_shower.SetTitle("Shower vs track extrapolation;x, [mm];y, [mm]")
    gr_shower.SetMarkerStyle(20)
    gr_line.Draw("same")
    gr_point.Draw("Psame")
    gr_shower.GetXaxis().SetRangeUser(-3000., 3000.)
    gr_shower.GetYaxis().SetRangeUser(-3000., 3000.)
    canvas.Update()
    input("wait")
    gr_fit = ROOT.TGraph()
    for i, (pos, t) in enumerate( zip(event.xyzECALHit, event.tECALHit) ):
        gr_fit.SetPoint(i, (pos - event.xyzTrackAtCalo).R(), t)
    gr_fit.Fit("pol1", "Q");
    gr_fit.SetLineStyle(20)
    gr_fit.Draw("AP")



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
               .Define("tof", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 0)")\
               .Define("slope", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 1)")
