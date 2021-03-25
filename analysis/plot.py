import ROOT
import numpy as np

ROOT.EnableImplicitMT(2)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetOptStat(1110)
ROOT.gInterpreter.Declare('#include "plot.hpp"')

def plot_pions():
    canvas = ROOT.TCanvas()
    df = ROOT.RDataFrame("pions", "/nfs/dust/ilc/user/dudarboh/analysis/pions.root")

    smearings = [0, 10, 30, 50, 100, 200, 300]
    methods = ["Closest", "Fastest", "Frank", "Cyl"]
    for s in smearings:
        for m in methods:
            df = df.Define("beta{1}{0}".format(s, m), "lengthCalo/(tof{1}{0}*SPEED_OF_LIGHT)".format(s, m))
            df = df.Define("betaCalib{1}{0}".format(s, m), "lengthCalo/(calibrateTOF(tof{1}{0}, nHitsCluster, \"{1}\")*SPEED_OF_LIGHT)".format(s, m))
            df = df.Define("mass{1}{0}".format(s, m), "mom / beta{1}{0} * sqrt(1. - beta{1}{0} * beta{1}{0}) * 1000".format(s, m))
            df = df.Define("massCalib{1}{0}".format(s, m), "mom / betaCalib{1}{0} * sqrt(1. - betaCalib{1}{0} * betaCalib{1}{0}) * 1000".format(s, m))

    histos = []
    for idx, s in enumerate(smearings):
        histos.append(df.Histo1D(("h{}".format(s), "{} ps".format(s),500, 0., 500.), "massCalibClosest{}".format(s)))
        histos[idx].SetLineColor(idx+1)
        histos[idx].SetLineWidth(2)

    histos[0].Draw()
    for i in range(1, 7):
        histos[i].Draw("same")

    canvas.BuildLegend()

    x = np.array([139.570183535, 493.677, 938.2720881629])
    line_pion = ROOT.TLine(x[0], 0., x[0], 1000.)
    line_pion.SetLineStyle(9)
    line_pion.Draw()

    # line_kaon = ROOT.TLine(x[1], 0., x[1], 1000.)
    # line_kaon.SetLineStyle(9)
    # line_kaon.Draw()
    #
    # line_proton = ROOT.TLine(x[2], 0., x[2], 1000.)
    # line_proton.SetLineStyle(9)
    # line_proton.Draw()

    # latex = ROOT.TLatex()
    # latex.SetTextSize(0.08)
    # latex.DrawLatex(x[0] + 10, 1000,"#pi^{#pm}");
    # latex.DrawLatex(x[1] + 10, 1000,"K^{#pm}");
    # latex.DrawLatex(x[2] + 10, 1000,"p^{#pm}");
    canvas.Update()
    input("wait")

def plot_photons():
    df = ROOT.RDataFrame("photons", "/nfs/dust/ilc/user/dudarboh/analysis/photons.root")
    # fit slices
    # fit slices
    # fit slices
    # canvas = ROOT.TCanvas()
    # canvas.Divide(2, 2)

    # df = df.Define("ratioClosest", "tofClosest0/tofTrue")\
    #        .Define("ratioFastest", "tofFastest0/tofTrue")\
    #        .Define("ratioFrank", "tofFrank0/tofTrue")\
    #        .Define("ratioCyl", "tofCyl0/tofTrue")\
    #
    # # canvas.cd(1)
    # h = df.Histo2D(("h"," Cylinder 5mm;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 200, 0.999, 1.001), "nHitsCluster","ratioCyl")
    # # h.Draw("colz")
    #
    # h.FitSlicesY(0,0,-1,100)
    #
    # # canvas.cd(2)
    # # h_chi2 = ROOT.gROOT.FindObject("h_chi2")
    # # h_chi2.Draw()
    #
    # # canvas.cd(3)
    # h_1 = ROOT.gROOT.FindObject("h_1")
    # h_1.Draw()
    #
    # # canvas.cd(4)
    # # h_2 = ROOT.gROOT.FindObject("h_2")
    # # h_2.Draw()
    # canvas.Update()
    # input("wait")

    canvas = ROOT.TCanvas()

    smearings = [0, 10, 30, 50, 100, 200, 300]
    methods = ["Closest", "Fastest", "Frank", "Cyl"]
    for s in smearings:
        for m in methods:
            # ratios
            df = df.Define("ratio{1}{0}".format(s, m), "tof{1}{0}/tofTrue".format(s, m))\
                   .Define("ratioCalib{1}{0}".format(s, m), "calibrateTOF(tof{1}{0}, nHitsCluster, \"{1}\")/tofTrue".format(s, m))\
            # deltas
            df = df.Define("dt{1}{0}".format(s, m), "tof{1}{0} - tofTrue".format(s, m))\
                   .Define("dtCalib{1}{0}".format(s, m), "calibrateTOF(tof{1}{0}, nHitsCluster, \"{1}\") - tofTrue".format(s, m))\


    h_before = df.Histo2D(("h_before"," before calibration;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 200, 0.999, 1.001), "nHitsCluster","ratioFrank0")
    h_before.FitSlicesY(0,0,-1,100)
    h_before_1 = ROOT.gROOT.FindObject("h_before_1")
    h_before_1.Draw()
    h_after = df.Histo2D(("h_after"," After calibration;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 200, 0.999, 1.001), "nHitsCluster","ratioCalibFrank0")
    h_after.FitSlicesY(0,0,-1,100)
    h_after_1 = ROOT.gROOT.FindObject("h_after_1")
    h_after_1.Draw("same")



    ################## 1D biases OF ALL METHODS CALIBRATED###############
    # h_closest = df.Histo1D(("h_closest"," Closest hit; #Delta TOF, [ns]; N_{PFO}", 400, -0.04, 0.04), "dtCalibClosest0")
    # h_fastest = df.Histo1D(("h_fastest"," Fastest hit; #Delta TOF, [ns]; N_{PFO}", 400, -0.04, 0.04), "dtCalibFastest0")
    # h_frank = df.Filter("ratioFrank>0").Histo1D(("h_frank"," Frank; #Delta TOF, [ns]; N_{PFO}", 400, -0.04, 0.04), "dtCalibFrank0")
    # h_cyl = df.Filter("ratioCyl>0").Histo1D(("h_cyl"," 5 mm cyl; #Delta TOF, [ns]; N_{PFO}", 400, -0.04, 0.04), "dtCalibCyl0")
    # h_closest.SetLineColor(1)
    # h_fastest.SetLineColor(2)
    # h_frank.SetLineColor(4)
    # h_cyl.SetLineColor(6)
    #
    # h_closest.SetLineWidth(3)
    # h_fastest.SetLineWidth(3)
    # h_frank.SetLineWidth(3)
    # h_cyl.SetLineWidth(3)
    #
    # h_closest.Draw()
    # h_fastest.Draw("same")
    # h_frank.Draw("same")
    # h_cyl.Draw("same")

    ################## Bias of 1 method vs smearing###############
    # h0 = df.Histo1D(("h0"," 0 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl0")
    # h10 = df.Histo1D(("h10"," 10 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl10")
    # h30 = df.Histo1D(("h30"," 30 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl30")
    # h50 = df.Histo1D(("h50"," 50 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl50")
    # h100 = df.Histo1D(("h100"," 100 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl100")
    # h200 = df.Histo1D(("h200"," 200 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl200")
    # h300 = df.Histo1D(("h300"," 300 ps; #Delta TOF, [ns]; N_{PFO}", 50000, -5., 5.), "dtCyl300")
    # histos = [h0, h10, h30, h50, h100, h200, h300]
    # for idx, histo in enumerate(histos):
    #     histo.SetLineColor(idx+1)
    #     histo.SetLineWidth(2)
    #     print(histo.GetStdDev())
    # histos[0].Draw()
    # for i in range(1, 7):
    #     histos[i].Draw("same")

    # BIAS BEFORE/AFTER CALLIBRATION
    # h_before = df.Histo1D(("h_before","Before calibration; #Delta TOF, [ns]; N_{PFO}", 400, -0.04, 0.04), "dtFrank0")
    # h_after = df.Histo1D(("h_after","After calibration; #Delta TOF, [ns]; N_{PFO}", 400, -0.04, 0.04), "dtCalibFrank0")
    # h_after.SetLineColor(2)
    # h_before.SetLineWidth(3)
    # h_after.SetLineWidth(3)
    # h_before.Draw()
    # h_after.Draw("sames")


    # RATIO PROFILES
    # h_closest = df.Histo2D(("h_closest"," Closest hit;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster","ratioClosest").ProfileX()
    # h_fastest = df.Histo2D(("h_fastest"," Fastest hit;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster", "ratioFastest").ProfileX()
    # h_frank = df.Filter("ratioFrank>0").Histo2D(("h_frank"," closest to track;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster", "ratioFrank").ProfileX()
    # h_cyl = df.Filter("ratioCyl>0").Histo2D(("h_cyl"," 5 mm cyl;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster", "ratioCyl").ProfileX()

    # BIAS PROFILES
    # h_closest = df.Histo2D(("h_closest"," Closest hit;n Hits in ECAL Cluster;#Delta TOF, [ns]", 100, 0, 100, 400, -0.04, 0.04), "nHitsCluster","dtClosest").ProfileX()
    # h_fastest = df.Histo2D(("h_fastest"," Fastest hit;n Hits in ECAL Cluster;#Delta TOF, [ns]", 100, 0, 100, 400, -0.04, 0.04), "nHitsCluster", "dtFastest").ProfileX()
    # h_frank = df.Histo2D(("h_frank"," closest to track;n Hits in ECAL Cluster;#Delta TOF, [ns]", 100, 0, 100, 400, -0.04, 0.04), "nHitsCluster", "dtFrank").ProfileX()
    # h_cyl = df.Histo2D(("h_cyl"," 5 mm cyl;n Hits in ECAL Cluster;#Delta TOF, [ns]", 100, 0, 100, 400, -0.04, 0.04), "nHitsCluster", "dtCyl").ProfileX()

    # RATIO CALIBRATED
    # h_closest = df.Histo2D(("h_closest"," Closest hit;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster","ratioCalibClosest")
    # h_fastest = df.Histo2D(("h_fastest"," Fastest hit;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster", "ratioCalibFastest")
    # h_frank = df.Filter("ratioFrank>0").Histo2D(("h_frank"," closest to track;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster", "ratioCalibFrank")
    # h_cyl = df.Filter("ratioCyl>0").Histo2D(("h_cyl"," 5 mm cyl;n Hits in ECAL Cluster;TOF_{reco}/TOF_{true}", 500, 0, 500, 100000, 0., 2.), "nHitsCluster", "ratioCalibCyl")
    #
    #


    # h_closest.Draw("colz")
    # canvas.Update()
    # input("wait")
    # h_fastest.Draw("colz")
    # canvas.Update()
    # input("wait")
    # h_frank.Draw("colz")
    # canvas.Update()
    # input("wait")
    # h_cyl.Draw("colz")
    # canvas.Update()
    # input("wait")

    canvas.BuildLegend()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.Update()
    input("wait")




# plot_photons()
# plot_pions()
