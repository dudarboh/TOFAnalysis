import ROOT
import sys
import numpy as np
ROOT.gInterpreter.Declare('#include "extract.hpp"')
canvas = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(10)


def draw_mass_spectrum():

    tree_name = "pTrackAtCalo_lengthTrackCalo_FrankAvg_100.0"
    ch = ROOT.TChain(tree_name)
    ch.Add("./final_roots/" + tree_name + ".root")
    df = ROOT.RDataFrame(ch)

    h1 = df.Filter("abs(PDG) == 211 || abs(PDG) == 321 || abs(PDG) == 2212").Define("m", "mom/beta * sqrt(1. - beta* beta) * 1000").Histo1D(("h1", "title", 1000, 0, 1000), "m")
    h1.SetTitle("TOF_{avg}; reconstructed mass, [MeV]; N PFOs")
    h1.SetLineWidth(3)
    h1.SetLineColor(2)
    h1.Draw()

    # tree_name = "pTrackAtCalo_lengthTrackCalo_Cyl_100.0"
    # ch = ROOT.TChain(tree_name)
    # ch.Add("./final_roots/" + tree_name + ".root")
    # df = ROOT.RDataFrame(ch)
    #
    # h11 = df.Filter("abs(PDG) == 211 || abs(PDG) == 321 || abs(PDG) == 2212").Define("m", "mom/beta * sqrt(1. - beta* beta) * 1000").Histo1D(("h11", "title", 1000, 0, 1000), "m")
    # h11.SetTitle("TOF_{cyl}; reconstructed mass, [MeV]; N PFOs")
    # h11.SetLineWidth(3)
    # h11.SetLineColor(11)
    # h11.Draw("same")


    tree_name = "pTrackAtCalo_lengthTrackCalo_Frank_100.0"
    ch = ROOT.TChain(tree_name)
    ch.Add("./final_roots/" + tree_name + ".root")
    df = ROOT.RDataFrame(ch)

    h2 = df.Filter("abs(PDG) == 211 || abs(PDG) == 321 || abs(PDG) == 2212").Define("m", "mom/beta * sqrt(1. - beta* beta) * 1000").Histo1D(("h2", "title", 1000, 0, 1000), "m")
    h2.SetTitle("TOF_{fit}; reconstructed mass, [MeV]; N PFOs")
    h2.SetLineWidth(3)
    h2.SetLineColor(6)
    h2.Draw("same")


    tree_name = "pTrackAtCalo_lengthTrackCalo_Closest_100.0"
    ch = ROOT.TChain(tree_name)
    ch.Add("./final_roots/" + tree_name + ".root")
    df = ROOT.RDataFrame(ch)

    h3 = df.Filter("abs(PDG) == 211 || abs(PDG) == 321 || abs(PDG) == 2212").Define("m", "mom/beta * sqrt(1. - beta* beta) * 1000").Histo1D(("h3", "title", 1000, 0, 1000), "m")
    h3.SetTitle("TOF_{closest}; reconstructed mass, [MeV]; N PFOs")
    h3.SetLineWidth(3)
    h3.SetLineColor(7)
    h3.Draw("same")

    tree_name = "pTrackAtCalo_lengthTrackCalo_Fastest_100.0"
    ch = ROOT.TChain(tree_name)
    ch.Add("./final_roots/" + tree_name + ".root")
    df = ROOT.RDataFrame(ch)

    h4 = df.Filter("abs(PDG) == 211 || abs(PDG) == 321 || abs(PDG) == 2212").Define("m", "mom/beta * sqrt(1. - beta* beta) * 1000").Histo1D(("h4", "title", 1000, 0, 1000), "m")
    h4.SetTitle("TOF_{fastest}; reconstructed mass, [MeV]; N PFOs")
    h4.SetLineWidth(3)
    h4.SetLineColor(8)
    h4.Draw("same")



    x = np.array([139.57018, 493.677, 938.2720881629])
    canvas.BuildLegend()
    lines = []
    for i in x:
        lines.append(ROOT.TLine(i, 0., i, 12000.))
        lines[-1].SetLineStyle(9)
        lines[-1].SetLineWidth(3)
        lines[-1].Draw()

    latex = ROOT.TLatex()
    latex.DrawLatex(x[0] - 10, 11000,"#pi");
    latex.DrawLatex(x[1] - 10, 11000,"K");
    latex.DrawLatex(x[2] - 10, 11000,"p");

    canvas.Update()
    input("wait")

def compare_masses():
    x = np.array([139.57018, 493.677, 938.2720881629])

    # ip lenIP ip avg
    y1 = np.array([1.22747e+02, 4.93306e+02, 9.46472e+02]) - x
    y1_err = np.array([1.55932e-01, 1.39194e-01, 1.16143e-01])
    gr1 = ROOT.TGraphErrors(3, x, y1, ROOT.nullptr, y1_err)
    gr1.Draw()
    gr1.SetTitle("l_{track}(#Omega_{IP}, #lambda_{IP}); m_{true}, [Mev]; m_{reco} - m_{true}, [MeV]")
    gr1.SetMarkerStyle(20)
    gr1.SetLineWidth(3)

    # ip lenCalo ip avg
    y2 = np.array([1.56571e+02, 5.06042e+02, 9.55616e+02]) - x
    y2_err = np.array([1.78989e-01, 2.08750e-01, 1.28547e-01])
    gr2 = ROOT.TGraphErrors(3, x, y2, ROOT.nullptr, y2_err)
    gr2.Draw("same")
    gr2.SetTitle("l_{track}(#Omega_{calo}, #lambda_{calo}); m_{true}, [Mev]; m_{reco} - m_{true}, [MeV]")
    gr2.SetLineColor(2)
    gr2.SetMarkerColor(2)
    gr2.SetMarkerStyle(20)
    gr2.SetLineWidth(3)

    lines = []
    for i in x:
        lines.append(ROOT.TLine(i, -20., i, 20.))
        lines[-1].SetLineStyle(9)
        lines[-1].SetLineWidth(3)
        lines[-1].Draw()

    latex = ROOT.TLatex()
    latex.DrawLatex(x[0] + 10, 7,"#pi^{#pm}");
    latex.DrawLatex(x[1] + 10, 7,"K^{#pm}");
    latex.DrawLatex(x[2] + 10, 7,"p^{#pm}");


    input("wait")

def compare_smearings():
    x = np.array([0, 10, 50., 100.])
    y_frank = np.array([5.85, 5.18, 3.16, 2.32])
    y_frankAvg = np.array([5.92, 5.77, 4.32, 3.3])

    gr1 = ROOT.TGraph(len(x), x, y_frank)
    gr2 = ROOT.TGraph(len(x), x, y_frankAvg)

    gr1.SetTitle("TOF_{fit}; time resolution, [ps]; momentum at 2 #sigma, [GeV]")
    gr1.Draw("APL")
    gr2.Draw("PLsame")
    gr2.SetLineColor(2)
    input("wait")

# draw_mass_spectrum()
# compare_masses()
compare_smearings()
