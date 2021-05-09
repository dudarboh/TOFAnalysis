import ROOT
import numpy as np
# ROOT.EnableImplicitMT(4)

def plot_spectrum():
    canvas = ROOT.TCanvas("name", "title", 800, 800)

    # methods = [
    # "pTrackAtIP_lengthTrackIP_Frank_0.0",
    # "pTrackAtIP_lengthTrackIP_FrankAvg_0.0",
    # "pTrackAtIP_lengthTrackIP_Closest_0.0",
    # ]

    # methods = [
    # "pTrackAtCalo_lengthTrackCalo_Frank_0.0",
    # "pTrackAtCalo_lengthTrackCalo_FrankAvg_0.0",
    # "pTrackAtCalo_lengthTrackCalo_Closest_0.0",
    # ]

    methods = [
    "pTrackAtIP_lengthTrackIP_FrankAvg_0.0",
    "pTrackAtCalo_lengthTrackCalo_FrankAvg_0.0",
    "pTrackAtCalo_lengthTrackCalo_FrankAvg_0.0_2.5GeVCut",
    "pTrackAtIP_lengthTrackIP_FrankAvg_0.0_2.5GeVCut",
    ]


    histos = []
    for i, m in enumerate(methods):
        if i == 2:
            df = ROOT.RDataFrame("pTrackAtCalo_lengthTrackCalo_FrankAvg_0.0", "./final_roots/" + m + ".root")
        elif i == 3:
            df = ROOT.RDataFrame("pTrackAtIP_lengthTrackIP_FrankAvg_0.0", "./final_roots/" + m + ".root")
        else:
            df = ROOT.RDataFrame(m, "./final_roots/" + m + ".root")
        h = df.Define("mass", "mom/beta * sqrt(1 - beta*beta)*1000").Histo1D(("h{}".format(i), "h{}".format(i), 1000, 0, 1000), "mass")
        h = df.Histo1D(("h{}".format(i), "h{}".format(i), 1000, 0, 10), "mom")
        histos.append(h)

    for i, h in enumerate(histos):
        print(h.GetEntries())
        h.Draw("" if i == 0 else "same")
        h.SetLineColor(i+1)


    canvas.BuildLegend()
    # x = np.array([139.57018, 493.677, 938.2720881629])
    # lines = []
    # for i in x:
    #     lines.append(ROOT.TLine(i, 0., i, 12000.))
    #     lines[-1].SetLineStyle(9)
    #     lines[-1].SetLineWidth(3)
    #     lines[-1].Draw()
    #
    # latex = ROOT.TLatex()
    # latex.DrawLatex(x[0] - 10, 11000,"#pi");
    # latex.DrawLatex(x[1] - 10, 11000,"K");
    # latex.DrawLatex(x[2] - 10, 11000,"p");

    canvas.Update()
    input("wait")

def plot_bias():
    canvas = ROOT.TCanvas("name", "title", 800, 800)

    x = np.array([139.57018, 493.677, 938.2720881629])
    bias_fit_ip = np.array([1.14871e+02, 4.84824e+02, 9.33769e+02]) - x
    bias_avg_ip = np.array([1.22656e+02, 4.93565e+02, 9.46643e+02]) - x
    bias_closest_ip = np.array([1.17873e+02, 4.87669e+02, 9.37008e+02]) - x
    bias_fit_calo = np.array([1.49989e+02, 4.96020e+02, 9.40080e+02]) - x
    bias_avg_calo = np.array([1.55833e+02, 5.04787e+02, 9.53494e+02]) - x
    bias_closest_calo = np.array([1.51774e+02, 4.98727e+02, 9.43518e+02]) - x

    sigma_fit_ip = np.array([1.18564e-01, 1.61536e-01, 1.27877e-01])
    sigma_avg_ip = np.array([1.43680e-01, 1.83707e-01, 1.39163e-01])
    sigma_closest_ip = np.array([8.62805e-02, 1.53674e-01, 9.14025e-02])
    sigma_fit_calo = np.array([1.20895e-01, 1.66740e-01, 9.22562e-02])
    sigma_avg_calo = np.array([1.79584e-01, 1.93135e-01, 1.21782e-01])
    sigma_closest_calo = np.array([8.43831e-02, 9.65479e-02, 9.73941e-02])

    biases = [bias_fit_ip, bias_avg_ip, bias_closest_ip, bias_fit_calo, bias_avg_calo, bias_closest_calo]

    graphs = []
    for i, b in enumerate(biases):
        graphs.append(ROOT.TGraph(3, x, b))
        graphs[i].SetLineColor(i+1)
        graphs[i].Draw("APL" if i == 0 else "PLsame")

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

plot_spectrum()
