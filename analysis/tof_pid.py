import ROOT
import numpy as np
# ROOT.EnableImplicitMT(4)

methods = [
"pTrackAtIP_lengthTrackIP_FrankAvg_0.0",
# "pTrackAtIP_lengthTrackIP_FrankAvg_10.0",
# "pTrackAtIP_lengthTrackIP_FrankAvg_50.0",
# "pTrackAtIP_lengthTrackIP_FrankAvg_100.0",
]

# methods = [
# "pTrackAtCalo_lengthTrackCalo_Frank_0.0",
# "pTrackAtCalo_lengthTrackCalo_Frank_10.0",
# "pTrackAtCalo_lengthTrackCalo_Frank_50.0",
# "pTrackAtCalo_lengthTrackCalo_Frank_100.0",
# ]


def beta_vs_p(tree_name):
    canvas = ROOT.TCanvas("img_beta_vs_p_{}".format(tree_name), "beta_vs_p_{}".format(tree_name) ,800, 800)
    df = ROOT.RDataFrame(tree_name, "./final_roots/" + tree_name + ".root")
    # gr_bg = df.Filter("abs(PDG) != 2212 && abs(PDG) != 321 && abs(PDG) != 211").Graph("mom", "beta")
    # gr_all = df.Graph("mom", "beta")
    gr_pion = df.Filter("abs(PDG) == 211").Graph("mom", "beta")
    gr_kaon = df.Filter("abs(PDG) == 321").Graph("mom", "beta")
    gr_proton = df.Filter("abs(PDG) == 2212").Graph("mom", "beta")


    # gr_all.Draw("AP")
    # gr_all.SetMarkerColor(1)
    # gr_all.SetMarkerStyle(1)
    # gr_all.SetTitle("Background;momentum, [GeV];#beta")
    # gr_all.SetLineWidth(1)
    # gr_all.GetXaxis().SetRangeUser(1., 10.)
    # gr_all.GetYaxis().SetRangeUser(0.8, 1.01)

    gr_pion.Draw("AP")
    gr_pion.SetMarkerColor(2)
    gr_pion.SetMarkerStyle(1)
    gr_pion.SetTitle("Pions;")
    gr_pion.SetLineColor(2)

    gr_kaon.Draw("Psame")
    gr_kaon.SetMarkerColor(3)
    gr_kaon.SetMarkerStyle(1)
    gr_kaon.SetTitle("Kaons")
    gr_kaon.SetLineColor(3)

    gr_proton.Draw("Psame")
    gr_proton.SetMarkerColor(4)
    gr_proton.SetMarkerStyle(1)
    gr_proton.SetTitle("Protons")
    gr_proton.SetLineColor(4)

    # gr_bg.Draw("Psame")
    # gr_bg.SetMarkerColor(1)
    # gr_bg.SetTitle("Background")
    # gr_bg.SetLineWidth(1)
    # gr_bg.SetLineColor(1)
    canvas.BuildLegend(0.5, 0.2, 0.8, 0.5)
    gr_pion.SetTitle("")
    canvas.Update()
    input("wait")

for m in methods:
    beta_vs_p(m)


def band_fit(tree_name, pdg=211):
    df = ROOT.RDataFrame(tree_name, "./final_roots/" + tree_name + ".root")

    if pdg == 211:
        h = df.Filter("abs(PDG) == {}".format(pdg)).Histo2D(("h{}".format(pdg), "PDG: {}".format(pdg), 90, 1., 10., 2000, 0.95, 1.05 ), "mom", "beta")
    elif pdg == 321:
        h = df.Filter("abs(PDG) == {}".format(pdg)).Histo2D(("h{}".format(pdg), "PDG: {}".format(pdg), 90, 1., 10., 500, 0.85, 1.05 ), "mom", "beta")
    elif pdg == 2212:
        h = df.Filter("abs(PDG) == {}".format(pdg)).Histo2D(("h{}".format(pdg), "PDG: {}".format(pdg), 90, 1., 10., 250, 0.7, 1.05 ), "mom", "beta")

    h.FitSlicesY()
    h_1 = ROOT.gROOT.FindObject("h{}_1".format(pdg))
    h_2 = ROOT.gROOT.FindObject("h{}_2".format(pdg))
    h_chi2 = ROOT.gROOT.FindObject("h{}_chi2".format(pdg))

    # canvas = ROOT.TCanvas("fit_slices_{}".format(pdg))
    # canvas.Divide(2, 2)
    # canvas.cd(1)
    # h.Draw("colz")
    # canvas.cd(2)
    # h_chi2.Draw()
    # canvas.cd(3)
    # h_1.Draw()
    # canvas.cd(4)
    # h_2.Draw()
    # canvas.Update()

    # Log Bins code. Might be needed later...
    # input("wait")
    # nBinsX = 100
    # nBinsY = 200
    # # log10(min/max momentum / GeV)
    # minBinX = 0
    # maxBinX = 1
    # # for TOF beta = v/c values from 0 to 1
    # minBinY = 0.7
    # maxBinY = 1.05
    # histbinsX = []
    # histbinsY = []
    # for i in range(nBinsX +1):
    #     histbinsX.append( pow( 10, minBinX + (maxBinX-minBinX)*i/nBinsX ) )
    # for i in range(nBinsY +1):
    #     histbinsY.append( minBinY+ (maxBinY-minBinY)*i/nBinsY )


    return h, h_1, h_2

def sep_power(h_mean1, h_mean2, h_std1, h_std2):
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(20)
    for i in range(1, h_mean1.GetNbinsX()+1):
        m1 = h_mean1.GetBinContent(i)
        m2 = h_mean2.GetBinContent(i)
        std1 = h_std1.GetBinContent(i)
        std2 = h_std2.GetBinContent(i)

        sep_power = abs(m1-m2) / np.sqrt( (std1*std1 + std2*std2)/2. )
        gr.SetPoint(i-1, h_mean1.GetBinCenter(i), sep_power)

    # gr.Draw("AP")
    # input("wait")
    return gr



def compare_slope(tree_name):
    canvas = ROOT.TCanvas("beta_vs_p_{}".format(tree_name))
    df = ROOT.RDataFrame(tree_name, "./final_roots/" + tree_name + ".root")
    # gr_bg = df.Filter("abs(PDG) != 2212 && abs(PDG) != 321 && abs(PDG) != 211").Graph("mom", "beta")
    h_pion = df.Filter("abs(PDG) == 211").Histo1D(("h_pion", "Pions", 1000, -0.01, 0.03), "slope")
    h_kaon = df.Filter("abs(PDG) == 321").Histo1D(("h_kaon", "Kaons", 1000, -0.01, 0.03), "slope")
    h_proton = df.Filter("abs(PDG) == 2212").Histo1D(("h_protons", "Protons", 1000, -0.01, 0.03), "slope")
    h_pion.Draw()
    h_kaon.Draw("same")
    h_kaon.SetLineColor(2)
    h_proton.Draw("same")
    h_proton.SetLineColor(3)
    input("sait")

graphs_pik = []
graphs_kp = []
for m in methods:
    # compare_slope(m)
    h_pion, h_pion_mean, h_pion_std = band_fit(m, pdg=211)
    h_kaon, h_kaon_mean, h_kaon_std = band_fit(m, pdg=321)
    h_proton, h_proton_mean, h_proton_std = band_fit(m, pdg=2212)
    graphs_pik.append( sep_power(h_pion_mean, h_kaon_mean, h_pion_std, h_kaon_std) )
    graphs_kp.append( sep_power(h_kaon_mean, h_proton_mean, h_kaon_std, h_proton_std) )





canvas = ROOT.TCanvas("sep_power_plots")
mg = ROOT.TMultiGraph()
mg.SetTitle("#pi^{#pm}/K^{#pm} separation power for different TOF methods; p, [GeV]; Sep. Power, [#sigma]")
# mg.SetTitle("K^{#pm}/p^{#pm} separation power for different TOF methods; p, [GeV]; Sep. Power, [#sigma]")
for i, gr in enumerate(graphs_pik):
    gr.SetTitle(methods[i])
    gr.SetLineColor(i+1)
    gr.SetMarkerColor(i+1)
    gr.SetLineWidth(1)
    mg.Add(gr, "PL")

mg.Draw("A")
canvas.SetGridx()
canvas.SetGridy()
canvas.BuildLegend(0.5, 0.2, 0.8, 0.5)
canvas.Update()
input("wait")
