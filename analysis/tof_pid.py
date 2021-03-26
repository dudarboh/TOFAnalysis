import ROOT
import numpy as np

ROOT.gStyle.SetPalette(1)

# # 2f_Z_hadronic data
ch = ROOT.TChain("TOFAnalysis")
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")
ROOT.gInterpreter.Declare('#include "extract.hpp"')



def get_beta(momentum="pTrackAtCalo", length="lengthTrackCalo", tof="Closest", smearing=0.):
    df = ROOT.RDataFrame(ch)
    df = df.Range(50000).Filter("nECALHits > 0")

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
        # FIXME: Change to average instead of the fit
        df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
               .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
               .Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
               .Define("tof", "fitFunc(tofHit[sel_frank], dToImpact[sel_frank], 0)")

    df = df.Define("beta", "{}/(tof * SPEED_OF_LIGHT)".format(length))
    return df

def beta_vs_p(df):
    canvas = ROOT.TCanvas("beta_vs_p")
    gr_pion = df.Filter("PDG == 211").Graph("mom", "beta")
    gr_pion.SetMarkerColor(2)
    gr_pion.SetTitle("Pions")
    gr_pion.SetLineWidth(5)
    gr_pion.SetLineColor(2)

    gr_kaon = df.Filter("PDG == 321").Graph("mom", "beta")
    gr_kaon.SetMarkerColor(3)
    gr_kaon.SetTitle("Kaons")
    gr_kaon.SetLineWidth(5)
    gr_kaon.SetLineColor(3)

    gr_proton = df.Filter("PDG == 2212").Graph("mom", "beta")
    gr_proton.SetMarkerColor(4)
    gr_proton.SetTitle("Protons")
    gr_proton.SetLineWidth(5)
    gr_proton.SetLineColor(4)

    gr_bg = df.Filter("PDG != 2212 && PDG != 321 && PDG != 211").Graph("mom", "beta")
    gr_bg.SetMarkerColor(1)
    gr_bg.SetTitle("Background")
    gr_bg.SetLineWidth(5)
    gr_bg.SetLineColor(1)

    gr_bg.Draw("AP")
    gr_bg.GetXaxis().SetRangeUser(1., 10.)
    gr_bg.GetYaxis().SetRangeUser(0.7, 1.1)

    gr_pion.Draw("Psame")
    gr_kaon.Draw("Psame")
    gr_proton.Draw("Psame")


    canvas.Update()
    return canvas, gr_bg, gr_pion, gr_kaon, gr_proton


def check_band_fit(df, pdg=211):
    canvas = ROOT.TCanvas("c{}".format(pdg))
    canvas.Divide(2, 2)
    h = df.Filter("abs(PDG) == {}".format(pdg)).Histo2D(("h{}".format(pdg), "PDG: {}".format(pdg), 100, 1., 10., 2000, 0.7, 1.1 ), "mom", "beta")

    canvas.cd(1)
    h.Draw("colz")
    h.FitSlicesY()
    # h.FitSlicesY(0,0,-1,100)

    canvas.cd(2)
    h_chi2 = ROOT.gROOT.FindObject("h{}_chi2".format(pdg))
    h_chi2.Draw()

    canvas.cd(3)
    h_1 = ROOT.gROOT.FindObject("h{}_1".format(pdg))
    h_1.Draw()

    canvas.cd(4)
    h_2 = ROOT.gROOT.FindObject("h{}_2".format(pdg))
    h_2.Draw()
    canvas.Update()
    return canvas, h, h_1, h_2

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


def sep_power(h_mean1, h_mean2, h_std1, h_std2):
    canvas = ROOT.TCanvas("sep_power")
    print(h_mean1.GetNbinsX())
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(20)
    for i in range(1, h_mean1.GetNbinsX()+1):
        m1 = h_mean1.GetBinContent(i)
        m2 = h_mean2.GetBinContent(i)
        std1 = h_std1.GetBinContent(i)
        std2 = h_std2.GetBinContent(i)

        sep_power = abs(m1-m2) / np.sqrt( (std1*std1 + std2*std2)/2. )
        gr.SetPoint(i, h_mean1.GetBinCenter(i), sep_power)

    gr.Draw("AP")
    return canvas, gr

df = get_beta()
c_beta_vs_p, gr_bg, gr_pion, gr_kaon, gr_proton = beta_vs_p(df)
c_pion, h_pion, h_pion_mean, h_pion_std = check_band_fit(df, pdg=211)
c_kaon, h_kaon, h_kaon_mean, h_kaon_std = check_band_fit(df, pdg=321)
c_proton, h_proton, h_proton_mean, h_proton_std = check_band_fit(df, pdg=2212)

c_sep_power, gr = sep_power(h_pion_mean, h_kaon_mean, h_pion_std, h_kaon_std)

input("wait")
