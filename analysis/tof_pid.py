import ROOT
import numpy as np
# ROOT.EnableImplicitMT(4)


df_0 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/0ps_250gev.root")
df_10 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/10ps_250gev.root")
df_30 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/30ps_250gev.root")
df_50 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/50ps_250gev.root")
df_100 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/100ps_250gev.root")

df_names = [0 , 10, 30, 50, 100]
dfs = [df_0, df_10, df_30, df_50, df_100]

beta_set = "track_length_set/(set_hit_time*SPEED_OF_LIGHT)"
beta_closest = "track_length_calo/(tof_closest*SPEED_OF_LIGHT)"
beta_fastest = "track_length_calo/(tof_fastest*SPEED_OF_LIGHT)"
beta_avg = "track_length_calo/(tof_frank_avg*SPEED_OF_LIGHT)"
beta_fit = "track_length_calo/(tof_frank_fit*SPEED_OF_LIGHT)"

beta_names =["set", "closest", "fastest", "avg", "fit"]
betas = [beta_set, beta_closest, beta_fastest, beta_avg, beta_fit]
pdgs = [211, 321, 2212]

canvas = ROOT.TCanvas()
graphs = {}

for i, df in zip(df_names, dfs):
    for j, beta in zip(beta_names, betas):
        for pdg in pdgs:
            gr = ROOT.TGraphErrors()
            gr.SetMarkerStyle(20)
            h = df.Filter("n_ecal_hits > 0 && abs(pos_fastest.z()) < 2000.")\
                 .Filter("abs(pdg) == {}".format(pdg))\
                 .Define("mom", "ts_calo_mom.r()")\
                 .Define("beta", beta)\
                 .Histo2D(("h", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")

            h.FitSlicesY()
            h_1 = ROOT.gROOT.FindObject("h_1")
            h_2 = ROOT.gROOT.FindObject("h_2")
            h_chi2 = ROOT.gROOT.FindObject("h_chi2")

            for p in range(1, h_1.GetNbinsX()+1):
                gr.SetPoint(p, h_1.GetBinCenter(i), h_1.GetBinContent(i))
                gr.SetPointError(p, 0., h_2.GetBinContent(i))
            graphs[i][j][pdg] = gr
        # sep_power = abs(m1-m2) / np.sqrt( (std1*std1 + std2*std2)/2. )



graphs[10]["set"][211].Draw("AP")
graphs[10]["set"][321].Draw("Psame")
graphs[10]["set"][2212].Draw("Asame")


canvas.SetGridx()
canvas.SetGridy()
canvas.Update()
input("wait")
