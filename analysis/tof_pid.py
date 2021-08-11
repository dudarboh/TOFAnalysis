import ROOT
import numpy as np
ROOT.EnableImplicitMT(2)

ROOT.gInterpreter.Declare('''

#define SPEED_OF_LIGHT 299.792458

''')

df_0 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/0ps.root")
df_10 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/10ps.root")
df_30 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/30ps.root")
df_50 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/50ps.root")
df_100 = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/100ps.root")

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
canvas.Divide(2, 1)
canvas.SetGridx()
canvas.SetGridy()
graphs = {}

for i, df in zip(df_names, dfs):
    for j, beta in zip(beta_names, betas):

        h_all = df.Filter("n_ecal_hits > 0 && abs(pos_fastest.z()) < 2000.")\
             .Define("mom", "ts_calo_mom.r()")\
             .Define("beta", beta)\
             .Histo2D(("h_all", "Method: {} {} ps; p (GeV);#beta".format(j, i), 30, 1., 10., 2000, 0.7, 1.05), "mom", "beta")

        canvas.cd(1)
        ROOT.gPad.SetLogz()
        h_all.Draw("colz")
        h_all.SetStats(0)

        for k, pdg in enumerate(pdgs):
            gr = ROOT.TGraphErrors()
            gr.SetMarkerStyle(20)
            gr.SetMarkerColor(k+1)
            gr.SetLineColor(k+1)

            h = df.Filter("n_ecal_hits > 0 && abs(pos_fastest.z()) < 2000.")\
                 .Filter("abs(pdg) == {}".format(pdg))\
                 .Define("mom", "ts_calo_mom.r()")\
                 .Define("beta", beta)\
                 .Histo2D(("h", "title;mom;beta", 30, 1., 10., 1000, 0.7, 1.05), "mom", "beta")

            h.FitSlicesY()
            h_1 = ROOT.gROOT.FindObject("h_1")
            h_2 = ROOT.gROOT.FindObject("h_2")
            h_chi2 = ROOT.gROOT.FindObject("h_chi2")

            for p in range(1, h_1.GetNbinsX()+1):
                gr.SetPoint(p-1, h_1.GetBinCenter(p), h_1.GetBinContent(p))
                gr.SetPointError(p-1, 0., h_2.GetBinContent(p))
            graphs[i, j, pdg] = gr

            canvas.cd(2)
            gr.Draw("APE" if pdg == 211 else "PEsame")
            gr.GetXaxis().SetRangeUser(1., 10.)
            gr.GetYaxis().SetRangeUser(0.7, 1.05)

        canvas.Update()
        canvas.Print("method_{}_{}_ps.png".format(j, i))
        input("wait")
        # sep_power = abs(m1-m2) / np.sqrt( (std1*std1 + std2*std2)/2. )



graphs[10, "set", 211].Draw("AP")
graphs[10, "set", 321].Draw("Psame")
graphs[10, "set", 2212].Draw("Psame")


canvas.SetGridx()
canvas.SetGridy()
canvas.Update()
input("wait")
