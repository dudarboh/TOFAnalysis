import ROOT

ROOT.gInterpreter.Declare('''

#include "TRandom3.h"

TRandom3 r;
#define SPEED_OF_LIGHT 299.792458
''')


ch = ROOT.TChain("SETAnalysis")

ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")

df = ROOT.RDataFrame(ch).Filter("n_ecal_hits > 0 && abs(ts_ecal_z0) < 1.")
df = df.Filter("abs(ts_ecal_pos.z()) > 2385.")
df = df.Define("mom", "ts_ecal_mom.r()")

canvas = ROOT.TCanvas()
canvas.SetGridx()
canvas.SetGridy()
ROOT.gStyle.SetMarkerStyle(0)
ROOT.gStyle.SetOptStat(10)

# hz = df.Histo1D(("hz", "TrackAtEcal.z0", 6000, -3000, 3000), "ts_ecal_z0")
# hz.Draw()
# input("wait")

# n_bad = df.Filter("track_length_ip < 0.").Count().GetValue()
# n_total = df.Count().GetValue()
#
# print("N phi fliped pfos: ", n_bad / n_total * 100, " %    ", n_bad, "/", n_total)


h1 = df.Define("dif", "track_length_ip - track_length_refit_z").Histo1D(("h1", "TDR; track length (mm); N pfo", 1000, -5000, 5000), "dif")
h1.Draw()
h2 = df.Define("dif", "track_length_calo - track_length_refit_z").Histo1D(("h2", "Dev", 1000, -5000, 5000), "dif")
h2.SetLineColor(2)
h2.Draw("same")
# h3 = df.Histo1D(("h3", "Refit1", 5200, -10000, 10000), "track_length_refit_tanL")
# h3.SetLineColor(4)
# h3.Draw("same")
# h4 = df.Histo1D(("h4", "Refit2", 5200, -10000, 10000), "track_length_refit_z")
# h4.SetLineColor(8)
# h4.Draw("same")





canvas.BuildLegend()
h1.SetTitle("Endcap ( |z| > 2385 mm )")
canvas.Update()
input("wait")
