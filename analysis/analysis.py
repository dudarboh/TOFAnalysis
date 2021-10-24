import ROOT

ROOT.gStyle.SetOptStat(10)
ROOT.gStyle.SetMarkerStyle(0)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)
ROOT.gStyle.SetNdivisions(512)
ROOT.gStyle.SetNdivisions(512, "Y")
canvas = ROOT.TCanvas()

ROOT.gInterpreter.Declare('''

#include "TRandom3.h"

TRandom3 r;
#define SPEED_OF_LIGHT 299.792458
''')


ch = ROOT.TChain("TOFAnalysis")

# ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
ch.Add("/afs/desy.de/user/d/dudarboh/TOFAnalysis/build/initial.root")
ch.AddFriend("fix=TOFAnalysis", "/afs/desy.de/user/d/dudarboh/TOFAnalysis/build/fixed.root")


df = ROOT.RDataFrame(ch)

h_diff = df.Define("dr_old", "(ts_ecal_pos - pos_closest).r()")\
          .Define("dr_new", "(fix.ts_ecal_pos - fix.pos_closest).r()")\
          .Define("dr_diff", "dr_new - dr_old")\
          .Histo1D(("h_diff", "DIFF;#Delta|#vec{r}_{track, ecal} - #vec{r}_{closest}| (mm); N PFO", 1000, -10000, 10000), "dr_diff")

# h_new = df.Histo1D(("h_new", "NEW;|#vec{r}_{track, ecal} - #vec{r}_{closest}| (mm); N PFO", 1000, 0, 10000), "dr_new")

h_diff.Draw()
# h_new.Draw("same")
# h_new.SetLineColor(ROOT.kRed + 1)


# h1 = df.Filter("std::abs( ts_ecal_pos.z() ) < 2385.").Histo1D(("h1", "Barrel;N ECal hits;N PFO", 150, 0, 150,), "n_ecal_hits" )
# h2 = df.Filter("std::abs( ts_ecal_pos.z() ) >= 2385.").Histo1D(("h2", "Endcap;N ECal hits;N PFO", 150, 0, 150,), "n_ecal_hits" )
# h1 = df.Filter("std::abs( ts_ecal_pos.z() ) < 2385.").Histo1D(("h1", "Barrel;Track length (mm);N PFO", 2500, 0, 20000,), "track_length_ecal" )
# h2 = df.Filter("std::abs( ts_ecal_pos.z() ) >= 2385.").Histo1D(("h2", "Endcap;Track length (mm);N PFO", 2500, 0, 20000,), "track_length_ecal" )

# h1 = df.Define("dr", "(ts_ecal_pos - pos_closest).r()").Filter("std::abs( ts_ecal_pos.z() ) < 2385.").Histo1D(("h1", "Barrel;|#vec{r}_{track, ecal} - #vec{r}_{closest}| (mm);N PFO", 1000, 0, 10000,), "dr" )
# h2 = df.Define("dr", "(ts_ecal_pos - pos_closest).r()").Filter("std::abs( ts_ecal_pos.z() ) >= 2385.").Histo1D(("h2", "Endcap;|#vec{r}_{track, ecal} - #vec{r}_{closest}| (mm);N PFO", 1000, 0, 10000,), "dr" )

# h1 = df.Histo1D(("h1", "Barrel;mom_hm_ecal (GeV);N PFO", 1000, 0, 100,), "mom_hm_ecal" )
# h2 = df.Define("mom", "ts_ecal_mom.r()").Histo1D(("h2", "Endcap;p ecal (GeV);N PFO", 1000, 0, 100,), "mom" )


# h2.SetLineColor(2)
# h1.Draw()
# h2.Draw("sames")

canvas.BuildLegend()
h_diff.SetTitle("")
canvas.Update()
input("wait")





df = df.Filter("n_ecal_hits > 0 && abs(ts_ecal_z0) < 1.")
df = df.Filter("abs(ts_ecal_pos.z()) > 2385.")
df = df.Define("mom", "ts_ecal_mom.r()")

canvas = ROOT.TCanvas()
canvas.SetGridx()
canvas.SetGridy()
ROOT.gStyle.SetMarkerStyle(0)

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
