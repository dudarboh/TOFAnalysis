import ROOT

# ROOT.gInterpreter.Declare('#include "extract.hpp"')

ROOT.gInterpreter.Declare('''

#include "TRandom3.h"

TRandom3 r;
#define SPEED_OF_LIGHT 299.792458
''')


canvas = ROOT.TCanvas()
canvas.SetGridx()

ch = ROOT.TChain("SETAnalysis")

ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/SET/250gev_prod.root")

df = ROOT.RDataFrame(ch).Filter("n_ecal_hits > 0 && abs(pos_fastest.z()) < 2000.")\
                        .Define("mom", "ts_calo_mom.r()")\

h_set = df.Define("beta", "track_length_set/(set_hit_time*SPEED_OF_LIGHT)")\
           .Histo2D(("h_set", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_set.SetStats(0)
h_set.SetMaximum(200)

h_closest = df.Define("beta", "track_length_calo/(tof_closest*SPEED_OF_LIGHT)")\
           .Histo2D(("h_closest", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_closest.SetStats(0)
h_closest.SetMaximum(200)

h_closest2 = df.Define("beta", "track_length_calo/((tof_closest+ (pos_closest - ts_calo_pos).r()/SPEED_OF_LIGHT )*SPEED_OF_LIGHT)")\
           .Histo2D(("h_closest2", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_closest2.SetStats(0)
h_closest2.SetMaximum(200)


h_closest_sim = df.Define("beta", "track_length_calo/(tof_closest_sim*SPEED_OF_LIGHT)")\
           .Histo2D(("h_closest_sim", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_closest_sim.SetStats(0)
h_closest_sim.SetMaximum(200)

h_fastest = df.Define("beta", "track_length_calo/(tof_fastest*SPEED_OF_LIGHT)")\
           .Histo2D(("h_fastest", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_fastest.SetStats(0)
h_fastest.SetMaximum(200)

h_avg = df.Define("beta", "track_length_calo/(tof_frank_avg*SPEED_OF_LIGHT)")\
           .Histo2D(("h_avg", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_avg.SetStats(0)
h_avg.SetMaximum(200)

h_fit = df.Define("beta", "track_length_calo/(tof_frank_fit*SPEED_OF_LIGHT)")\
           .Histo2D(("h_fit", "title;mom;beta", 200, 1., 10., 500, 0.95, 1.01), "mom", "beta")
h_fit.SetStats(0)
h_fit.SetMaximum(200)


h_set.Draw("colz")
canvas.Update()
# input("wait")


h_closest.Draw("colz")
canvas.Update()
# input("wait")

h_closest2.Draw("colz")
canvas.Update()
# input("wait")


h_closest_sim.Draw("colz")
canvas.Update()
# input("wait")

h_fastest.Draw("colz")
canvas.Update()
input("wait")

h_avg.Draw("colz")
canvas.Update()
input("wait")

h_fit.Draw("colz")
canvas.Update()
input("wait")
