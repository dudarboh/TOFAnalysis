import ROOT

# ROOT.gInterpreter.Declare('#include "extract.hpp"')

ROOT.gInterpreter.Declare('''

#include "TRandom3.h"

TRandom3 r;

''')


canvas = ROOT.TCanvas()

ch = ROOT.TChain("SETAnalysis")

ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/SET/final1.root")

df = ROOT.RDataFrame(ch)
df1 = df.Filter("ts_last_mom.rho() > 5.")\
       .Define("theta", "ts_last_mom.theta()")\
       .Define("z", "ts_last_pos.z()")\
       .Define("n_hits", "n_set_hits + r.Gaus(0, 0.05)")

gr = df1.Graph("z", "n_hits")
# gr.Draw("AP")
gr.SetMarkerStyle(1)

h2 = df.Filter("abs(ts_calo_pos.z()) < 2200")\
        .Histo1D(("h", "h", 2, 0, 2),"n_set_hits")
h2.Draw()
# h1 = df.Histo1D(("h1", "All SimHits; z (mm); n hits", 4600, -2300, 2300), "z")
# h2 = df.Filter("is_passed != 0").Histo1D(("h2", "Matched to SpacePoint", 4600, -2300, 2300), "z")
#
#
# h1.Draw()
# h2.SetLineColor(2)
# h2.Draw("same")
canvas.Update()
input("wait")
