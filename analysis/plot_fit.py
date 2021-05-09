import ROOT
ROOT.gInterpreter.Declare('#include "extract.hpp"')

ch = ROOT.TChain("TOFAnalysis")
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")



smearing = 0.



df = ROOT.RDataFrame(ch)
df = df.Filter("nECALHits > 0 && abs(xyzCluster.Z()) < 2200. && abs(xyzVtxMC.R()) < .5").Range(1)
df = df.Define("tofHit", "tofHit(tECALHit, {})".format(smearing))
df = df.Define("dToImpact", "dToImpact(xyzECALHit, xyzTrackAtCalo)")\
       .Define("dToLine", "dToLine(xyzECALHit, xyzTrackAtCalo, pTrackAtCalo)")\
       .Define("sel_frank", "selectHits(dToLine, layerECALHit, true, 10, 9999.)")\
       .Define("plot_fit", "plot_fit(tofHit[sel_frank], dToImpact[sel_frank])")

h = df.Histo1D("plot_fit")
h.Draw()
