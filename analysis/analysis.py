import ROOT
import time
import cProfile, pstats, io
from pstats import SortKey
import numpy as np

# This all for 3D plots
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt

# calib_kaons data
# df = ROOT.RDataFrame("ana_tree", "/nfs/dust/ilc/user/dudarboh/final_files/calib_kaons.root")

# 2f_Z_hadronic data
df = ROOT.RDataFrame("ana_tree", "/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result[1-3].root")

def tof_analysis():
    # Include all TOF algorithms written on c++, Hit struct which helps and speed of light const
    ROOT.gInterpreter.Declare('#include "analysis.h"')


    # test =  df.Filter("nCalHits > 0 && abs(d0) < 10 && abs(z0) < 20 && p < 2.5").Define("red_chi2", "chi2/ndf").Histo1D(("test","test", 500, 0., 5.), "red_chi2")
    # test.Draw()
    # input("wait")

    # All basic background rejecting cuts
    df1 = df.Filter("nCalHits > 0 && abs(d0) < 10 && abs(z0) < 20 && p < 2.5 && chi2/ndf < 1.01")\
    .Define("layer_cut", "layerCalHit < 10")\
    .Define("l", "layerCalHit[layer_cut] ")\
    .Filter("l.size() > 0")\
    .Define("t", "tCalHit[layer_cut]")\
    .Define("r", "dToRefPointCalHit[layer_cut]")\
    .Define("d", "dToLineCalHit[layer_cut]")\
    .Define("tof", "tof_fit(t, r, d, l, 0)")\
    .Define("betaCalo", "tof_fit(t, r, d, l, 1)")\
    .Define("lenIP", "abs((phi-phiCalState)/omega)*sqrt(1. + tanL*tanL)")\
    .Define("lenCalo", "abs((phi-phiCalState)/omegaCalState)*sqrt(1. + tanLCalState*tanLCalState)")\
    .Define("lenAVG", "(lenIP+lenCalo)/2.")\

    h1 = df1.Define("beta", "lenCalo/(tof * c)")\
    .Define("mass", "p / beta * sqrt(1. - beta * beta) * 1000")\
    .Histo1D(("h1","CalState Length", 2000, 0., 1000.), "mass")

    h2 = df1.Define("beta", "lenAVG/(tof * c)")\
    .Define("mass", "p / beta * sqrt(1. - beta * beta) * 1000")\
    .Histo1D(("h2","Avg len", 2000, 0., 1000.), "mass")

    h3 = df1.Define("beta", "lenCalo/(tof * c)").Define("betaIP", "lenIP/(tof * c)")\
    .Define("mass", " (p/betaIP -dEdX*lenCalo) * sqrt(1. - beta * beta) * 1000")\
    .Histo1D(("h3"," p/betaIP", 2000, 0., 1000.), "mass")

    h4 = df1.Define("beta", "lenCalo/(tof * c)").Define("betaIP", "lenIP/(tof * c)")\
    .Define("mass", " (p/betaIP -dEdX*lenIP) * sqrt(1. - beta * beta) * 1000")\
    .Histo1D(("h4"," p/betaIP len IP", 2000, 0., 1000.), "mass")


    # h1.Draw()
    # input("wait")
    # h2.Draw()
    # input("wait")
    # h3.Draw()
    # input("wait")
    # h4.Draw()
    # input("wait")

    histos = [h1, h2, h3, h4]


    canvas = ROOT.TCanvas()
    hs = ROOT.THStack()

    for i, h in enumerate(histos):
        hs.Add(h.GetPtr())
        h.SetLineColor(i+1)

    hs.Draw("nostack")

    hs.SetTitle("Different lengths;mass, [MeV];N_{PFO}")
    canvas.BuildLegend()
    canvas.SetGridx()
    canvas.SetGridy()
    # Kaon line
    line1 = ROOT.TLine(493.677, 0., 493.677, 2500.)
    line1.SetLineWidth(3)
    line1.SetLineColor(6)
    line1.Draw()
    # Pion line
    line2 = ROOT.TLine(139.570, 0., 139.570, 6800.)
    line2.SetLineWidth(3)
    line2.SetLineColor(6)
    line2.Draw()

    # Proton line
    line3 = ROOT.TLine(938.272, 0., 938.272, 2000.)
    line3.SetLineWidth(3)
    line3.SetLineColor(6)
    line3.Draw()

    canvas.Update()
    input("wait")


def test_bias():
    canvas = ROOT.TCanvas()
    x = np.array([139.570183535, 493.677, 938.2720881629])
    mg = ROOT.TMultiGraph()

    # IP Fit
    # y = (np.array([113.066, 483.685, 9.33090e+02]) - x)
    # IP Avg
    # y = (np.array([1.18410e+02, 4.90807e+02, 9.44576e+02]) - x)
    # IP Closest
    # y = (np.array([114.629, 4.85649e+02, 9.36287e+02]) - x)

    # Fit
    y1_prev = (np.array([1.44045e+02, 4.96367e+02, 9.42326e+02]) - x)
    y1 = (np.array([1.45008e+02, 4.96367e+02, 9.42326e+02]) - x)

    gr_fit2 = ROOT.TGraphErrors(3, x, y)
    gr_fit2.SetTitle("TOF Fit (Calo)")
    gr_fit2.SetMarkerStyle(22)
    gr_fit2.SetMarkerColor(1)
    gr_fit2.SetLineColor(1)
    gr_fit2.SetLineStyle(10)
    mg.Add(gr_fit2)

    # Avg
    y = (np.array([1.48365e+02, 5.03939e+02, 9.54064e+02]) - x)
    gr_avg2 = ROOT.TGraphErrors(3, x, y)
    gr_avg2.SetTitle("TOF avg (Calo)")
    gr_avg2.SetMarkerStyle(22)
    gr_avg2.SetMarkerColor(2)
    gr_avg2.SetLineColor(2)
    gr_avg2.SetLineStyle(10)
    mg.Add(gr_avg2)

    # Closest
    y = (np.array([1.45594e+02, 4.98747e+02, 9.45545e+02]) - x)
    gr_closest2 = ROOT.TGraphErrors(3, x, y)
    gr_closest2.SetTitle("TOF closest (Calo)")
    gr_closest2.SetMarkerStyle(22)
    gr_closest2.SetMarkerColor(3)
    gr_closest2.SetLineColor(3)
    gr_closest2.SetLineStyle(10)
    mg.Add(gr_closest2)


    mg.Draw("APL")
    canvas.BuildLegend()
    mg.SetTitle("Bias vs mass;m_{true}, [MeV]; m_{reco}-m_{true}, [MeV]")
    canvas.Update()
    input("wait")


pr = cProfile.Profile()
pr.enable()

tof_analysis()
# test_bias()




pr.disable()
ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
ps.print_stats()
