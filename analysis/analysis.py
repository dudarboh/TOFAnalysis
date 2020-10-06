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

# ROOT.EnableImplicitMT();


# c in mm/ns
c = 299.792458

# calib_kaons data
# df = ROOT.RDataFrame("ana_tree", "/nfs/dust/ilc/user/dudarboh/final_files/calib_kaons.root")

# 2f_Z_hadronic data
df = ROOT.RDataFrame("ana_tree", "/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic_new/*.root")

def tof_analysis():
    # Include all TOF algorithms written on c++, Hit struct which helps and speed of light const
    ROOT.gInterpreter.Declare('#include "analysis.h"')

    filtered_df = df.Filter("nCalHits > 0 && abs(d0) < 10 && abs(z0) < 20 && p < 2.5")\
    .Define("layerNewCalHit", "layerCalHit[layerCalHit < 10] ").Filter("layerNewCalHit.size() > 0")\
    .Define("tNewCalHit", "tCalHit[layerCalHit < 10]")\
    .Define("dToRefPointNewCalHit", "dToRefPointCalHit[layerCalHit < 10]")\
    .Define("dToLineNewCalHit", "dToLineCalHit[layerCalHit < 10]")\
    .Define("length_IP", "abs((-phiCalState+phi)/omega)*sqrt(1. + tanL*tanL)")

    df_closest = filtered_df.Define("tof", "tof_closest(tNewCalHit, dToRefPointNewCalHit)")\
    .Define("beta", "length/(tof * c)")\
    .Define("mass", "p / beta * sqrt(1. - beta * beta) * 1000")

    df_avg = filtered_df.Define("tof", "tof_avg(tNewCalHit, dToRefPointNewCalHit, dToLineNewCalHit, layerNewCalHit)")\
    .Define("beta", "length/(tof * c)")\
    .Define("mass", "p / beta * sqrt(1. - beta * beta) * 1000")

    df_fit = filtered_df.Define("tof", "tof_fit(tNewCalHit, dToRefPointNewCalHit, dToLineNewCalHit, layerNewCalHit)")\
    .Define("beta", "length/(tof * c)")\
    .Define("mass", "p / beta * sqrt(1. - beta * beta) * 1000")

    h1 = df_fit.Histo1D(("h1","TOF fit", 1000, 0., 1000.), "mass")
    h1.SetLineColor(1)
    h2 = df_avg.Histo1D(("h2","TOF avg", 1000, 0., 1000.), "mass")
    h2.SetLineColor(2)
    h3 = df_closest.Histo1D(("h3","TOF closest", 1000, 0., 1000.), "mass")
    h3.SetLineColor(3)

    canvas = ROOT.TCanvas()
    hs = ROOT.THStack()
    hs.Add(h1.GetPtr())
    hs.Add(h2.GetPtr())
    hs.Add(h3.GetPtr())
    hs.Draw("nostack")
    hs.SetTitle("TOF estimators (#Omega_{IP}, #lambda_{IP});mass, [MeV];N_{PFO}")
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
    #Fit algo
    canvas = ROOT.TCanvas()
    x = np.array([139.570183535, 493.677, 938.2720881629])
    mg = ROOT.TMultiGraph()

    y = (np.array([113.066, 483.685, 9.33090e+02]) - x)
    gr_fit = ROOT.TGraphErrors(3, x, y)
    gr_fit.SetTitle("TOF Fit (IP)")
    gr_fit.SetMarkerColor(1)
    gr_fit.SetLineColor(1)
    mg.Add(gr_fit)

    # Avg algo
    y = (np.array([1.18410e+02, 4.90807e+02, 9.44576e+02]) - x)
    gr_avg = ROOT.TGraphErrors(3, x, y)
    gr_avg.SetTitle("TOF avg (IP)")
    gr_avg.SetMarkerColor(2)
    gr_avg.SetLineColor(2)
    mg.Add(gr_avg)

    # Closest algo
    y = (np.array([114.629, 4.85649e+02, 9.36287e+02]) - x)
    gr_closest = ROOT.TGraphErrors(3, x, y)
    gr_closest.SetTitle("TOF closest (IP)")
    gr_closest.SetMarkerColor(3)
    gr_closest.SetLineColor(3)
    mg.Add(gr_closest)

    # Fit
    y = (np.array([1.44045e+02, 4.96367e+02, 9.42326e+02]) - x)
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

# tof_analysis()
test_bias()




pr.disable()
ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
ps.print_stats()
