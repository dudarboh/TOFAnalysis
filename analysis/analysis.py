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
ch = ROOT.TChain("PFO")
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result1.root")

ch1 = ROOT.TChain("Cluster")
ch1.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result1.root")

ch2 = ROOT.TChain("PionTrack")
ch2.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result1.root")

ch3 = ROOT.TChain("TrackerHits")
ch3.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result1.root")

ch4 = ROOT.TChain("ECALHits")
ch4.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result1.root")

ch.AddFriend(ch1, "cluster")
ch.AddFriend(ch2, "piFit")
ch.AddFriend(ch3, "tr")
ch.AddFriend(ch4, "cal")

d = ROOT.RDataFrame(ch)
ROOT.gInterpreter.Declare('#include "analysis.h"')

def photon_tof():
    # Cuts on nice photons
    d1 = d.Filter("cal.nHits > 10 \
                  && PDG[0] == 22\
                  && nMC == 1\
                  && isBackscatter[0] == 0\
                  && isDecayedInTracker[0] == 0\
                  && abs(cluster.z) < 2200.\
                  && abs(xMC[0]) < 1.\
                  && abs(yMC[0]) < 1.\
                  && abs(zMC[0]) < 1.")

    # Cuts for nice hits in calo for fit
    d2 = d1.Define("caloHitCut", "cal.layer < 10 && cal.t > 1.e-3")\
    .Define("layerCal", "cal.layer[caloHitCut]").Filter("layerCal.size() > 0")\
    .Define("xCal", "cal.x[caloHitCut]")\
    .Define("yCal", "cal.y[caloHitCut]")\
    .Define("zCal", "cal.z[caloHitCut]")\
    .Define("tCal", "cal.t[caloHitCut]")

    # Line parameters
    d3 = d2.Define("pxLine", "cluster.x - xMC[0]")\
    .Define("pyLine", "cluster.y - yMC[0]")\
    .Define("pzLine", "cluster.z - zMC[0]")



    # Plane parameters
    d4 = d3.Define("phiLine", "atan2(pyLine,pxLine)").Define("thetaLine", "atan2( sqrt(pxLine*pxLine + pyLine*pyLine), pzLine)")\
    .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")

    d5 = d4.Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, pxLine, pyLine, pzLine, cluster.x, cluster.y, cluster.z)")

    # Impact parameters
    d6 = d5.Define("xImpact", "tLineParam*pxLine + cluster.x")\
        .Define("yImpact", "tLineParam*pyLine + cluster.y")\
        .Define("zImpact", "tLineParam*pzLine + cluster.z")\

    # TOF from 10 closest hits
    d7 = d6.Define("dToLine", "dToLine(xCal, yCal, zCal, xImpact, yImpact, zImpact, phiLine, thetaLine)")\
    .Define("dToRef", "dToRef(xCal, yCal, zCal, xImpact, yImpact, zImpact)")\
    .Define("tofFit", "tof_fit(tCal, dToRef, dToLine, layerCal)")\

    # TOF from track length to impact point
    d8 = d7.Define("tofLen", "sqrt((xImpact-xMC[0])*(xImpact-xMC[0]) + (yImpact-yMC[0])*(yImpact-yMC[0]) + (zImpact-zMC[0])*(zImpact-zMC[0]))/c")
    # d5 = d4.Define("tofLen", "sqrt((xImpact)*(xImpact) + (yImpact)*(yImpact) + (zImpact)*(zImpact))/c")

    h1 = d8.Define("diff", "tofFit - tofLen").Histo1D(("h1","TOF_{fit} - TOF_{len} (Cluster pos)", 500, -.25, .25), "diff")

    # Line parameters
    d3 = d2.Define("rCluster", "sqrt(cluster.x*cluster.x + cluster.y*cluster.y + cluster.z*cluster.z)")\
        .Define("pxLine", "getLineDir(rCluster, cluster.phi, cluster.theta, 0)")\
        .Define("pyLine", "getLineDir(rCluster, cluster.phi, cluster.theta, 1)")\
        .Define("pzLine", "getLineDir(rCluster, cluster.phi, cluster.theta, 2)")\

    # Plane parameters
    d4 = d3.Define("nx", "getPlaneDir(cluster.phi, 0)").Define("ny", "getPlaneDir(cluster.phi, 1)").Define("nz", "getPlaneDir(cluster.phi, 2)")

    d5 = d4.Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, pxLine, pyLine, pzLine, cluster.x, cluster.y, cluster.z)")

    # Impact parameters
    d6 = d5.Define("xImpact", "tLineParam*pxLine + cluster.x")\
        .Define("yImpact", "tLineParam*pyLine + cluster.y")\
        .Define("zImpact", "tLineParam*pzLine + cluster.z")\

    # TOF from 10 closest hits
    d7 = d6.Define("dToLine", "dToLine(xCal, yCal, zCal, xImpact, yImpact, zImpact, cluster.phi, cluster.theta)")\
    .Define("dToRef", "dToRef(xCal, yCal, zCal, xImpact, yImpact, zImpact)")\
    .Define("tofFit", "tof_fit(tCal, dToRef, dToLine, layerCal)")\

    # TOF from track length to impact point
    d8 = d7.Define("tofLen", "sqrt((xImpact-xMC[0])*(xImpact-xMC[0]) + (yImpact-yMC[0])*(yImpact-yMC[0]) + (zImpact-zMC[0])*(zImpact-zMC[0]))/c")
    # d5 = d4.Define("tofLen", "sqrt((xImpact)*(xImpact) + (yImpact)*(yImpact) + (zImpact)*(zImpact))/c")

    h2 = d8.Define("diff", "tofFit - tofLen").Histo1D(("h2","TOF_{fit} - TOF_{len} (Cluster dir)", 500, -.25, .25), "diff")

    h1.Draw()
    h2.Draw("same")
    input("wait")

def tof_analysis():
    c = ROOT.TCanvas()
    # Good PFO cuts
    d1 = d.Filter("cal.nHits > 0")

    # Cut on calorimeter layer and time
    d2 = d1.Define("caloHitCut", "cal.layer < 10 && cal.t > 1.e-3")\
    .Define("layerCal", "cal.layer[caloHitCut]").Filter("cal.layer.size() > 0")\
    .Define("xCal", "cal.x[caloHitCut]")\
    .Define("yCal", "cal.y[caloHitCut]")\
    .Define("zCal", "cal.z[caloHitCut]")\
    .Define("tCal", "cal.t[caloHitCut]")
    # Define calorimeter hit parameters
    # d3 = d2.Define("dToLine", "dToLine(xCal, yCal, zCal, pFit.xRefCalState, pFit.yRefCalState, pFit.zRefCalState, pFit.phiCalState, atan(1./pFit.tanLCalState))")\
    # .Define("dToRef", "dToRef(xCal, yCal, zCal, pFit.xRefCalState, pFit.yRefCalState, pFit.zRefCalState)")\
    # .Define("tof", "tof_fit(tCal, dToRef, dToLine, layerCal)")\
    # .Define("length", "abs((pFit.phi-pFit.phiCalState)/pFit.omegaCalState)*sqrt(1. + pFit.tanLCalState*pFit.tanLCalState)")\
    # .Define("beta", "length/(tof * c)")\
    #
    # h2 = d3.Define("mass", "pFit.pCalState / beta * sqrt(1. - beta * beta) * 1000")\
    # .Histo1D(("h2","Calo momentum", 2000, 0., 1000.), "mass")

    h1 = d2.Define("diff", "piFit.tanL - piFit.tanLFirstHit")\
    .Histo1D(("h1","tan(#lambda_{IP}) - tan(#lambda_{first})", 2000, -0.001, 0.001), "diff")

    h1.Draw()
    # h2.Draw("same")
    # h2.SetLineColor(2)

    # h4.Draw()
    # input("wait")

    # hs = ROOT.THStack()
    # hs.Add(h1.GetPtr())
    # hs.Add(h2.GetPtr())
    # hs.Draw("nostack")
    c.BuildLegend()
    c.Update()
    input("wait")


def test_bias():
    canvas = ROOT.TCanvas()
    x = np.array([139.570183535, 493.677, 938.2720881629])
    x_err = np.array([0., 0.013, 0.])
    mg = ROOT.TMultiGraph()
    mg.SetTitle("Bias vs mass (fit mass corrected);m_{true}, [MeV]; m_{reco}-m_{true}, [MeV]")

    # IP Fit OLD and wrong results
    # y = (np.array([113.066, 483.685, 9.33090e+02]) - x)
    # IP Avg
    # y = (np.array([1.18410e+02, 4.90807e+02, 9.44576e+02]) - x)
    # IP Closest
    # y = (np.array([114.629, 4.85649e+02, 9.36287e+02]) - x)

    # P from IP. Reference. NO nSETHits cut
    # y = np.array([1.44266e+02, 4.96488e+02, 9.42166e+02]) - x
    # y = np.array([1.44266e+02, 4.96488e+02, 9.42166e+02]) - x
    # y_err = np.array([1.05738e-01, 1.15293e-01, 9.78054e-02])
    # gr00 = ROOT.TGraphErrors(3, x, y, x_err, y_err)
    # gr00.SetTitle("p_{IP} w/o N_{SET} cut")
    # gr00.SetMarkerStyle(20)
    # gr00.SetLineStyle(10)
    # gr00.SetMarkerColor(1)
    # gr00.SetLineColor(1)
    # mg.Add(gr00)
    # mg.Draw("APL")

    # P from Calo. NO nSETHits cut
    # pion fit
    y = np.array([1.43536e+02, 4.94994e+02, 9.39755e+02]) - x
    # kaon fit
    y = np.array([1.43682e+02, 4.94791e+02, 9.39536e+02]) - x
    # proton fit
    y = np.array([1.45416e+02, 4.95307e+02, 9.39818e+02]) - x
    #mixed fit
    y = np.array([1.43536e+02, 4.94791e+02, 9.39818e+02]) - x

    y_err = np.array([1.03603e-01, 1.22941e-01, 1.00797e-01])
    gr11 = ROOT.TGraphErrors(3, x, y, x_err, y_err)
    gr11.SetTitle("p_{calo} w/o N_{SET} cut")
    gr11.SetMarkerStyle(20)
    gr11.SetLineStyle(10)
    gr11.SetMarkerColor(2)
    gr11.SetLineColor(2)
    mg.Add(gr11)
    mg.Draw("APL")

    mg.GetXaxis().SetNdivisions(510)
    mg.GetYaxis().SetNdivisions(520)
    mg.SetMaximum(10.)
    mg.SetMinimum(-11.)

    canvas.Update()
    input("wait")

    # pion fit
    y = np.array([1.29531e+02, 4.87609e+02, 9.33610e+02]) - x
    # kaon fit
    y = np.array([1.29719e+02, 4.87284e+02, 9.33730e+02]) - x
    # proton fit
    y = np.array([1.30076e+02, 4.87430e+02, 9.33759e+02]) - x
    #mixed fit
    y = np.array([1.29531e+02, 4.87284e+02, 9.33759e+02]) - x

    y_err = np.array([1.70929e-01, 2.18768e-01, 1.35535e-01])
    gr3 = ROOT.TGraphErrors(3, x, y, x_err, y_err)
    gr3.SetTitle("SET")
    gr3.SetMarkerStyle(20)
    gr3.SetLineStyle(1)
    gr3.SetMarkerColor(4)
    gr3.SetLineColor(4)
    mg.Add(gr3)
    mg.Draw("APL")
    canvas.Update()
    input("wait")

    # P from Calo
    # pion fit
    y = np.array([1.47467e+02, 4.94567e+02, 9.39285e+02]) - x
    # kaon fit
    y = np.array([1.47637e+02, 4.94406e+02, 9.39199e+02]) - x
    # proton fit
    y = np.array([1.48266e+02, 4.94793e+02, 9.39285e+02]) - x
    #mixed fit
    y = np.array([1.47467e+02, 4.94406e+02, 9.39285e+02]) - x

    y_err = np.array([1.63851e-01, 1.49119e-01, 1.12972e-01])
    gr1 = ROOT.TGraphErrors(3, x, y, x_err, y_err)
    gr1.SetTitle("p_{calo} w/ N_{SET} cut")
    gr1.SetMarkerStyle(20)
    gr1.SetLineStyle(1)
    gr1.SetMarkerColor(2)
    gr1.SetLineColor(2)
    mg.Add(gr1)

    mg.Draw("APL")
    canvas.BuildLegend()
    canvas.Update()
    input("wait")


pr = cProfile.Profile()
pr.enable()

photon_tof()
# tof_analysis()
# test_bias()


pr.disable()
ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
ps.print_stats()
