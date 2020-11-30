import ROOT
import time
import numpy as np

ROOT.gStyle.SetPalette(1)

# 2f_Z_hadronic data
ch = ROOT.TChain("PFO")
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch1 = ROOT.TChain("Cluster")
ch1.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch2 = ROOT.TChain("PionTrack")
ch2.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch3 = ROOT.TChain("TrackerHits")
ch3.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch4 = ROOT.TChain("ECALHits")
ch4.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

ch.AddFriend(ch1, "cluster")
ch.AddFriend(ch2, "piFit")
ch.AddFriend(ch3, "tr")
ch.AddFriend(ch4, "cal")

df_initial = ROOT.RDataFrame(ch)
ROOT.gInterpreter.Declare('#include "analysis.hpp"')

def good_photons(df):
    res = df.Filter("nMC == 1 && PDG[0] == 22 && cal.nHits > 0 \
    && isBackscatter[0] == 0 && isDecayedInTracker[0] == 0\
    && abs(cluster.z) < 2200. && abs(xMC[0]) < .5\
    && abs(yMC[0]) < .5 && abs(zMC[0]) < .5")
    return res;

def use_mc(df):
    # Calculations of track length and direction and intersection based on true info
    res = df.Define("phiLine", "atan2(pyMC[0], pxMC[0])").Define("thetaLine", "atan2( sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0]), pzMC[0])")\
    .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")\
    .Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, pxMC[0], pyMC[0], pzMC[0], xMC[0], yMC[0], zMC[0])")\
    .Define("xImpact", "tLineParam*pxMC[0] + xMC[0]")\
    .Define("yImpact", "tLineParam*pyMC[0] + yMC[0]")\
    .Define("zImpact", "tLineParam*pzMC[0] + zMC[0]")\
    .Define("tSmeared", "smear_time(cal.t)")\
    .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact, pxMC[0], pyMC[0], pzMC[0])")\
    .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact)")\
    .Define("tofFit", "tof_fit(cal.t, dToRef, dToLine, cal.layer)")\
    .Define("tofFitCut", "tof_fit(cal.t, dToRef, dToLine, cal.layer, 10, 4.)")\
    .Define("tofFitCyl", "tof_fit_cyl(cal.t, dToRef, dToLine, cal.layer)")\
    .Define("tofFastest", "tof_fastest(cal.t, dToRef)")\
    .Define("tofClosest", "tof_closest(cal.t, dToRef)")\
    .Define("tofAvg", "tof_avg(cal.t, dToRef, dToLine, cal.layer)")\
    .Define("tofFitSmeared", "tof_fit(tSmeared, dToRef, dToLine, cal.layer)")\
    .Define("tofFastestSmeared", "tof_fastest(tSmeared, dToRef)")\
    .Define("tofClosestSmeared", "tof_closest(tSmeared, dToRef)")\
    .Define("tofAvgSmeared", "tof_avg(tSmeared, dToRef, dToLine, cal.layer)")\
    .Define("tofLen", "sqrt((xImpact-xMC[0])*(xImpact-xMC[0]) + (yImpact-yMC[0])*(yImpact-yMC[0]) + (zImpact-zMC[0])*(zImpact-zMC[0]))/c")\
    .Define("len", "sqrt((xImpact-xMC[0])*(xImpact-xMC[0]) + (yImpact-yMC[0])*(yImpact-yMC[0]) + (zImpact-zMC[0])*(zImpact-zMC[0]))")\
    .Define("mom", "sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0] + pzMC[0]*pzMC[0])")\
    .Define("dTfit", "tofFit - tofLen")\
    .Define("dTfitCut", "tofFitCut - tofLen")\
    .Define("dTFitCyl", "tofFitCyl - tofLen")\
    .Define("dTfastest", "tofFastest - tofLen")\
    .Define("dTclosest", "tofClosest - tofLen")\
    .Define("dTavg", "tofAvg - tofLen")\
    .Define("dTfitSmeared", "tofFitSmeared - tofLen")\
    .Define("dTfastestSmeared", "tofFastestSmeared - tofLen")\
    .Define("dTclosestSmeared", "tofClosestSmeared - tofLen")\
    .Define("dTavgSmeared", "tofAvgSmeared - tofLen")\
    .Define("dToRef_fastest", "dToRef_fastest(cal.t, dToRef)")\
    .Define("dToLine_fastest", "dToLine_fastest(cal.t, dToLine)")\
    .Define("layer_fastest", "layer_fastest(cal.t, cal.layer)")\
    .Define("tHitFit", "time_hit_fit(cal.t, dToLine, cal.layer)")\
    .Define("dToLineHitFit", "dToLine_hit_fit(dToLine, cal.layer)")\
    .Define("dToRefHitFit", "dToRef_hit_fit(dToRef, dToLine, cal.layer)")\
    .Define("layerHitFit", "layer_hit_fit(dToLine, cal.layer)")

    return res


def use_cluster_position(df):
    # Calculations of track length and direction and intersection based on cluster position
    res = df.Define("phiLine", "atan2(cluster.y, cluster.x)").Define("thetaLine", "atan2( sqrt(cluster.x*cluster.x + cluster.y*cluster.y), cluster.z)")\
    .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")\
    .Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, cluster.x, cluster.y, cluster.z, cluster.x, cluster.y, cluster.z)")\
    .Define("xImpact", "tLineParam*cluster.x + cluster.x")\
    .Define("yImpact", "tLineParam*cluster.y + cluster.y")\
    .Define("zImpact", "tLineParam*cluster.z + cluster.z")\
    .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact, phiLine, thetaLine)")\
    .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact)")\
    .Define("tofFit", "tof_fit(cal.t, dToRef, dToLine, cal.layer, 10)")\
    .Define("tofLen", "sqrt((xImpact)*(xImpact) + (yImpact)*(yImpact) + (zImpact)*(zImpact))/c")\
    .Define("mom", "sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0] + pzMC[0]*pzMC[0])")\
    .Define("dT", "tofFit - tofLen")
    return res


def use_cluster_direction(df):
    # Calculations of track length and direction and intersection based on cluster direction
    res = df.Define("phiLine", "cluster.phi").Define("thetaLine", "cluster.theta")\
    .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")\
    .Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, sin(thetaLine)*cos(phiLine), sin(thetaLine)*sin(phiLine), cos(thetaLine), cluster.x, cluster.y, cluster.z)")\
    .Define("xImpact", "tLineParam*sin(thetaLine)*cos(phiLine) + cluster.x")\
    .Define("yImpact", "tLineParam*sin(thetaLine)*sin(phiLine) + cluster.y")\
    .Define("zImpact", "tLineParam*cos(thetaLine) + cluster.z")\
    .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact, phiLine, thetaLine)")\
    .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact)")\
    .Define("tofFit", "tof_fit(cal.t, dToRef, dToLine, cal.layer, 10)")\
    .Define("tofLen", "sqrt((xImpact)*(xImpact) + (yImpact)*(yImpact) + (zImpact)*(zImpact))/c")\
    .Define("mom", "sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0] + pzMC[0]*pzMC[0])")\
    .Define("dT", "tofFit - tofLen")
    return res


def analysis():
    canvas = ROOT.TCanvas()
    df = good_photons(df_initial).Range(300000)
    df = use_mc(df)
    # Define all posible histograms you may want to study

    # h_tof_fastest = df.Histo1D(("h_tof_fastest","TOF fastest; TOF, [ns]; N_{PFO}", 600, 5., 11.), "tofFastest")
    # h_tof_closest = df.Histo1D(("h_tof_closest","TOF closest; TOF, [ns]; N_{PFO}", 600, 5., 11.), "tofClosest")

    # h_tof_fastest_bias = df.Histo1D(("h_tof_fastest_bias","TOF fastest; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTfastest")
    # h_tof_closest_bias = df.Histo1D(("h_tof_closest_bias","TOF closest; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTclosest")
    h_tof_fit_bias = df.Histo1D(("h_tof_fit_bias","TOF fit (hit per layer in 10 layers); TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTfit")
    # h_tof_avg_bias = df.Histo1D(("h_tof_avg_bias","TOF avg; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTavg")
    # h_tof_fastest_bias = df.Histo1D(("h_tof_fastest_bias","TOF fastest; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.5, .5), "dTfastestSmeared")
    # h_tof_closest_bias = df.Histo1D(("h_tof_closest_bias","TOF closest; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.5, .5), "dTclosestSmeared")
    # h_tof_fit_bias = df.Histo1D(("h_tof_fit_bias","TOF fit; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.5, .5), "dTfitSmeared")
    # h_tof_avg_bias = df.Histo1D(("h_tof_avg_bias","TOF avg; TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.5, .5), "dTavgSmeared")

    h_tof_fit_cyl_bias = df.Histo1D(("h_tof_fit_cyl_bias","TOF fit (all hits in d<4 mm); TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTFitCyl")
    h_tof_fit_cut_bias = df.Histo1D(("h_tof_fit_cut_bias","TOF fit (hit per layer in 10 layers + d<4 mm cut); TOF - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTfitCut")


    # prof_tof_fit_vs_p = df.Histo2D(("prof_tof_fit_vs_p","Profile TOF bias vs p; p, [GeV]; TOF_{fit} - #frac{l}{c}, [ns]", 1000, 0., 100., 1000, -.05, .05), "mom", "dTfit").ProfileX()
    # prof_tof_fit_vs_nHits = df.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Histo2D(("prof_tof_fit_vs_nHits","Profile TOF bias vs nCaloHits; nCaloHits; TOF_{fit} - #frac{l}{c}, [ns]", 150, 0., 150, 1000, -.05, .05), "nCaloHits", "dTfit")
    # prof_tof_fit_vs_theta = df.Histo2D(("prof_tof_fit_vs_theta","Profile TOF bias vs #theta; #theta; TOF_{fit} - #frac{l}{c}, [ns]", 1000, -6.28, 6.28, 1000, -.05, .05), "thetaLine", "dTfit")
    # prof_p_vs_nHits = df.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Histo2D(("prof_p_vs_nHits","Profile p vs nCaloHits; nCaloHits; p, [GeV]", 150, 0., 150, 1000, 0., 100.), "nCaloHits", "mom")
    # h_tof_bias_8_12 = df.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Filter("8 <= nCaloHits && nCaloHits <= 12").Histo1D(("h_tof_bias_8_12","N_{calo, hits} [8-12]; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "dTfit")
    # h_tof_bias_28_32 = df.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Filter("28 <= nCaloHits && nCaloHits <= 32").Histo1D(("h_tof_bias_28_32","N_{calo, hits} [28-32]; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "dTfit")
    # h_tof_bias_68_72 = df.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Filter("68 <= nCaloHits && nCaloHits <= 72").Histo1D(("h_tof_bias_68_72","N_{calo, hits} [68-72]; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "dTfit")
    # h_tof_len = df.Histo1D(("h_tof_len","TOF len; #frac{l}{c}, [ns]; N_{PFO}", 600, 5., 11.), "tofLen")
    # prof_tof_fastest_vs_phi = df.Histo2D(("prof_tof_fastest_vs_phi","Profile TOF fastest bias vs #phi; #phi, [rad]; TOF_{fastest} - #frac{l}{c}, [ns]", 100, -6.28, 6.28, 1000, -.5, .5), "phiLine", "dTfastest").ProfileX()
    # prof_tof_fastest_vs_theta = df.Histo2D(("prof_tof_fastest_vs_theta","Profile TOF fastest bias vs #theta; #theta, [rad]; TOF_{fastest} - #frac{l}{c}, [ns]", 300, 0, 3.14, 1000, -.5, .5), "thetaLine", "dTfastest").ProfileX()
    # prof_tof_fastest_vs_dToRef = df.Histo2D(("prof_tof_fastest_vs_dToRef","Profile TOF fastest bias vs dToRef; dToRef, [mm]; TOF_{fastest} - #frac{l}{c}, [ns]", 1000, 0, 300., 1000, -2., 2.), "dToRef_fastest", "dTfastest").ProfileX()
    # h_dToRef_fastest = df.Histo1D(("h_dToRef_fastest","dToRef fastest; dToRef, [mm]; N_{PFO}", 1000, 0, 300.), "dToRef_fastest")
    # h_layer_fastest = df.Histo1D(("h_layer_fastest","layer fastest; layer; N_{PFO}", 30, 0, 30), "layer_fastest")
    # h_dToRef_layers_fastest = df.Histo2D(("h_dToRef_layers_fastest","layer fastest; layer; dToRef, [mm]", 30, 0, 30, 400, 0, 200.), "layer_fastest", "dToRef_fastest")
    # h_layer_nCaloHits_fastest = df.Histo2D(("h_dToRef_layers_fastest","layer fastest; layer; nCaloHits (all)", 30, 0, 30, 300, 0, 300), "layer_fastest", "cal.nHits")
    # prof_layer_nCaloHits_fastest = df.Histo2D(("prof_layer_nCaloHits_fastest","layer fastest; layer fastest; nCaloHits (all)", 30, 0, 30, 300, 0, 300), "layer_fastest", "cal.nHits").ProfileX()
    # prof_layer_p_fastest = df.Histo2D(("prof_layer_p_fastest","layer fastest; p (GeV); layer fastest", 1000, 0, 100, 30, 0, 30), "mom", "layer_fastest").ProfileX()
    # h_tof_fit_bias_layer_fastest = df.Histo2D(("h_tof_fit_bias_layer_fastest","TOF fit bias vs layer of fastest TOF; layer; TOF fit bias, [nm]", 30, 0, 30, 1000, -0.5, 0.5), "layer_fastest", "dTfit")
    # h_tof_fastest_bias_smeared = df.Histo1D(("h_tof_fastest_bias_smeared","TOF fastest smeared; TOF_{fastest, smeared} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.4, .4), "dTfastestSmeared")
    # h_tof_closest_bias_smeared = df.Histo1D(("h_tof_closest_bias_smeared","TOF closest smeared; TOF_{closest, smeared} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.4, .4), "dTclosestSmeared")
    # h_tof_fit_bias_smeared = df.Histo1D(("h_tof_fit_bias_smeared","TOF fit smeared; TOF_{fit, smeared} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.4, .4), "dTfitSmeared")
    # # print_failed_fits = df.Filter("dTfit < -0.01").Define("print_status", "fit_plot(cal.t, dToRef, dToLine, cal.layer)").Histo1D(("name", "title", 2, 0,1), "print_status")
    #
    # h_dToLine_fit = df.Histo1D(("h_dToLine_fit","All PFO; distance to track line, [mm]; N_{hits}", 1000, 0., 20.), "dToLineHitFit")
    # h_dToLine_fit_wrong = df.Filter("dTfit < -0.01").Histo1D(("h_dToLine_fit_wrong","Biased PFO; distance to track line, [mm]; N_{hits}", 1000, 0., 20.), "dToLineHitFit")
    #
    # h_t_fit = df.Histo1D(("h_t_fit","All PFO; TOF, [ns]; N_{hits}", 1000, 5., 11.), "tHitFit")
    # h_t_fit_wrong = df.Filter("dTfit < -0.01").Histo1D(("h_t_fit_wrong","Biased PFO; TOF, [ns]; N_{hits}", 1000, 5., 11.), "tHitFit")
    #
    # h_dToRef_fit = df.Filter("dTfit >= -0.01").Histo1D(("h_dToRef_fit","Good PFO; distance to Impact point, [mm]; N_{hits}", 1000, 0., 200.), "dToRefHitFit")
    # h_dToRef_fit_wrong = df.Filter("dTfit < -0.01").Histo1D(("h_dToRef_fit_wrong","Biased PFO; distance to Impact point, [mm]; N_{hits}", 1000, 0., 200.), "dToRefHitFit")
    #
    # h_t_vs_dToRef_fit = df.Filter("dTfit >= -0.01").Histo2D(("h_t_vs_dToRef_fit","Good PFO; distance to Impact point, [mm]; TOF, [ns]", 1000, 0., 200., 1000, 5., 11.), "dToRefHitFit", "tHitFit")
    # h_t_vs_dToRef_fit_wrong = df.Filter("dTfit < -0.01").Histo2D(("h_t_vs_dToRef_fit_wrong","Biased PFO; distance to Impact point, [mm]; TOF, [ns]", 1000, 0., 200., 1000, 5., 11.), "dToRefHitFit", "tHitFit")
    #
    # h_tof_fit_cyl_bias = df.Histo1D(("h_tof_fit_cyl_bias","TOF fit d<4 mm; TOF_{fit, cylinder} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTFitCyl")
    # h_tof_fit_cut_bias = df.Histo1D(("h_tof_fit_cut_bias","TOF fit (10layers) + d<4 mm ; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.1, .1), "dTfitCut")

    # 1, 4, 6, 8
    h_tof_fit_bias.Draw()
    h_tof_fit_bias.SetLineColor(6)

    h_tof_fit_cyl_bias.Draw("same")
    h_tof_fit_cyl_bias.SetLineColor(6)
    h_tof_fit_cyl_bias.SetLineStyle(2)
    h_tof_fit_cut_bias.Draw("same")
    h_tof_fit_cut_bias.SetLineColor(6)
    h_tof_fit_cut_bias.SetLineStyle(8)


    canvas.BuildLegend()
    canvas.SetGridx()
    canvas.SetGridy()
    h_tof_fit_bias.SetTitle("Different fit methods with cylinder cut")
    canvas.Update()
    input("wait")




analysis()
