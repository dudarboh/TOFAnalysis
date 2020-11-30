import ROOT
import time
import numpy as np

# calib_kaons data
# d = ROOT.RDataFrame("ana_tree", "/nfs/dust/ilc/user/dudarboh/final_files/calib_kaons.root")

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

df = ROOT.RDataFrame(ch)
ROOT.gInterpreter.Declare('#include "analysis.hpp"')

def good_pfos(df):
    df_res = df.Filter("nMC == 1\
                  && cal.nHits > 0\
                  && isBackscatter[0] == 0\
                  && abs(cluster.z) < 2200.\
                  && charge != 0.")
    return df_res

def good_kaons(df):
    df_res = df.Filter("nMC == 1\
                  && abs(PDG[0]) == 321\
                  && isBackscatter[0] == 0\
                  && abs(cluster.z) < 2200.\
                  && abs(xMC[0]) < .5\
                  && abs(yMC[0]) < .5\
                  && abs(zMC[0]) < .5")
    return df_res

def good_photons(df):
    df_res = df.Filter("nMC == 1\
                  && PDG[0] == 22\
                  && isBackscatter[0] == 0\
                  && isDecayedInTracker[0] == 0\
                  && abs(cluster.z) < 2200.\
                  && abs(xMC[0]) < .5\
                  && abs(yMC[0]) < .5\
                  && abs(zMC[0]) < .5")
    return df_res


def photons(df):
    # Cuts on nice photons
    start_time = time.time()
    canvas = ROOT.TCanvas()
    hs = ROOT.THStack()
    hs.SetTitle("TOF bias for photons;TOF_{fit} - TOF_{len}, [ns];N_{PFO}")
    df1 = good_photons(df).Range(1000000)

    # Calculations based on true info
    df2 = df1.Define("phiLine", "atan2(pyMC[0], pxMC[0])").Define("thetaLine", "atan2( sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0]), pzMC[0])")\
    .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")\
    .Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, pxMC[0], pyMC[0], pzMC[0], xMC[0], yMC[0], zMC[0])")\
    .Define("xImpact", "tLineParam*pxMC[0] + xMC[0]")\
    .Define("yImpact", "tLineParam*pyMC[0] + yMC[0]")\
    .Define("zImpact", "tLineParam*pzMC[0] + zMC[0]")\
    .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact, pxMC[0], pyMC[0], pzMC[0])")\
    .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact)")\
    .Define("tofFit", "tof_fit(cal.t, dToRef, dToLine, cal.layer, 10)")\
    .Define("tofLen", "sqrt((xImpact-xMC[0])*(xImpact-xMC[0]) + (yImpact-yMC[0])*(yImpact-yMC[0]) + (zImpact-zMC[0])*(zImpact-zMC[0]))/c")\
    .Define("len", "sqrt((xImpact-xMC[0])*(xImpact-xMC[0]) + (yImpact-yMC[0])*(yImpact-yMC[0]) + (zImpact-zMC[0])*(zImpact-zMC[0]))")

    # h2 = df2.Define("diff", "tofFit - tofLen").Histo1D(("h2","TOF bias for photons; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "diff")
    # h2.Draw()
    # hs.Add(h2.GetPtr())

    # dT vs p profile
    # h2 = df2.Define("diff", "tofFit - tofLen").Define("mom", "sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0] + pzMC[0]*pzMC[0])")\
    #         .Histo2D(("h2","p vs TOF bias Profile;p, [GeV];TOF_{fit} - TOF_{len}, [ns]", 1000, 0., 100., 1000, -.05, .05), "mom", "diff")

    # dT vs N calo hits profile
    # h2 = df2.Define("diff", "tofFit - tofLen").Define("nCaloHits", "nCaloHits(cal.layer, 10)")\
    #         .Histo2D(("h2","nCaloHits vs TOF bias Profile;nCaloHits;TOF_{fit} - TOF_{len}, [ns]", 150, 0., 150, 1000, -.05, .05), "nCaloHits", "diff")

    # theta profile
    # h2 = df2.Define("diff", "tofFit - tofLen")\
    #         .Histo2D(("h2","#theta vs TOF bias Profile;#theta ;TOF_{fit} - TOF_{len}, [ns]", 1000, -6.28, 6.28, 1000, -.05, .05), "thetaLine", "diff")

    # p vs N calo hits profile
    # h2 = df2.Define("mom", "sqrt(pxMC[0]*pxMC[0] + pyMC[0]*pyMC[0] + pzMC[0]*pzMC[0])").Define("nCaloHits", "nCaloHits(cal.layer, 10)")\
    #         .Histo2D(("h2","p vs nCaloHits Profile; nCaloHits; p, [GeV]", 150, 0., 150, 1000, 0., 100.), "nCaloHits", "mom")
    # prof = h2.ProfileX()
    # prof.Draw()

    # 4 histo from a profile with different dT
    h2 = df2.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Filter("8 <= nCaloHits && nCaloHits <= 12").Define("diff", "tofFit - tofLen").Histo1D(("h2"," N_{calo, hits} [8-12] ; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "diff")
    h3 = df2.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Filter("28 <= nCaloHits && nCaloHits <= 32").Define("diff", "tofFit - tofLen").Histo1D(("h3"," N_{calo, hits} [28-32] ; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "diff")
    h4 = df2.Define("nCaloHits", "nCaloHits(cal.layer, 10)").Filter("68 <= nCaloHits && nCaloHits <= 72").Define("diff", "tofFit - tofLen").Histo1D(("h4"," N_{calo, hits} [68-72] ; TOF_{fit} - #frac{l}{c}, [ns]; N_{PFO}", 1000, -.05, .05), "diff")

    h2.Scale(1./h2.Integral())
    h2.Draw()

    h3.Scale(1./h3.Integral())
    h3.Draw("histo same")
    h3.SetLineColor(2)

    h4.Scale(1./h4.Integral())
    h4.Draw("histo same")
    h4.SetLineColor(3)

    # Calculations based on cluster position
    # df3 = df1.Define("phiLine", "atan2(cluster.y, cluster.x)").Define("thetaLine", "atan2( sqrt(cluster.x*cluster.x + cluster.y*cluster.y), cluster.z)")\
    # .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")\
    # .Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, cluster.x, cluster.y, cluster.z, cluster.x, cluster.y, cluster.z)")\
    # .Define("xImpact", "tLineParam*cluster.x + cluster.x")\
    # .Define("yImpact", "tLineParam*cluster.y + cluster.y")\
    # .Define("zImpact", "tLineParam*cluster.z + cluster.z")\
    # .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact, phiLine, thetaLine)")\
    # .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact)")\
    # .Define("tofFit", "tof_fit(cal.t, dToRef, dToLine, cal.layer, 10)")\
    # .Define("tofLen", "sqrt((xImpact)*(xImpact) + (yImpact)*(yImpact) + (zImpact)*(zImpact))/c")
    #
    # h3 = df3.Define("diff", "tofFit - tofLen").Histo1D(("h3","Cluster pos", 1000, -.05, .05), "diff")
    # h3.SetLineColor(2)
    # hs.Add(h3.GetPtr())
    #
    # # Calculations based on cluster direction
    # df4 = df1.Define("phiLine", "cluster.phi").Define("thetaLine", "cluster.theta")\
    # .Define("nx", "getPlaneDir(phiLine, 0)").Define("ny", "getPlaneDir(phiLine, 1)").Define("nz", "getPlaneDir(phiLine, 2)")\
    # .Define("tLineParam", "intersection(nx, ny, nz, nx, ny, nz, sin(thetaLine)*cos(phiLine), sin(thetaLine)*sin(phiLine), cos(thetaLine), cluster.x, cluster.y, cluster.z)")\
    # .Define("xImpact", "tLineParam*sin(thetaLine)*cos(phiLine) + cluster.x")\
    # .Define("yImpact", "tLineParam*sin(thetaLine)*sin(phiLine) + cluster.y")\
    # .Define("zImpact", "tLineParam*cos(thetaLine) + cluster.z")\
    # .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact, phiLine, thetaLine)")\
    # .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, xImpact, yImpact, zImpact)")\
    # .Define("tofFit", "tof_fit(cal.t, dToRef, dToLine, cal.layer, 10)")\
    # .Define("tofLen", "sqrt((xImpact)*(xImpact) + (yImpact)*(yImpact) + (zImpact)*(zImpact))/c")
    #
    # h4 = df4.Define("diff", "tofFit - tofLen").Histo1D(("h4","Cluster direction", 1000, -.05, .05), "diff")
    # hs.Add(h4.GetPtr())
    # h4.SetLineColor(3)




    # hs.Draw("nostack")
    canvas.BuildLegend()
    canvas.Update()
    input("wait")


def kaons(df):
    # Cuts on nice kaons
    start_time = time.time()
    df1 = good_kaons(df).Range(10000)

    # histo = df1.Define("length", "abs((piFit.phi-piFit.phiCalState)/piFit.omegaCalState)*sqrt(1. + piFit.tanLCalState*piFit.tanLCalState)")\
    # .Define("momentum", "sqrt(piFit.px*piFit.px + piFit.py*piFit.py + piFit.pz*piFit.pz)")\
    # .Histo2D(("histo", "title", 500, 0., 15., 500, 1700., 5000.), "momentum", "length")

    # histo = df1.Define("momentumTrack", "sqrt(piFit.px*piFit.px + piFit.py*piFit.py + piFit.pz*piFit.pz)")\
    # .Define("momentumPFO", "sqrt(px*px + py*py + pz*pz)").Define("diff", "momentumPFO - momentumTrack")\
    # .Histo1D(("histo", "title", 1000, -15., 15.), "diff")

    histo = df1.Histo1D(("histo", "title", 5000, -500., 500.), "cal.x")


    histo.Draw()
    # print("Execution time: {} for {} entries".format(time.time() - start_time, df1.Count().GetValue()))
    input("wait")


def pfos(df):
    # Cuts on nice photons
    start_time = time.time()
    canvas = ROOT.TCanvas()
    df1 = good_pfos(df).Range(1000)

    h_pip_lip_avg = df1.Define("mom", "sqrt(piFit.px*piFit.px + piFit.py*piFit.py + piFit.pz*piFit.pz)")\
                    .Define("length", "abs(( piFit.phi - piFit.phiCalState ) / piFit.omega) * sqrt(1. + piFit.tanL*piFit.tanL )")\
                    .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState)")\
                    .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState, piFit.pxCalState, piFit.pyCalState, piFit.pzCalState)")\
                    .Define("tof", "tof_avg(cal.t, dToRef, dToLine, cal.layer)")\
                    .Define("beta", "length/(tof * c)")\
                    .Define("mass", "mom / beta * sqrt(1. - beta * beta) * 1000")\
                    .Histo1D(("h_pip_lip_avg","title; m, [MeV]; N_{PFOs}", 2000, 0., 1000.), "mass")

    h_pcal_lcal_avg = df1.Define("mom", "sqrt(piFit.pxCalState*piFit.pxCalState + piFit.pyCalState*piFit.pyCalState + piFit.pzCalState*piFit.pzCalState)")\
                         .Define("length", "abs((piFit.phi-piFit.phiCalState)/piFit.omegaCalState)*sqrt(1. + piFit.tanLCalState*piFit.tanLCalState)")\
                         .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState)")\
                         .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState, piFit.pxCalState, piFit.pyCalState, piFit.pzCalState)")\
                         .Define("tof", "tof_avg(cal.t, dToRef, dToLine, cal.layer)")\
                         .Define("beta", "length/(tof * c)")\
                         .Define("mass", "mom / beta * sqrt(1. - beta * beta) * 1000")\
                         .Histo1D(("h_pcal_lcal_avg","title; m, [MeV]; N_{PFOs}", 2000, 0., 1000.), "mass")
    h_pcal_lcal_avg.SetLineColor(2)

    h_pcal_lcal_fit = df1.Define("mom", "sqrt(piFit.pxCalState*piFit.pxCalState + piFit.pyCalState*piFit.pyCalState + piFit.pzCalState*piFit.pzCalState)")\
                         .Define("length", "abs((piFit.phi-piFit.phiCalState)/piFit.omegaCalState)*sqrt(1. + piFit.tanLCalState*piFit.tanLCalState)")\
                         .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState)")\
                         .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState, piFit.pxCalState, piFit.pyCalState, piFit.pzCalState)")\
                         .Define("tof", "tof_fit(cal.t, dToRef, dToLine, cal.layer)")\
                         .Define("beta", "length/(tof * c)")\
                         .Define("mass", "mom / beta * sqrt(1. - beta * beta) * 1000")\
                         .Histo1D(("h_pcal_lcal_fit","title; m, [MeV]; N_{PFOs}", 2000, 0., 1000.), "mass")
    h_pcal_lcal_fit.SetLineColor(3)

    h_pip_lip_fit = df1.Define("mom", "sqrt(piFit.px*piFit.px + piFit.py*piFit.py + piFit.pz*piFit.pz)")\
                         .Define("length", "abs((piFit.phi-piFit.phiCalState)/piFit.omega)*sqrt(1. + piFit.tanL*piFit.tanL)")\
                         .Define("dToRef", "dToRef(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState)")\
                         .Define("dToLine", "dToLine(cal.x, cal.y, cal.z, piFit.xRefCalState, piFit.yRefCalState, piFit.zRefCalState, piFit.pxCalState, piFit.pyCalState, piFit.pzCalState)")\
                         .Define("tof", "tof_fit(cal.t, dToRef, dToLine, cal.layer)")\
                         .Define("beta", "length/(tof * c)")\
                         .Define("mass", "mom / beta * sqrt(1. - beta * beta) * 1000")\
                         .Histo1D(("h_pip_lip_fit","title; m, [MeV]; N_{PFOs}", 2000, 0., 1000.), "mass")
    h_pip_lip_fit.SetLineColor(4)


    h_pip_lip_avg.Draw()
    h_pcal_lcal_avg.Draw("same")
    h_pcal_lcal_fit.Draw("same")
    h_pip_lip_fit.Draw("same")
    # h_pip_lip_avg.SetMaximum(10000)







    x = np.array([139.570183535, 493.677, 938.2720881629])
    line_pion = ROOT.TLine(x[0], 0., x[0], 1000.)
    line_pion.SetLineStyle(9)
    line_pion.Draw()

    line_kaon = ROOT.TLine(x[1], 0., x[1], 1000.)
    line_kaon.SetLineStyle(9)
    line_kaon.Draw()

    line_proton = ROOT.TLine(x[2], 0., x[2], 1000.)
    line_proton.SetLineStyle(9)
    line_proton.Draw()

    latex = ROOT.TLatex()
    latex.SetTextSize(0.08)
    latex.DrawLatex(x[0] + 10, 1000,"#pi^{#pm}");
    latex.DrawLatex(x[1] + 10, 1000,"K^{#pm}");
    latex.DrawLatex(x[2] + 10, 1000,"p^{#pm}");
    canvas.Update()
    input("wait")


def bias_plots():
    canvas = ROOT.TCanvas()
    x = np.array([139.570183535, 493.677, 938.2720881629])
    # p_ip, len_ip, TOF avg
    # y_pip_lip_avg = np.array([1.18410e+02, 4.90807e+02, 9.44576e+02]) - x
    y_pip_lip_avg = np.array([118., 492., 945.8]) - x
    gr_pip_lip_avg = ROOT.TGraph(3, x, y_pip_lip_avg)
    gr_pip_lip_avg.SetTitle("p_{IP}, #Omega_{IP}, #lambda_{IP}, TOF_{avg}")
    gr_pip_lip_avg.SetMarkerStyle(20)
    gr_pip_lip_avg.SetMinimum(-25)
    gr_pip_lip_avg.SetMaximum(16)
    gr_pip_lip_avg.Draw("APL")

    # y_pcal_lcal_avg = np.array([112.2, 503.46, 952.15]) - x
    # gr_pcal_lcal_avg = ROOT.TGraph(3, x, y_pcal_lcal_avg)
    # gr_pcal_lcal_avg.SetTitle("p_{cal}, #Omega_{cal}, #lambda_{cal}, TOF_{avg}")
    # gr_pcal_lcal_avg.SetMarkerColor(3)
    # gr_pcal_lcal_avg.SetLineColor(3)
    # gr_pcal_lcal_avg.SetMarkerStyle(21)
    # gr_pcal_lcal_avg.SetMarkerSize(1.5)
    # gr_pcal_lcal_avg.Draw("PLsame")
    #
    # y_pip_lip_fit = np.array([151.5, 484., 932.7]) - x
    # gr_pip_lip_fit = ROOT.TGraph(3, x, y_pip_lip_fit)
    # gr_pip_lip_fit.SetTitle("p_{IP}, #Omega_{IP}, #lambda_{IP}, TOF_{fit}")
    # gr_pip_lip_fit.SetMarkerColor(4)
    # gr_pip_lip_fit.SetLineColor(4)
    # gr_pip_lip_fit.SetMarkerStyle(22)
    # gr_pip_lip_fit.SetMarkerSize(1.5)
    # gr_pip_lip_fit.Draw("PLsame")

    # y_pcal_lcal_fit = np.array([144.7, 494.697, 939.06]) - x
    # gr_pcal_lcal_fit = ROOT.TGraph(3, x, y_pcal_lcal_fit)
    # gr_pcal_lcal_fit.SetTitle("p_{cal}, #Omega_{cal}, #lambda_{cal}, TOF_{fit}")
    # gr_pcal_lcal_fit.SetMarkerColor(8)
    # gr_pcal_lcal_fit.SetLineColor(8)
    # gr_pcal_lcal_fit.SetLineWidth(3)
    # gr_pcal_lcal_fit.SetMarkerStyle(23)
    # gr_pcal_lcal_fit.SetMarkerSize(1.5)
    # gr_pcal_lcal_fit.Draw("PLsame")


    canvas.BuildLegend()
    gr_pip_lip_avg.SetTitle("Mass bias for #pi^{#pm}, K^{#pm}, p; m_{true}, [MeV]; m_{reco} - m_{true}, [MeV]")








    line_pion = ROOT.TLine(x[0], gr_pip_lip_avg.GetMinimum(), x[0], gr_pip_lip_avg.GetMaximum())
    line_pion.SetLineStyle(9)
    line_pion.Draw()

    line_kaon = ROOT.TLine(x[1], gr_pip_lip_avg.GetMinimum(), x[1], gr_pip_lip_avg.GetMaximum())
    line_kaon.SetLineStyle(9)
    line_kaon.Draw()

    line_proton = ROOT.TLine(x[2], gr_pip_lip_avg.GetMinimum(), x[2], gr_pip_lip_avg.GetMaximum())
    line_proton.SetLineStyle(9)
    line_proton.Draw()

    latex = ROOT.TLatex()
    latex.SetTextSize(0.08)
    latex.DrawLatex(x[0] + 10, gr_pip_lip_avg.GetMaximum() - 10,"#pi^{#pm}");
    latex.DrawLatex(x[1] + 10, gr_pip_lip_avg.GetMaximum() - 10,"K^{#pm}");
    latex.DrawLatex(x[2] + 10, gr_pip_lip_avg.GetMaximum() - 10,"p^{#pm}");

    canvas.SetGridy()
    canvas.Update()
    input("wait")









photons(df)
# kaons(df)
# pfos(df)

# bias_plots()
