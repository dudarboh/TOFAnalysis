import ROOT
import numpy as np
import time
import math
ROOT.EnableImplicitMT(2)

ROOT.gInterpreter.Declare('''

#define SPEED_OF_LIGHT 299.792458

double getOldMass(double phi1, double phi2, double omega, double tanL, double tof, double mom){

    double dPhi = std::abs(phi2 - phi1);
    // if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    double trackLength = dPhi / std::abs(omega) * std::sqrt(1. + tanL*tanL);
    double beta = trackLength / (tof * SPEED_OF_LIGHT);

    if (1 - beta*beta < 0) return -999.;
    double mass = mom / beta * std::sqrt(1 - beta*beta);
    return mass;
}

double getNewMass1(double mom, double beta){
    if (beta > 1.)return -999.;
    return mom / beta * std::sqrt(1. - beta*beta);
}

double getNewMass2(double mom, double beta){
    if (beta > 1.)return -999.;
    return std::sqrt(2.*mom*mom*(1./beta - 1.));
}


''')

ROOT.gStyle.SetNdivisions(512)
ROOT.gStyle.SetNdivisions(512, "Y")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetMarkerStyle(20)
ROOT.gStyle.SetMarkerSize(1.2)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)
ROOT.gStyle.SetPalette(57)

canvas_glob = ROOT.TCanvas()
test_canvas = ROOT.TCanvas()

# Get true beta vs p curves

def get_curves(df, mom_str, method, n_mom_bins=30, to_draw=True):
    graphs = {}
    gr_diff = {}
    gr_sep = {}
    pdg_line = {}
    pdgs = [211, 321, 2212]
    colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kGreen+2]
    m_pdg = {211 : 0.13957039, 321 : 0.493677, 2212 : 0.938272088}

    if to_draw:
        canvas = ROOT.TCanvas()
        canvas.Divide(2, 2)
        for i in range(1, 5):
            canvas.cd(i)
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()


    # 2D histo
    h_all = df.Histo2D(("h_all", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 0., 15., 200, -0.1, 1.3), mom_str, method)
    if to_draw:
        canvas.cd(1)
        ROOT.gPad.SetLogz()
        h_all.SetStats(1)
        h_all.Draw("colz")

    for k, pdg in enumerate(pdgs):
        if (pdg == 211):
            h = df.Filter("abs(pdg) == {}".format(pdg))\
            .Histo2D(("h", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 0., 15., 200, -0.1, 1.3), mom_str, method)
        elif (pdg == 321):
            h = df.Filter("abs(pdg) == {}".format(pdg))\
            .Histo2D(("h", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 0., 15., 200, -0.1, 1.3), mom_str, method)
        elif (pdg == 2212):
            h = df.Filter("abs(pdg) == {}".format(pdg))\
            .Histo2D(("h", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 0., 15., 200, -0.1, 1.3), mom_str, method)


        graphs[pdg] = ROOT.TGraphErrors()
        graphs[pdg].SetMarkerStyle(20)
        graphs[pdg].SetMarkerColor(colors[k])
        graphs[pdg].SetLineColor(colors[k])

        gr_diff[pdg] = ROOT.TGraphErrors()
        gr_diff[pdg].SetMarkerStyle(20)
        gr_diff[pdg].SetMarkerColor(colors[k])
        gr_diff[pdg].SetLineColor(colors[k])
        if pdg == 211:
            graphs[pdg].SetTitle("; p (GeV); mass (GeV)")
            gr_diff[pdg].SetTitle("#pi^{#pm}; p (GeV); mass - m_{PDG} (GeV)")
        elif pdg == 321:
            graphs[pdg].SetTitle("; p (GeV); mass (GeV)")
            gr_diff[pdg].SetTitle("K^{#pm}; p (GeV); mass - m_{PDG} (GeV)")
        elif pdg == 2212:
            graphs[pdg].SetTitle("; p (GeV); mass (GeV)")
            gr_diff[pdg].SetTitle("p; p (GeV); mass - m_{PDG} (GeV)")

        # *************Fit in Slices section***********
        points = 0
        mom_points = []
        for bin in range(1, n_mom_bins + 1):
            hp = h.ProjectionY("hp", bin, bin)
            if hp.GetEntries() < 100:
                continue
            hp.SetStats(0)
            max_bin = hp.GetXaxis().GetBinCenter( hp.GetMaximumBin() )

            func = ROOT.TF1("func", "gaus", max_bin-0.5, max_bin+0.5)
            func.SetParameter(1, max_bin)
            func.SetParLimits(1, max_bin-0.1, max_bin+0.1)
            func.SetParameter(2, 0.03)
            # func.SetParLimits(2, 0., 0.2)
            func.SetNpx(1000)
            hp.Fit(func, "QRN")
            mom =  h.GetXaxis().GetBinCenter( bin )
            mean_fit = func.GetParameter(1)
            sigma_fit = func.GetParameter(2)

            graphs[pdg].SetPoint( points, mom, mean_fit )
            graphs[pdg].SetPointError( points, 0., sigma_fit )
            gr_diff[pdg].SetPoint( points, mom, mean_fit - m_pdg[pdg] )
            # gr_diff[pdg].SetPointError( points, 0., sigma_fit )
            mom_points.append(mom)
            points += 1
        if to_draw:
            graphs[pdg].GetXaxis().SetRangeUser(0., 15.)
            graphs[pdg].GetYaxis().SetRangeUser(-0.1, 1.3)
            gr_diff[pdg].GetXaxis().SetRangeUser(0., 15.)
            gr_diff[pdg].GetYaxis().SetRangeUser(-0.2, 0.2)
            canvas.cd(2)
            if pdg == 211:
                graphs[pdg].Draw("APE")
            pdg_line[pdg] = ROOT.TLine(0., m_pdg[pdg], 15., m_pdg[pdg])
            pdg_line[pdg].SetLineColor(4)
            pdg_line[pdg].Draw()
            graphs[pdg].Draw("PEsame")
            test_canvas.cd()
            if pdg == 211:
                graphs[pdg].Draw("APE")
            pdg_line[pdg] = ROOT.TLine(0., m_pdg[pdg], 15., m_pdg[pdg])
            pdg_line[pdg].SetLineColor(4)
            pdg_line[pdg].Draw()
            graphs[pdg].Draw("PEsame")

            canvas.cd(4)
            gr_diff[pdg].Draw("APE0" if pdg == 211 else "PE0same")
    if to_draw:
        canvas.cd(4)
        ROOT.gPad.BuildLegend(0.25, 0.7, .35, .9)
        gr_diff[211].SetTitle("")



    # *************Get separation power plots***********
    gr_sep["pik"] = ROOT.TGraphErrors()
    gr_sep["pik"].SetMarkerStyle(20)
    gr_sep["pik"].SetMarkerColor(12)
    gr_sep["pik"].SetLineColor(12)
    gr_sep["pik"].SetTitle("#pi/K; p (GeV); Sep. power")
    gr_sep["kp"] = ROOT.TGraphErrors()
    gr_sep["kp"].SetMarkerStyle(20)
    gr_sep["kp"].SetMarkerColor(46)
    gr_sep["kp"].SetLineColor(46)
    gr_sep["kp"].SetTitle("K/p; p (GeV); Sep. power")

    for p in range(points):
        mean_pi = graphs[211].GetPointY(p)
        sigma_pi = graphs[211].GetErrorY(p)
        mean_k = graphs[321].GetPointY(p)
        sigma_k = graphs[321].GetErrorY(p)
        mean_p = graphs[2212].GetPointY(p)
        sigma_p = graphs[2212].GetErrorY(p)

        sep_power_pik = abs(mean_pi - mean_k) / np.sqrt( (sigma_pi*sigma_pi + sigma_k*sigma_k)/2. )
        sep_power_kp = abs(mean_k - mean_p) / np.sqrt( (sigma_k*sigma_k + sigma_p*sigma_p)/2. )

        gr_sep["pik"].SetPoint( p, mom_points[p] , sep_power_pik )
        gr_sep["kp"].SetPoint( p, mom_points[p] , sep_power_kp )

        gr_sep["pik"].GetXaxis().SetRangeUser(0., 15.)
        gr_sep["pik"].GetYaxis().SetRangeUser(0., 15.)
    if to_draw:
        canvas.cd(3)
        gr_sep["pik"].Draw("APL")
        gr_sep["kp"].Draw("PLsame")
        ROOT.gPad.BuildLegend(0.7, 0.7, .9, .9)
        gr_sep["pik"].SetTitle("")
        for i in range(1, 5):
            canvas.cd(i)
            ROOT.gPad.Update()

        # canvas.Print("method_{}_{}_ps.png".format(m, s))
        print("Finished with method", method)
        input("wait")

    return graphs, gr_diff, gr_sep




df = ROOT.RDataFrame("TOFAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/final.root")
df = df.Filter("n_ecal_hits > 0 && abs(ts_ecal_pos.z()) < 2385. && (ts_ecal_pos - pos_closest).r() < 4000.")\
        .Define("mom_ecal", "ts_ecal_mom.r()")\

df = df.Define("mass_ecal", "getOldMass(ts_ip_phi, ts_ecal_phi, ts_ecal_omega, ts_ecal_tanL, tof_closest_0, mom_ecal)")

resolutions = [0, 10, 30, 50, 100]

histos = []
for i,res in enumerate(resolutions):
    df = df.Define("beta_{}".format(res), "track_length_ecal/(tof_closest_{}*SPEED_OF_LIGHT)".format(res))\
           .Define("mass_{}".format(res), "getNewMass1(mom_hm_ecal, beta_{})".format(res))\

#     histos.append( df.Histo1D(("h_{}".format(res), "{} ps; mass (MeV); N pfo".format(res), 1000, 0, 1000), "mass_{}".format(res)) )
#     histos[-1].Draw("" if i == 0 else "same")
#     histos[-1].SetLineColor(ROOT.kRed + 4 - i )
#
# lines = []
# for i in np.array([0.13957039, 0.493677, 0.938272088])*1000.:
#     lines.append(ROOT.TLine(i, 0., i, 15000.))
#     lines[-1].SetLineColor(4)
#     lines[-1].SetLineStyle(9)
#     lines[-1].SetLineWidth(3)
#     lines[-1].Draw()
#
# canvas.Update()
# input("wait")

gr_refit, gr_diff_refit, gr_sep_refit_0 = get_curves(df, "mom_ecal", "mass_ecal", n_mom_bins=30, to_draw=True)
gr_refit, gr_diff_refit, gr_sep_refit_0 = get_curves(df, "mom_hm_ecal", "mass_0", n_mom_bins=30, to_draw=True)
# gr_refit, gr_diff_refit, gr_sep_refit_0 = get_curves(df, "mom_ecal", "mass4_0", n_mom_bins=30, to_draw=True)


# gr_ecal, gr_diff_ecal, gr_sep_ecal = get_curves(df, "mom_ecal", "mass_ecal", n_mom_bins=30, to_draw=True)
# gr_refit, gr_diff_refit, gr_sep_refit_0 = get_curves(df, "mom_hm_ecal", "mass_0", n_mom_bins=30, to_draw=False)
# gr_refit, gr_diff_refit, gr_sep_refit_10 = get_curves(df, "mom_hm_ecal", "mass_10", n_mom_bins=30, to_draw=False)
# gr_refit, gr_diff_refit, gr_sep_refit_30 = get_curves(df, "mom_hm_ecal", "mass_30", n_mom_bins=30, to_draw=False)
# gr_refit, gr_diff_refit, gr_sep_refit_50 = get_curves(df, "mom_hm_ecal", "mass_50", n_mom_bins=30, to_draw=False)
# gr_refit, gr_diff_refit, gr_sep_refit_100 = get_curves(df, "mom_hm_ecal", "mass_100", n_mom_bins=30, to_draw=False)
#
# canvas_glob.cd()
# gr_sep_refit_0["pik"].Draw("APL")
# gr_sep_refit_0["pik"].SetLineColor(ROOT.kRed + 4)
#
# gr_sep_refit_10["pik"].Draw("PLsame")
# gr_sep_refit_10["pik"].SetLineColor(ROOT.kRed + 3)
#
# gr_sep_refit_30["pik"].Draw("PLsame")
# gr_sep_refit_30["pik"].SetLineColor(ROOT.kRed + 2)
#
# gr_sep_refit_50["pik"].Draw("PLsame")
# gr_sep_refit_50["pik"].SetLineColor(ROOT.kRed + 1)
#
# gr_sep_refit_100["pik"].Draw("PLsame")
# gr_sep_refit_100["pik"].SetLineColor(ROOT.kRed)

#
# gr_sep_ecal["kp"].Draw("PLsame")
# gr_sep_ecal["kp"].SetLineColor(4)
# gr_sep_refit["kp"].Draw("PLsame")
# gr_sep_refit["kp"].SetLineColor(6)
#
# canvas_glob.Update()
# input("wait")

def plot_matrix(df, mom_string, mass_string):
    gr, gr_diff, gr_sep = get_curves(df, mom_string, mass_string, n_mom_bins=30, to_draw=False)
    n_points = gr[211].GetN()

    matrix = np.zeros((n_points, 4, 4))

    mom_bins = [ gr[211].GetPointX(i) for i in range(n_points) ]
    mom_bins = np.array(mom_bins)
    data = df.Filter("{0} > 0. && {0} <= 15.".format(mom_string)).AsNumpy(["pdg", mom_string, mass_string])

    for i, (pdg, mom, mass) in enumerate( zip(data["pdg"], data[mom_string], data[mass_string]) ):
        # if i > 200000:
        #     break
        if i%10000 == 0:
            print("Event", i)


        # closeset momentum bin to fill matrix
        closest_mom_bin = int(np.argmin( abs(mom_bins - mom) ))

        # true particle bin
        if abs(pdg) == 211:
            true_pfo = 0
        elif abs(pdg) == 321:
            true_pfo = 1
        elif abs(pdg) == 2212:
            true_pfo = 2
        else:
            true_pfo = 3

        # track length 0, cant reconstruct
        if mass == 0 :
            matrix[closest_mom_bin, true_pfo, 3] += 1
            continue

        if mass < 0:
            matrix[closest_mom_bin, true_pfo, 0] += 1
            continue

        d_to_pi = abs( mass - gr[211].Eval(mom) ) / gr[211].GetErrorY( closest_mom_bin )
        d_to_k = abs( mass - gr[321].Eval(mom) ) / gr[321].GetErrorY( closest_mom_bin )
        d_to_p = abs( mass - gr[2212].Eval(mom) ) / gr[2212].GetErrorY( closest_mom_bin )
        reco_pfo = int( np.argmin(np.array([d_to_pi, d_to_k, d_to_p])) )

        matrix[closest_mom_bin, true_pfo, reco_pfo] += 1


    # Get Full matrix into histo for the cross check
    h = ROOT.TH2F("name", ";;", 4, 0, 4, 4, 0, 4)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h.SetMarkerSize(2)
    h.GetXaxis().SetBinLabel(1, "#pi^{#pm}_{true}")
    h.GetXaxis().SetBinLabel(2, "K^{#pm}_{true}")
    h.GetXaxis().SetBinLabel(3, "p_{true}")
    h.GetXaxis().SetBinLabel(4, "other")

    h.GetYaxis().SetBinLabel(1, "#pi^{#pm}_{reco}")
    h.GetYaxis().SetBinLabel(2, "K^{#pm}_{reco}")
    h.GetYaxis().SetBinLabel(3, "p_{reco}")
    h.GetYaxis().SetBinLabel(4, "undefined")

    mom_summed = np.sum(matrix, axis=0)
    for i in range(4):
        for j in range(4):
            h.SetBinContent(i+1, j+1, mom_summed[i, j])

    for c in range(1, 5):
        norm = h.GetBinContent(c, 1) + h.GetBinContent(c, 2) + h.GetBinContent(c, 3) + h.GetBinContent(c, 4)
        for r in range(1, 5):
            h.SetBinContent(c, r, h.GetBinContent(c, r) / norm)
    # h.Draw("colz text")
    # input("wait")
    eff_pi = matrix[:, 0, 0] / (matrix[:, 0, 0] + matrix[:, 0, 1] + matrix[:, 0, 2] + matrix[:, 0, 3])
    eff_k = matrix[:, 1, 1] / (matrix[:, 1, 0] + matrix[:, 1, 1] + matrix[:, 1, 2] + matrix[:, 1, 3])
    eff_p = matrix[:, 2, 2] / (matrix[:, 2, 0] + matrix[:, 2, 1] + matrix[:, 2, 2] + matrix[:, 2, 3])

    miss_id_pi = (matrix[:, 1, 0] + matrix[:, 2, 0] + matrix[:, 3, 0]) / (matrix[:, 0, 0] + matrix[:, 1, 0] + matrix[:, 2, 0] + matrix[:, 3, 0])
    miss_id_k = (matrix[:, 0, 1] + matrix[:, 2, 1] + matrix[:, 3, 1]) / (matrix[:, 0, 1] + matrix[:, 1, 1] + matrix[:, 2, 1] + matrix[:, 3, 1])
    miss_id_p = (matrix[:, 0, 2] + matrix[:, 1, 2] + matrix[:, 3, 2]) / (matrix[:, 0, 2] + matrix[:, 1, 2] + matrix[:, 2, 2] + matrix[:, 3, 2])

    return mom_bins, eff_pi, eff_k, eff_p, miss_id_pi, miss_id_k, miss_id_p



# mom_bins_old, eff_pi_old, eff_k_old, eff_p_old, miss_id_pi_old, miss_id_k_old, miss_id_p_old = plot_matrix(df, "mom_ecal", "mass_ecal")
mom_bins_new_0, eff_pi_new_0, eff_k_new_0, eff_p_new_0, miss_id_pi_new_0, miss_id_k_new_0, miss_id_p_new_0 = plot_matrix(df, "mom_hm_ecal", "mass_0")
mom_bins_new_10, eff_pi_new_10, eff_k_new_10, eff_p_new_10, miss_id_pi_new_10, miss_id_k_new_10, miss_id_p_new_10 = plot_matrix(df, "mom_hm_ecal", "mass_10")
mom_bins_new_30, eff_pi_new_30, eff_k_new_30, eff_p_new_30, miss_id_pi_new_30, miss_id_k_new_30, miss_id_p_new_30 = plot_matrix(df, "mom_hm_ecal", "mass_30")
mom_bins_new_50, eff_pi_new_50, eff_k_new_50, eff_p_new_50, miss_id_pi_new_50, miss_id_k_new_50, miss_id_p_new_50 = plot_matrix(df, "mom_hm_ecal", "mass_50")
mom_bins_new_100, eff_pi_new_100, eff_k_new_100, eff_p_new_100, miss_id_pi_new_100, miss_id_k_new_100, miss_id_p_new_100 = plot_matrix(df, "mom_hm_ecal", "mass_100")


gr_0 = ROOT.TGraph(len(mom_bins_new_0), mom_bins_new_0, miss_id_pi_new_0)
gr_0.Draw("APL")
gr_0.SetLineColor(ROOT.kRed + 4)
gr_0.SetMarkerColor(ROOT.kRed + 4)
gr_0.SetTitle("0 ps")
gr_10 = ROOT.TGraph(len(mom_bins_new_10), mom_bins_new_10, miss_id_pi_new_10)
gr_10.SetLineColor(ROOT.kRed + 3)
gr_10.SetMarkerColor(ROOT.kRed + 3)
gr_10.SetTitle("10 ps")
gr_10.Draw("PLsame")
gr_30 = ROOT.TGraph(len(mom_bins_new_30), mom_bins_new_30, miss_id_pi_new_30)
gr_30.SetLineColor(ROOT.kRed + 2)
gr_30.SetMarkerColor(ROOT.kRed + 2)
gr_30.SetTitle("30 ps")
gr_30.Draw("PLsame")

gr_50 = ROOT.TGraph(len(mom_bins_new_50), mom_bins_new_50, miss_id_pi_new_50)
gr_50.SetLineColor(ROOT.kRed + 1)
gr_50.SetMarkerColor(ROOT.kRed + 1)
gr_50.SetTitle("50 ps")
gr_50.Draw("PLsame")

gr_100 = ROOT.TGraph(len(mom_bins_new_100), mom_bins_new_100, miss_id_pi_new_100)
gr_100.SetLineColor(ROOT.kRed)
gr_100.SetMarkerColor(ROOT.kRed)
gr_100.SetTitle("100 ps")
gr_100.Draw("PLsame")
input("wait")


gr_0 = ROOT.TGraph(len(mom_bins_new_0), mom_bins_new_0, miss_id_k_new_0)
gr_0.Draw("APL")
gr_0.SetLineColor(ROOT.kRed + 4)
gr_0.SetMarkerColor(ROOT.kRed + 4)
gr_0.SetTitle("0 ps")
gr_10 = ROOT.TGraph(len(mom_bins_new_10), mom_bins_new_10, miss_id_k_new_10)
gr_10.SetLineColor(ROOT.kRed + 3)
gr_10.SetMarkerColor(ROOT.kRed + 3)
gr_10.SetTitle("10 ps")
gr_10.Draw("PLsame")
gr_30 = ROOT.TGraph(len(mom_bins_new_30), mom_bins_new_30, miss_id_k_new_30)
gr_30.SetLineColor(ROOT.kRed + 2)
gr_30.SetMarkerColor(ROOT.kRed + 2)
gr_30.SetTitle("30 ps")
gr_30.Draw("PLsame")

gr_50 = ROOT.TGraph(len(mom_bins_new_50), mom_bins_new_50, miss_id_k_new_50)
gr_50.SetLineColor(ROOT.kRed + 1)
gr_50.SetMarkerColor(ROOT.kRed + 1)
gr_50.SetTitle("50 ps")
gr_50.Draw("PLsame")

gr_100 = ROOT.TGraph(len(mom_bins_new_100), mom_bins_new_100, miss_id_k_new_100)
gr_100.SetLineColor(ROOT.kRed)
gr_100.SetMarkerColor(ROOT.kRed)
gr_100.SetTitle("100 ps")
gr_100.Draw("PLsame")
input("wait")

gr_0 = ROOT.TGraph(len(mom_bins_new_0), mom_bins_new_0, miss_id_p_new_0)
gr_0.Draw("APL")
gr_0.SetLineColor(ROOT.kRed + 4)
gr_0.SetMarkerColor(ROOT.kRed + 4)
gr_0.SetTitle("0 ps")
gr_10 = ROOT.TGraph(len(mom_bins_new_10), mom_bins_new_10, miss_id_p_new_10)
gr_10.SetLineColor(ROOT.kRed + 3)
gr_10.SetMarkerColor(ROOT.kRed + 3)
gr_10.SetTitle("10 ps")
gr_10.Draw("PLsame")
gr_30 = ROOT.TGraph(len(mom_bins_new_30), mom_bins_new_30, miss_id_p_new_30)
gr_30.SetLineColor(ROOT.kRed + 2)
gr_30.SetMarkerColor(ROOT.kRed + 2)
gr_30.SetTitle("30 ps")
gr_30.Draw("PLsame")

gr_50 = ROOT.TGraph(len(mom_bins_new_50), mom_bins_new_50, miss_id_p_new_50)
gr_50.SetLineColor(ROOT.kRed + 1)
gr_50.SetMarkerColor(ROOT.kRed + 1)
gr_50.SetTitle("50 ps")
gr_50.Draw("PLsame")

gr_100 = ROOT.TGraph(len(mom_bins_new_100), mom_bins_new_100, miss_id_p_new_100)
gr_100.SetLineColor(ROOT.kRed)
gr_100.SetMarkerColor(ROOT.kRed)
gr_100.SetTitle("100 ps")
gr_100.Draw("PLsame")
input("wait")



# gr_pi_new = ROOT.TGraph(len(mom_bins_new), mom_bins_new, eff_pi_new)
# gr_pi_new.Draw("PLsame")
# gr_pi_new.SetLineColor(4)
# gr_pi_new.SetMarkerColor(6)
# gr_pi_new.SetTitle("#pi^{#pm} ID (new)")
# gr_k_new = ROOT.TGraph(len(mom_bins_new), mom_bins_new, eff_k_new)
# gr_k_new.SetLineColor(5)
# gr_k_new.SetMarkerColor(6)
# gr_k_new.SetTitle("K^{#pm} ID (new)")
# gr_k_new.Draw("PLsame")
# gr_p_new = ROOT.TGraph(len(mom_bins_new), mom_bins_new, eff_p_new)
# gr_p_new.SetLineColor(6)
# gr_p_new.SetMarkerColor(6)
# gr_p_new.SetTitle("p ID (new)")
# gr_p_new.Draw("PLsame")


input("wait")
