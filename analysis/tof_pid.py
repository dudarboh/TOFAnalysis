import ROOT
import numpy as np
import time
ROOT.EnableImplicitMT(2)

ROOT.gInterpreter.Declare('''

#include "Math/ProbFunc.h"
#define SPEED_OF_LIGHT 299.792458

double prob_pid(double mom, double beta, TGraphErrors gr){
    std::vector<double> mom_arr {1.15, 1.45, 1.75, 2.05, 2.35, 2.65,
    2.95, 3.25, 3.55, 3.85, 4.15, 4.45,
    4.75, 5.05, 5.35, 5.65, 5.95, 6.25,
    6.55, 6.85, 7.15, 7.45, 7.75, 8.05,
    8.35, 8.65, 8.95, 9.25, 9.55, 9.85};

    int closest_idx = std::min_element( mom_arr.begin(), mom_arr.end(), [&mom](double a, double b) { return std::abs(a - mom) < std::abs(b - mom) ;} ) - mom_arr.begin();

    cout<<"*************************"<<endl;
    cout<<"Momentum: "<<mom<<endl;
    cout<<"Closest mom:"<<mom_arr[closest_idx]<<endl;
    cout<<"Beta:"<<beta<<endl;
    cout<<"True Beta:"<<gr.Eval(mom)<<endl;


    // distance to the pion line in sigmas
    double d_to_line = std::abs( beta - gr.Eval(mom) ) / gr.GetErrorY(closest_idx);
    cout<<"Delta beta: "<< std::abs( beta - gr.Eval(mom) )<<endl;
    cout<<"Sigma:      "<< gr.GetErrorY(closest_idx)<<endl;
    cout<<"d in sigmas :"<< d_to_line<<endl;
    double prob = std::abs(ROOT::Math::normal_cdf(d_to_line) - ROOT::Math::normal_cdf( - d_to_line) ) * 100.;
    std::cout<<prob<<endl;
    return prob;
}

''')

ROOT.gStyle.SetNdivisions(512)
# ROOT.gStyle.SetNdivisions(512, "Y")

# Get true beta vs p curves
def get_true_curves():
    f_true = {}
    masses = [0.13957039, 0.493677, 0.938272088]
    pdgs = [211, 321, 2212]


    # values to check what happens to sep power plots with perferct dl
    l0 = 2000. # mm
    dt = 0.010 # ns
    n_bins = 30
    canvas = ROOT.TCanvas()
    true_graphs = {}
    h = df.Histo1D(("h", "title", n_bins, 1., 10.), "mom")
    for k, (mass, pdg) in enumerate( zip(masses, pdgs)):
        true_graphs[pdg] = ROOT.TGraphErrors()
        f_true[pdg] = ROOT.TF1("f_true_{}".format(pdg), "sqrt(1. - {0}*{0}/(x*x) )".format(mass), 1., 10.)
        f_true[pdg].SetNpx(1000)
        f_true[pdg].SetLineStyle(2)
        f_true[pdg].SetLineColor(k+1)

        for bin in range(1, n_bins + 1):
            mom = h.GetXaxis().GetBinCenter( bin )
            beta = f_true[pdg].Eval(mom)
            time = l0 / (beta*299.792458)
            true_graphs[pdg].SetPoint( bin-1, mom, beta )
            true_graphs[pdg].SetPointError( bin-1, 0., beta * dt / time )

    gr_true_sep = {}
    gr_true_sep["pik"] = ROOT.TGraphErrors()
    gr_true_sep["pik"].SetMarkerStyle(20)
    gr_true_sep["pik"].SetMarkerColor(12)
    gr_true_sep["pik"].SetLineColor(12)
    gr_true_sep["pik"].SetTitle("#pi/K; p (GeV); Sep. power")
    gr_true_sep["kp"] = ROOT.TGraphErrors()
    gr_true_sep["kp"].SetMarkerStyle(20)
    gr_true_sep["kp"].SetMarkerColor(46)
    gr_true_sep["kp"].SetLineColor(46)
    gr_true_sep["kp"].SetTitle("K/p; p (GeV); Sep. power")

    for bin in range(1, n_bins + 1):
        mom = h.GetXaxis().GetBinCenter( bin )
        mean_pi = true_graphs[211].GetPointY(bin -1)
        sigma_pi = true_graphs[211].GetErrorY(bin - 1)
        mean_k = true_graphs[321].GetPointY(bin -1)
        sigma_k = true_graphs[321].GetErrorY(bin - 1)
        mean_p = true_graphs[2212].GetPointY(bin -1)
        sigma_p = true_graphs[2212].GetErrorY(bin - 1)

        sep_power_pik = abs(mean_pi - mean_k) / np.sqrt( (sigma_pi*sigma_pi + sigma_k*sigma_k)/2. )
        sep_power_kp = abs(mean_k - mean_p) / np.sqrt( (sigma_k*sigma_k + sigma_p*sigma_p)/2. )

        gr_true_sep["pik"].SetPoint( bin-1, mom , sep_power_pik )
        gr_true_sep["kp"].SetPoint( bin-1, mom , sep_power_kp )

        gr_true_sep["pik"].GetYaxis().SetRangeUser(0., 5.)
        gr_true_sep["kp"].GetYaxis().SetRangeUser(0., 5.)

    gr_true_sep["pik"].Draw("APL")
    gr_true_sep["kp"].Draw("PLsame")
    ROOT.gPad.BuildLegend(0.7, 0.7, .9, .9)
    # input("wait")

    return f_true


def get_curves(df, m, s, f_true, n_mom_bins=30, to_draw=True):
    graphs = {}
    gr_diff = {}
    gr_sep = {}
    pdgs = [211, 321, 2212]

    if to_draw:
        canvas = ROOT.TCanvas()
        canvas.Divide(2, 2)
        for i in range(1, 5):
            canvas.cd(i)
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()

    # 2D histo
    h_all = df.Filter("abs(pdg) == 211")\
              .Histo2D(("h_all", "Method: {} {} ps; p (GeV);#beta".format(m, s), n_mom_bins, 1., 10., 2000, 0.7, 1.05), "mom", "beta_{}_{}".format(m, s))
    if to_draw:
        canvas.cd(1)
        ROOT.gPad.SetLogz()
        h_all.SetStats(0)
        h_all.Draw("colz")

    for k, pdg in enumerate(pdgs):
        h = df.Filter("abs(pdg) == {}".format(pdg))\
        .Histo2D(("h", "Method: {} {} ps; p (GeV);#beta".format(m, s), n_mom_bins, 1., 10., 2000, 0.7, 1.05), "mom", "beta_{}_{}".format(m, s))

        graphs[pdg] = ROOT.TGraphErrors()
        graphs[pdg].SetMarkerStyle(20)
        graphs[pdg].SetMarkerColor(k+1)
        graphs[pdg].SetLineColor(k+1)

        gr_diff[pdg] = ROOT.TGraphErrors()
        gr_diff[pdg].SetMarkerStyle(20)
        gr_diff[pdg].SetMarkerColor(k+1)
        gr_diff[pdg].SetLineColor(k+1)
        if pdg == 211:
            graphs[pdg].SetTitle("#pi^{#pm}; p (GeV); #beta")
            gr_diff[pdg].SetTitle("#pi^{#pm}; p (GeV); #beta - #beta_{true}")
        elif pdg == 321:
            graphs[pdg].SetTitle("K^{#pm}; p (GeV); #beta")
            gr_diff[pdg].SetTitle("K^{#pm}; p (GeV); #beta - #beta_{true}")
        elif pdg == 2212:
            graphs[pdg].SetTitle("p; p (GeV); #beta")
            gr_diff[pdg].SetTitle("p; p (GeV); #beta - #beta_{true}")

        # *************Fit in Slices section***********
        for bin in range(1, n_mom_bins + 1):
            hp = h.ProjectionY("hp", bin, bin)
            hp.SetStats(0)
            max_bin = hp.GetXaxis().GetBinCenter( hp.GetMaximumBin() )

            func = ROOT.TF1("func", "gaus", max_bin-0.05, max_bin+0.05)
            func.SetParameter(1, max_bin)
            func.SetParLimits(1, max_bin-0.05, max_bin+0.05)
            func.SetParameter(2, 0.03)
            func.SetParLimits(2, 0., 0.1)
            func.SetNpx(1000)
            hp.Fit(func, "QRN")
            mom =  h.GetXaxis().GetBinCenter( bin )
            mean_fit = func.GetParameter(1)
            sigma_fit = func.GetParameter(2)
            graphs[pdg].SetPoint( bin-1, mom, mean_fit )
            graphs[pdg].SetPointError( bin-1, 0., sigma_fit )

            gr_diff[pdg].SetPoint( bin-1, mom, mean_fit - f_true[pdg].Eval(mom) )
            gr_diff[pdg].SetPointError( bin-1, 0., sigma_fit )

            graphs[pdg].GetXaxis().SetRangeUser(1., 10.)
            graphs[pdg].GetYaxis().SetRangeUser(0.7, 1.05)
            gr_diff[pdg].GetXaxis().SetRangeUser(1., 10.)
            gr_diff[pdg].GetYaxis().SetRangeUser(-0.02, 0.02)
        if to_draw:
            canvas.cd(2)
            graphs[pdg].Draw("APE" if pdg == 211 else "PEsame")
            f_true[pdg].Draw("same")

            canvas.cd(4)
            gr_diff[pdg].Draw("APE" if pdg == 211 else "PEsame")
            ROOT.gPad.BuildLegend(0.7, 0.7, .9, .9)


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

    for bin in range(1, n_mom_bins + 1):
        mom = h.GetXaxis().GetBinCenter( bin )
        mean_pi = graphs[211].GetPointY(bin -1)
        sigma_pi = graphs[211].GetErrorY(bin - 1)
        mean_k = graphs[321].GetPointY(bin -1)
        sigma_k = graphs[321].GetErrorY(bin - 1)
        mean_p = graphs[2212].GetPointY(bin -1)
        sigma_p = graphs[2212].GetErrorY(bin - 1)

        sep_power_pik = abs(mean_pi - mean_k) / np.sqrt( (sigma_pi*sigma_pi + sigma_k*sigma_k)/2. )
        sep_power_kp = abs(mean_k - mean_p) / np.sqrt( (sigma_k*sigma_k + sigma_p*sigma_p)/2. )

        gr_sep["pik"].SetPoint( bin-1, mom , sep_power_pik )
        gr_sep["kp"].SetPoint( bin-1, mom , sep_power_kp )

        gr_sep["pik"].GetYaxis().SetRangeUser(0., 5.)
        gr_sep["kp"].GetYaxis().SetRangeUser(0., 5.)
    if to_draw:
        canvas.cd(3)
        gr_sep["pik"].Draw("APL")
        gr_sep["kp"].Draw("PLsame")
        ROOT.gPad.BuildLegend(0.7, 0.7, .9, .9)

        for i in range(1, 5):
            canvas.cd(i)
            ROOT.gPad.Update()

        # canvas.Print("method_{}_{}_ps.png".format(m, s))
        print("Finished with method", m, s, "ps")
        input("wait")

    return graphs, gr_diff, gr_sep



# Variables for the loops
smearings = [0 , 10, 30, 50, 100]
methods =["set", "set_avg", "closest", "fastest", "frank_fit", "frank_avg"]
pdgs = [211, 321, 2212]

# Filter and calculate betas
df = ROOT.RDataFrame("SETAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
df = df.Filter("n_ecal_hits > 0 && abs(pos_closest.z()) < 2000.")\
       .Define("mom", "ts_calo_mom.r()")

for s in smearings:
    df = df.Define("beta_set_{}".format(s), "track_length_set/(tof_set_front_{}*SPEED_OF_LIGHT)".format(s))\
           .Define("beta_closest_{}".format(s), "track_length_calo/(tof_closest_{}*SPEED_OF_LIGHT)".format(s))\
           .Define("beta_fastest_{}".format(s), "track_length_calo/(tof_fastest_{}*SPEED_OF_LIGHT)".format(s))\
           .Define("beta_frank_fit_{}".format(s), "track_length_calo/(tof_frank_fit_{}*SPEED_OF_LIGHT)".format(s))\
           .Define("beta_frank_avg_{}".format(s), "track_length_calo/(tof_frank_avg_{}*SPEED_OF_LIGHT)".format(s))\
           .Define("beta_set_avg_{}".format(s), "track_length_set/(0.5*( tof_set_front_{0}+tof_set_back_{0} )*SPEED_OF_LIGHT)".format(s))\



canvas2 = ROOT.TCanvas()
canvas2.SetGridx()
canvas2.SetGridy()

# n_pions = df.Filter("abs(pdg) == 211").Count().GetValue()
# n_kaons = df.Filter("abs(pdg) == 321").Count().GetValue()
# n_protons = df.Filter("abs(pdg) == 2212").Count().GetValue()
# n_tot = n_pions + n_kaons + n_protons
# print(n_pions/n_tot*100., n_kaons/n_tot*100., n_protons/n_tot*100.)
# h = ROOT.TH1F("h", "", 3, 0, 3)
# h.Fill("#pi^{#pm}", n_pions)
# h.Fill("K^{#pm}", n_kaons)
# h.Fill("protons", n_protons)

# h = df.Histo1D(("h_tr_len", "title;x;y", 1000, 1750, 4000), "track_length_calo")
# h.Draw("histo")
# input("wait")

f_true = get_true_curves()
graphs = {}
gr_diff = {}
gr_sep = {}

for s in smearings:
    for m_idx, m in enumerate(methods):
        graphs[m, s], gr_diff[m, s], gr_sep[m, s] = get_curves(df, m , s, f_true, 30, True)

        #     canvas2.cd()
        #     gr_sep[m, s]["pik"].Draw("APL" if m_idx == 0 else "PLsame")
        #     gr_sep[m, s]["pik"].SetLineColor(m_idx+1)
        #     gr_sep[m, s]["pik"].SetMarkerColor(m_idx+1)
        #     gr_sep[m, s]["pik"].SetTitle(m)
        #     canvas2.Update()
        # canvas2.BuildLegend()

        # ROOT.gInterpreter.Declare("""
        #     TGraphErrors gr_{0}_{1}_pi = TGraphErrors(TPython::Exec("graphs[ '{0}' , {1}][211] ") );
        #     TGraphErrors gr_{0}_{1}_k = TGraphErrors(TPython::Exec("graphs[ '{0}' , {1}][321] ") );
        #     TGraphErrors gr_{0}_{1}_p = TGraphErrors(TPython::Exec("graphs[ '{0}' , {1}][2212] ") );
        # """.format(m, s)         )
        # df = df.Define("prob_{}_{}_pi".format(m, s), "prob_pid(mom, beta_{0}_{1}, gr_{0}_{1}_pi)".format(m, s))
        # test = df.Histo1D(("test", "title", 500, 0., 1.), "prob_{}_{}_pi".format(m, s))
        # test.Draw()

        input("wait")
