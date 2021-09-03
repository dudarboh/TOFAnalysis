import ROOT
import numpy as np
import time
import math
ROOT.EnableImplicitMT(2)

ROOT.gInterpreter.Declare('''

#define SPEED_OF_LIGHT 299.792458

''')

ROOT.gStyle.SetNdivisions(512)
ROOT.gStyle.SetNdivisions(512, "Y")
ROOT.gStyle.SetOptStat(10)

# Get true beta vs p curves

def get_curves(df, mom_str, method, n_mom_bins=50, to_draw=True):
    graphs = {}
    gr_diff = {}
    gr_sep = {}
    pdg_line = {}
    pdgs = [211, 321, 2212]
    m_pdg = {211 : 0.13957039, 321 : 0.493677, 2212 : 0.938272088}

    if to_draw:
        canvas = ROOT.TCanvas()
        canvas.Divide(2, 2)
        for i in range(1, 5):
            canvas.cd(i)
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()


    # 2D histo
    h_all = df.Histo2D(("h_all", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 1., 10., 200, -0.1, 1.3), mom_str, method)
    if to_draw:
        canvas.cd(1)
        ROOT.gPad.SetLogz()
        h_all.SetStats(1)
        h_all.Draw("colz")

    for k, pdg in enumerate(pdgs):
        if (pdg == 211):
            h = df.Filter("abs(pdg) == {}".format(pdg))\
            .Histo2D(("h", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 1., 10., 500, -0.1, 1.3), mom_str, method)
        elif (pdg == 321):
            h = df.Filter("abs(pdg) == {}".format(pdg))\
            .Histo2D(("h", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 1., 10., 500, -0.1, 1.3), mom_str, method)
        elif (pdg == 2212):
            h = df.Filter("abs(pdg) == {}".format(pdg))\
            .Histo2D(("h", "Method: {} ; p (GeV); mass (GeV)".format(method), n_mom_bins, 1., 10., 500, -0.1, 1.3), mom_str, method)


        graphs[pdg] = ROOT.TGraphErrors()
        graphs[pdg].SetMarkerStyle(20)
        graphs[pdg].SetMarkerColor(k+1)
        graphs[pdg].SetLineColor(k+1)

        gr_diff[pdg] = ROOT.TGraphErrors()
        gr_diff[pdg].SetMarkerStyle(20)
        gr_diff[pdg].SetMarkerColor(k+1)
        gr_diff[pdg].SetLineColor(k+1)
        if pdg == 211:
            graphs[pdg].SetTitle("#pi^{#pm}; p (GeV); mass (GeV)")
            gr_diff[pdg].SetTitle("#pi^{#pm}; p (GeV); mass - m_{PDG} (GeV)")
        elif pdg == 321:
            graphs[pdg].SetTitle("K^{#pm}; p (GeV); mass (GeV)")
            gr_diff[pdg].SetTitle("K^{#pm}; p (GeV); mass - m_{PDG} (GeV)")
        elif pdg == 2212:
            graphs[pdg].SetTitle("p; p (GeV); mass (GeV)")
            gr_diff[pdg].SetTitle("p; p (GeV); mass - m_{PDG} (GeV)")

        # *************Fit in Slices section***********
        points = 0
        mom_points = []
        for bin in range(1, n_mom_bins + 1):
            hp = h.ProjectionY("hp", bin, bin)
            if hp.GetEntries() < 10:
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
            graphs[pdg].GetXaxis().SetRangeUser(1., 10.)
            graphs[pdg].GetYaxis().SetRangeUser(-0.1, 1.3)
            gr_diff[pdg].GetXaxis().SetRangeUser(1., 10.)
            gr_diff[pdg].GetYaxis().SetRangeUser(-0.2, 0.2)
            canvas.cd(2)
            if pdg == 211:
                graphs[pdg].Draw("APE")
            pdg_line[pdg] = ROOT.TLine(1., m_pdg[pdg], 10., m_pdg[pdg])
            pdg_line[pdg].SetLineColor(4)
            pdg_line[pdg].Draw()
            graphs[pdg].Draw("PEsame")
            canvas.cd(4)
            gr_diff[pdg].Draw("APE0" if pdg == 211 else "PE0same")
    if to_draw:
        canvas.cd(4)
        ROOT.gPad.BuildLegend(0.25, 0.7, .35, .9)



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

        gr_sep["pik"].GetXaxis().SetRangeUser(1., 10.)
        gr_sep["pik"].GetYaxis().SetRangeUser(0., 5.)
    if to_draw:
        canvas.cd(3)
        gr_sep["pik"].Draw("APL")
        gr_sep["kp"].Draw("PLsame")
        ROOT.gPad.BuildLegend(0.7, 0.7, .9, .9)

        for i in range(1, 5):
            canvas.cd(i)
            ROOT.gPad.Update()

        # canvas.Print("method_{}_{}_ps.png".format(m, s))
        print("Finished with method", method)
        input("wait")

    return graphs, gr_diff, gr_sep




df = ROOT.RDataFrame("TOFAnalysis", "/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
df = df.Filter("n_ecal_hits > 0 && abs(ts_ecal_pos.z()) < 2385. && abs(ts_ecal_z0) < 1.")\
        .Define("mom_ip", "ts_ip_mom.r()")\
        .Define("mom_ecal", "ts_ecal_mom.r()")\
        .Define("mom_tanL", "sqrt(mom2_hm_TanL)")\
        .Define("mom_dz", "sqrt(mom2_hm_dz)")

df = df.Define("beta_ip", "track_length_ip/(tof_closest_0*SPEED_OF_LIGHT)")\
       .Define("beta_ecal", "track_length_calo/(tof_closest_0*SPEED_OF_LIGHT)")\
       .Define("mass_ip", "mom_ip / beta_ip * sqrt(1. - beta_ip*beta_ip)")\
       .Define("mass_ecal", "mom_ecal / beta_ecal * sqrt(1. - beta_ecal*beta_ecal)")\
       .Define("mass_tanL", "sqrt(2. * mom2_hm_TanL * (SPEED_OF_LIGHT*tof_closest_0/track_length_refit_tanL - 1.))")\
       .Define("mass_z", "sqrt(2. * mom2_hm_dz * (SPEED_OF_LIGHT*tof_closest_0/track_length_refit_z - 1.))")



# gr_ip, gr_diff_ip, gr_sep_ip = get_curves(df, "mom_ip", "mass_ip", n_mom_bins=50, to_draw=False)
# gr_ecal, gr_diff_ecal, gr_sep_ecal = get_curves(df, "mom_ecal", "mass_ecal", n_mom_bins=50, to_draw=False)
# gr_tanL, gr_diff_tanL, gr_sep_tanL = get_curves(df, "mom_tanL", "mass_tanL", n_mom_bins=50, to_draw=False)
# gr_z, gr_diff_z, gr_sep_z = get_curves(df, "mom_dz", "mass_z", n_mom_bins=50, to_draw=False)


def plot_ip_matrix():
    gr_ip, gr_diff_ip, gr_sep_ip = get_curves(df, "mom_ip", "mass_ip", n_mom_bins=50, to_draw=False)
    file = ROOT.TFile("/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
    tree = file.TOFAnalysis

    speed_of_light = 299.792458
    n_points = gr_ip[211].GetN()
    mom_bins = []
    for i in range(n_points):
        mom_bins.append(gr_ip[211].GetPointX(i))
    mom_bins = np.array(mom_bins)

    h = ROOT.TH2F("name", "title; x; y", 4, 0, 4, 4, 0, 4)
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

    for i, pfo in enumerate(tree):

        # if i > 200000:
            # break
        if i%10000 == 0:
            print("Event", i)

        if not (pfo.n_ecal_hits > 0 and abs(pfo.ts_ecal_pos.z()) < 2385. and abs(pfo.ts_ecal_z0) < 1.):
            continue

        pdg = pfo.pdg
        if abs(pdg) == 211:
            fillx = 0
        elif abs(pdg) == 321:
            fillx = 1
        elif abs(pdg) == 2212:
            fillx = 2
        else:
            fillx = 3
        mom_ip = pfo.ts_ip_mom.r()
        if mom_ip > 10. or mom_ip < 1.:
            continue
        min_idx = int(np.argmin( abs(mom_bins - mom_ip) ))

        tof_ip = pfo.tof_closest_0

        length_ip = pfo.track_length_ip
        beta_ip = length_ip / (tof_ip*speed_of_light)

        if 1. - beta_ip*beta_ip < 0:
            h.Fill(fillx, 3)
            continue
        mass_reco = mom_ip / beta_ip * math.sqrt(1. - beta_ip*beta_ip)


        prob_pi = 1. - math.erf( abs( mass_reco - gr_ip[211].Eval(mom_ip) ) / gr_ip[211].GetErrorY( min_idx ) /math.sqrt(2))
        prob_k = 1. - math.erf( abs( mass_reco - gr_ip[321].Eval(mom_ip) ) / gr_ip[321].GetErrorY( min_idx ) /math.sqrt(2))
        prob_p = 1. - math.erf( abs( mass_reco - gr_ip[2212].Eval(mom_ip) ) / gr_ip[2212].GetErrorY( min_idx )/math.sqrt(2))
        if prob_pi == 0 and prob_k == 0 and prob_p == 0:
            h.Fill(fillx, 3)
            continue
        filly = int( np.argmax(np.array([prob_pi, prob_k, prob_p])) )

        h.Fill(fillx, filly)


        # print("********new pfo************")
        # print("PDG", pdg)
        # print("Momentum", mom_ip)
        # print("Momentum bin", min_idx)
        # print("tof", tof_ip)
        # print("length", length_ip)
        # print("beta", beta_ip)
        # print("mass reco", mass_reco)
        # print("")
        # print("diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) ))
        # print("diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) ))
        # print("diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) ))
        # print("")
        # print("sigma_diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("")
        # print("Prob pi", prob_pi)
        # print("Prob k", prob_k)
        # print("Prob p", prob_p)
        # input("wait")
    h.Scale(1./h.GetEntries())
    h.Draw("colz text")

    eff = {}
    pur = {}
    names = ["total true pi", "total true K", "total true protons", "total true others"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(i, j)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  efficiency", h.GetBinContent(i, i)/n_total * 100)
        eff[names[i-1]] = h.GetBinContent(i, i)/n_total

    names = ["total reco pi", "total reco K", "total reco protons", "total undefined"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(j, i)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  purity", h.GetBinContent(i, i)/n_total * 100)
        pur[names[i-1]] = h.GetBinContent(i, i) / n_total


    input("wait")


def plot_ecal_matrix():
    gr, gr_diff, gr_sep = get_curves(df, "mom_ecal", "mass_ecal", n_mom_bins=50, to_draw=False)
    file = ROOT.TFile("/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
    tree = file.TOFAnalysis

    speed_of_light = 299.792458
    n_points = gr[211].GetN()
    mom_bins = []
    for i in range(n_points):
        mom_bins.append(gr[211].GetPointX(i))
    mom_bins = np.array(mom_bins)

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

    for i, pfo in enumerate(tree):

        # if i > 200000:
            # break
        if i%10000 == 0:
            print("Event", i)

        if not (pfo.n_ecal_hits > 0 and abs(pfo.ts_ecal_pos.z()) < 2385. and abs(pfo.ts_ecal_z0) < 1.):
            continue
        mom = pfo.ts_ecal_mom.r()
        if mom > 10. or mom < 1.:
            continue

        pdg = pfo.pdg
        if abs(pdg) == 211:
            fillx = 0
        elif abs(pdg) == 321:
            fillx = 1
        elif abs(pdg) == 2212:
            fillx = 2
        else:
            fillx = 3
        min_idx = int(np.argmin( abs(mom_bins - mom) ))

        tof = pfo.tof_closest_0

        length = pfo.track_length_calo
        beta = length / (tof*speed_of_light)

        if 1. - beta*beta < 0:
            h.Fill(fillx, 3)
            continue
        mass_reco = mom / beta * math.sqrt(1. - beta*beta)


        prob_pi = 1. - math.erf( abs( mass_reco - gr[211].Eval(mom) ) / gr[211].GetErrorY( min_idx ) /math.sqrt(2))
        prob_k = 1. - math.erf( abs( mass_reco - gr[321].Eval(mom) ) / gr[321].GetErrorY( min_idx ) /math.sqrt(2))
        prob_p = 1. - math.erf( abs( mass_reco - gr[2212].Eval(mom) ) / gr[2212].GetErrorY( min_idx )/math.sqrt(2))
        if prob_pi == 0 and prob_k == 0 and prob_p == 0:
            h.Fill(fillx, 3)
            continue
        filly = int( np.argmax(np.array([prob_pi, prob_k, prob_p])) )

        h.Fill(fillx, filly)


        # print("********new pfo************")
        # print("PDG", pdg)
        # print("Momentum", mom_ip)
        # print("Momentum bin", min_idx)
        # print("tof", tof_ip)
        # print("length", length_ip)
        # print("beta", beta_ip)
        # print("mass reco", mass_reco)
        # print("")
        # print("diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) ))
        # print("diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) ))
        # print("diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) ))
        # print("")
        # print("sigma_diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("")
        # print("Prob pi", prob_pi)
        # print("Prob k", prob_k)
        # print("Prob p", prob_p)
        # input("wait")
    h.Scale(1./h.GetEntries())
    h.Draw("colz text")

    eff = {}
    pur = {}
    names = ["total true pi", "total true K", "total true protons", "total true others"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(i, j)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  efficiency", h.GetBinContent(i, i)/n_total * 100)
        eff[names[i-1]] = h.GetBinContent(i, i)/n_total

    names = ["total reco pi", "total reco K", "total reco protons", "total undefined"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(j, i)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  purity", h.GetBinContent(i, i)/n_total * 100)
        pur[names[i-1]] = h.GetBinContent(i, i) / n_total


    input("wait")


def plot_tanL_matrix():
    gr, gr_diff, gr_sep = get_curves(df, "mom_tanL", "mass_tanL", n_mom_bins=50, to_draw=False)
    file = ROOT.TFile("/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
    tree = file.TOFAnalysis

    speed_of_light = 299.792458
    n_points = gr[211].GetN()
    mom_bins = []
    for i in range(n_points):
        mom_bins.append(gr[211].GetPointX(i))
    mom_bins = np.array(mom_bins)

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

    for i, pfo in enumerate(tree):

        # if i > 200000:
            # break
        if i%10000 == 0:
            print("Event", i)

        if not (pfo.n_ecal_hits > 0 and abs(pfo.ts_ecal_pos.z()) < 2385. and abs(pfo.ts_ecal_z0) < 1.):
            continue
        mom = math.sqrt(pfo.mom2_hm_TanL)
        if mom > 10. or mom < 1.:
            continue

        pdg = pfo.pdg
        if abs(pdg) == 211:
            fillx = 0
        elif abs(pdg) == 321:
            fillx = 1
        elif abs(pdg) == 2212:
            fillx = 2
        else:
            fillx = 3
        min_idx = int(np.argmin( abs(mom_bins - mom) ))

        tof = pfo.tof_closest_0
        length = pfo.track_length_refit_tanL


        if 2. * pfo.mom2_hm_TanL * (speed_of_light*tof/length - 1.) < 0:
            h.Fill(fillx, 3)
            continue
        mass_reco = math.sqrt(2. * pfo.mom2_hm_TanL * (speed_of_light*tof/length - 1.))


        prob_pi = 1. - math.erf( abs( mass_reco - gr[211].Eval(mom) ) / gr[211].GetErrorY( min_idx ) /math.sqrt(2))
        prob_k = 1. - math.erf( abs( mass_reco - gr[321].Eval(mom) ) / gr[321].GetErrorY( min_idx ) /math.sqrt(2))
        prob_p = 1. - math.erf( abs( mass_reco - gr[2212].Eval(mom) ) / gr[2212].GetErrorY( min_idx )/math.sqrt(2))
        if prob_pi == 0 and prob_k == 0 and prob_p == 0:
            h.Fill(fillx, 3)
            continue
        filly = int( np.argmax(np.array([prob_pi, prob_k, prob_p])) )

        h.Fill(fillx, filly)


        # print("********new pfo************")
        # print("PDG", pdg)
        # print("Momentum", mom_ip)
        # print("Momentum bin", min_idx)
        # print("tof", tof_ip)
        # print("length", length_ip)
        # print("beta", beta_ip)
        # print("mass reco", mass_reco)
        # print("")
        # print("diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) ))
        # print("diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) ))
        # print("diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) ))
        # print("")
        # print("sigma_diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("")
        # print("Prob pi", prob_pi)
        # print("Prob k", prob_k)
        # print("Prob p", prob_p)
        # input("wait")
    h.Scale(1./h.GetEntries())
    h.Draw("colz text")

    eff = {}
    pur = {}
    names = ["total true pi", "total true K", "total true protons", "total true others"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(i, j)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  efficiency", h.GetBinContent(i, i)/n_total * 100)
        eff[names[i-1]] = h.GetBinContent(i, i)/n_total

    names = ["total reco pi", "total reco K", "total reco protons", "total undefined"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(j, i)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  purity", h.GetBinContent(i, i)/n_total * 100)
        pur[names[i-1]] = h.GetBinContent(i, i) / n_total


    input("wait")


def plot_dz_matrix():
    gr, gr_diff, gr_sep = get_curves(df, "mom_dz", "mass_z", n_mom_bins=50, to_draw=False)
    file = ROOT.TFile("/nfs/dust/ilc/user/dudarboh/final_files/SET/final.root")
    tree = file.TOFAnalysis

    speed_of_light = 299.792458
    n_points = gr[211].GetN()
    mom_bins = []
    for i in range(n_points):
        mom_bins.append(gr[211].GetPointX(i))
    mom_bins = np.array(mom_bins)

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

    for i, pfo in enumerate(tree):

        # if i > 200000:
            # break
        if i%10000 == 0:
            print("Event", i)

        if not (pfo.n_ecal_hits > 0 and abs(pfo.ts_ecal_pos.z()) < 2385. and abs(pfo.ts_ecal_z0) < 1.):
            continue
        mom = math.sqrt(pfo.mom2_hm_dz)
        if mom > 10. or mom < 1.:
            continue

        pdg = pfo.pdg
        if abs(pdg) == 211:
            fillx = 0
        elif abs(pdg) == 321:
            fillx = 1
        elif abs(pdg) == 2212:
            fillx = 2
        else:
            fillx = 3
        min_idx = int(np.argmin( abs(mom_bins - mom) ))

        tof = pfo.tof_closest_0
        length = pfo.track_length_refit_z


        if 2. * pfo.mom2_hm_dz * (speed_of_light*tof/length - 1.) < 0:
            h.Fill(fillx, 3)
            continue
        mass_reco = math.sqrt(2. * pfo.mom2_hm_dz * (speed_of_light*tof/length - 1.))


        prob_pi = 1. - math.erf( abs( mass_reco - gr[211].Eval(mom) ) / gr[211].GetErrorY( min_idx ) /math.sqrt(2))
        prob_k = 1. - math.erf( abs( mass_reco - gr[321].Eval(mom) ) / gr[321].GetErrorY( min_idx ) /math.sqrt(2))
        prob_p = 1. - math.erf( abs( mass_reco - gr[2212].Eval(mom) ) / gr[2212].GetErrorY( min_idx )/math.sqrt(2))
        if prob_pi == 0 and prob_k == 0 and prob_p == 0:
            h.Fill(fillx, 3)
            continue
        filly = int( np.argmax(np.array([prob_pi, prob_k, prob_p])) )

        h.Fill(fillx, filly)


        # print("********new pfo************")
        # print("PDG", pdg)
        # print("Momentum", mom_ip)
        # print("Momentum bin", min_idx)
        # print("tof", tof_ip)
        # print("length", length_ip)
        # print("beta", beta_ip)
        # print("mass reco", mass_reco)
        # print("")
        # print("diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) ))
        # print("diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) ))
        # print("diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) ))
        # print("")
        # print("sigma_diff to pion", abs( mass_reco - gr_ip[211].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to kaon", abs( mass_reco - gr_ip[321].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("sigma_diff to proton", abs( mass_reco - gr_ip[2212].Eval(mom_ip) )/ gr_ip[211].GetErrorY( min_idx ) )
        # print("")
        # print("Prob pi", prob_pi)
        # print("Prob k", prob_k)
        # print("Prob p", prob_p)
        # input("wait")
    h.Scale(1./h.GetEntries())
    h.Draw("colz text")

    eff = {}
    pur = {}
    names = ["total true pi", "total true K", "total true protons", "total true others"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(i, j)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  efficiency", h.GetBinContent(i, i)/n_total * 100)
        eff[names[i-1]] = h.GetBinContent(i, i)/n_total

    names = ["total reco pi", "total reco K", "total reco protons", "total undefined"]
    for i in range(1, 5):
        n_total = 0.
        for j in range(1, 5):
            n_total += h.GetBinContent(j, i)
        print(names[i-1], h.GetBinContent(i, i),"/",n_total, "  purity", h.GetBinContent(i, i)/n_total * 100)
        pur[names[i-1]] = h.GetBinContent(i, i) / n_total


    input("wait")


plot_dz_matrix()
