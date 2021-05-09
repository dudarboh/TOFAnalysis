import ROOT
import numpy as np
from array import array

ROOT.EnableImplicitMT(2)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetOptStat(1110)


def time_res_bias(m_true=139.5701835353535353535, l_track=2000., momentum=1000., resolution=0.01):
    c_light = 299.792458
    t_true = l_track / c_light * np.sqrt(1 + (m_true / momentum)**2) # (6.7359462 ns)
    func_str = "{0}*x/(sqrt(x^2+{1}^2)*sqrt(2*pi)*{2})*exp(-({0}*sqrt(x^2+{1}^2) - {3})^2/(2*{2}^2))".format(l_track / (momentum * c_light), momentum, resolution, t_true)
    f = ROOT.TF1("f", func_str, 0., 500.)
    f.SetNpx(1000000)
    return f.GetMaximumX() - m_true


def len_res_bias(m_true=139.5701835353535353535, l_track=2000., momentum=1000., resolution=0.01):
    c_light = 299.792458
    t_true = l_track / c_light * np.sqrt(1 + (m_true / momentum)**2) # (6.7359462 ns)
    func_str = "{0}*x/((x^2+{1}^2)*sqrt(x^2+{1}^2)*sqrt(2*pi)*{2})*exp(-({0}/sqrt(x^2+{1}^2) - {3})^2/(2*{2}^2))".format(t_true * c_light * momentum, momentum, resolution, l_track)
    f = ROOT.TF1("f", func_str, 0., 500.)
    f.SetNpx(1000000)
    return f.GetMaximumX() - m_true


def smear_bias_vs_mom():
    canvas = ROOT.TCanvas()
    # for time resolution plots
    res = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]  # ns
    # for track length resolution plots
    # res = [0.01, 0.1, 1., 3., 5.]  # mm

    graphs = [ROOT.TGraph() for i in range(len(res))]
    mom = np.linspace(0.8, 5., 100)*1000
    for i, s in enumerate(res):
        for j, m in enumerate(mom):
            print("mom", m, "smear", s)
            bias = time_res_bias(l_track=2500., momentum=m, resolution=s)
            # bias = len_res_bias(l_track=2500., momentum=m, resolution=s)
            graphs[i].SetPoint(j, m, bias)

    mg = ROOT.TMultiGraph()
    for i, g in enumerate(graphs):
        mg.Add(g)
        g.SetLineColor(i+1)
        g.SetLineWidth(2)
        g.SetTitle("{0} ps; momentum, [MeV]; #Delta m, [MeV]".format(int(res[i]*1000)))
        # g.SetTitle("{0} um; momentum, [MeV]; #Delta m, [MeV]".format(int(res[i]*1000)))

    mg.Draw("AL")
    canvas.BuildLegend()

    canvas.Update()
    input("wait")



def math_time_resolution():
    canvas = ROOT.TCanvas()
    c_light = 299.792458
    m_pi = 139.5701835353535353535  # MeV
    l_track = 2000.  # mm
    momentum =  10000.  # MeV
    t_true = l_track / c_light * np.sqrt(1 + (m_pi / momentum)**2) # (6.7359462 ns)
    time_resolution = [0.001, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3]  # ns
    mass_str = ["{0}*x/(sqrt(x^2+{1}^2)*sqrt(2*pi)*{2})*exp(-({0}*sqrt(x^2+{1}^2) - {3})^2/(2*{2}^2))".format(l_track / (momentum * c_light), momentum, res, t_true) for res in time_resolution]
    funcs = [ROOT.TF1("f{}".format(i), m, 0., 500.) for i, m in enumerate(mass_str)]

    for i, f in enumerate(funcs):
        f.SetNpx(1000)
        f.SetLineColor(i+1)
        f.SetLineWidth(2)
        f.SetTitle("{0} ps time resolution, Mass peak shift {1} MeV".format(int(time_resolution[i]*1000), round((f.GetMaximumX() - m_pi), 4)))
        f.Draw("" if i == 0 else "same")


    x = np.array([139.570183535, 493.677, 938.2720881629])
    line_pion = ROOT.TLine(x[0], 0., x[0], 1000.)
    line_pion.SetLineStyle(9)
    line_pion.SetLineWidth(3)

    canvas.BuildLegend()
    line_pion.Draw()
    canvas.Update()
    input("wait")


def math_length_resolution():
    canvas = ROOT.TCanvas()
    c_light = 299.792458
    m_pi = 139.5701835353535353535  # MeV
    tof = 6.7359462 # ns
    momentum = 1000.  # MeV
    l_true = tof * c_light * momentum / np.sqrt(m_pi*m_pi + momentum*momentum)  # mm (2000 mm)
    l_resolution = [0.01, 0.1, 1., 3., 5.]  # mm
    mass_str = ["{0}*x/((x^2+{1}^2)*sqrt(x^2+{1}^2)*sqrt(2*pi)*{2})*exp(-({0}/sqrt(x^2+{1}^2) - {3})^2/(2*{2}^2))".format(tof * c_light * momentum, momentum, res, l_true) for res in l_resolution]
    funcs = [ROOT.TF1("f{}".format(i), m, 0., 500.) for i, m in enumerate(mass_str)]

    for i, f in enumerate(funcs):
        f.SetNpx(10000)
        f.SetLineColor(i+1)
        f.SetLineWidth(2)
        f.SetTitle("{0} um track length resolution, Mass peak shift {1} MeV".format(int(l_resolution[i]*1000), round((f.GetMaximumX() - m_pi), 4)))
        f.Draw("" if i == 0 else "same")


    x = np.array([139.570183535, 493.677, 938.2720881629])
    line_pion = ROOT.TLine(x[0], 0., x[0], 1000.)
    line_pion.SetLineStyle(9)
    line_pion.SetLineWidth(3)

    canvas.BuildLegend()
    line_pion.Draw()
    canvas.Update()
    input("wait")


def math_mom_resolution():
    canvas = ROOT.TCanvas()
    c_light = 299.792458 # mm/ns
    m_pi = 139.5701835353535353535  # MeV
    tof = 6.7359462 # ns
    length = 2000.  # mm (1 GeV)
    mom_true = m_pi / np.sqrt((tof*c_light/length)**2 - 1)  # MeV
    mom_resolution = [0.01, 0.02, 0.05, 0.1, 1.]  # MeV
    mass_str = ["{0}/(sqrt(2*pi)*{1})*exp(-({0}*x - {2})^2/(2*{1}^2))".format(1. / np.sqrt((tof*c_light/length)**2 - 1), res, mom_true) for res in mom_resolution]
    funcs = [ROOT.TF1("f{}".format(i), m, 0., 500.) for i, m in enumerate(mass_str)]

    for i, f in enumerate(funcs):
        f.SetNpx(100000)
        f.SetLineColor(i+1)
        f.SetLineWidth(2)
        f.SetTitle("{0} keV track length resolution, Mass peak shift {1} MeV".format(int(mom_resolution[i]*1000), round((f.GetMaximumX() - m_pi), 4)))
        f.Draw("" if i == 0 else "same")


    x = np.array([139.570183535, 493.677, 938.2720881629])
    line_pion = ROOT.TLine(x[0], 0., x[0], 1000.)
    line_pion.SetLineStyle(9)
    line_pion.SetLineWidth(3)

    canvas.BuildLegend()
    line_pion.Draw()
    canvas.Update()
    input("wait")


def tof_resolution():
    canvas = ROOT.TCanvas()

    x = array("d", [0, 10, 30, 50, 100, 200, 300])
    std_dev_closest = np.array([0.04426776784226744, 0.045369404059289506, 0.053212150877549345, 0.06676857388766645, 0.10911738826687646, 0.2050308280317226, 0.3031871280341813])*1000
    std_dev_closest = array("d", std_dev_closest)
    gr_closest = ROOT.TGraph(len(x), x, std_dev_closest)
    gr_closest.SetTitle("Closest; ECAL hit resolution, [ps]; TOF_{reco} resolution, [ps]")
    gr_closest.SetLineColor(1)
    gr_closest.Draw("P")

    std_dev_fastest = np.array([0.016595086385795942, 0.01980792944166161, 0.030630083143368857, 0.04221194400365668, 0.07299083139913864, 0.13558136135349277, 0.19941955902618647])*1000
    std_dev_fastest = array("d", std_dev_fastest)
    gr_fastest = ROOT.TGraph(len(x), x, std_dev_fastest)
    gr_fastest.SetTitle("Fastest; ECAL hit resolution; TOF_{reco} resolution")
    gr_fastest.SetLineColor(2)
    gr_fastest.Draw("same")

    std_dev_frank = np.array([0.09994011369599361, 0.10506129344351821, 0.13096954444819614, 0.16713110169447057, 0.2426765213961192, 0.3507702741227225, 0.46157181713825324])*1000
    std_dev_frank = array("d", std_dev_frank)
    gr_frank = ROOT.TGraph(len(x), x, std_dev_frank)
    gr_frank.SetTitle("Frank; ECAL hit resolution; TOF_{reco} resolution")
    gr_frank.SetLineColor(3)
    gr_frank.Draw("same")

    std_dev_cyl = np.array([0.045420428142173085, 0.06247360765592981, 0.11726147541634689, 0.1709884477097624, 0.28043845580795507, 0.42770916866506364, 0.5554273182284295])*1000
    std_dev_cyl = array("d", std_dev_cyl)
    gr_cyl = ROOT.TGraph(len(x), x, std_dev_cyl)
    gr_cyl.SetTitle("5 mm Cylinder; ECAL hit resolution; TOF_{reco} resolution")
    gr_cyl.SetLineColor(4)
    gr_cyl.Draw("same")

    canvas.BuildLegend()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.Update()
    input("wait")

# math_time_resolution()
# math_length_resolution()
# math_mom_resolution()
smear_bias_vs_mom()
