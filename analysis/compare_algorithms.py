import ROOT

def compare1():
    methods = ["Frank's average",
                "Closest",
                "Fit all hits",
                "#splitline{Fit all hits}{refit if #chi^{2}>2}",
                "#splitline{Fit closest to line}{in 10 layers}",
                "#splitline{Fit 4 mm cyl}{in 10 layers}",
                "#splitline{Fit 4 mm cyl}{in 30 layers}",
                "#splitline{Fit 15 mm cyl}{in 10 layers}",
                "#splitline{Fit 15 mm cyl}{#splitline{in 10 layers}{refit if #chi^{2} > 2}}"


                ]

    means_histo = [0.005359, 0.001287, 0.003334, 0.001595, -0.002332, -0.0004862, -0.0008527, 0.003215, 0.00007909]
    sigmas_histo = [0.006956, 0.004817, 0.01296, 0.009424, -0.009301, 0.006654, 0.006453, 0.01127, 0.007433]

    means_fit = [1.55624e-03, 6.14139e-04, 3.74024e-03, 3.55577e-03, -8.67148e-04, -2.18261e-04, -4.08056e-04, 1.51402e-03, 1.08629e-03]
    sigmas_fit = [1.68825e-03, 2.69909e-03, 7.95629e-03, 6.57819e-03, 2.47641e-03, 2.27118e-03, 2.38832e-03, 5.23011e-03, 4.49907e-03]


    # methods = ["Frank's average",
    #             "Closest",
    #             "Fit all hits",
    #             "#splitline{Fit closest to line}{in 10 layers}",
    #             "#splitline{Fit 4 mm cyl}{in 10 layers}",
    #             "#splitline{Fit 4 mm cyl}{in 30 layers}"]
    #
    # means_histo = [0.005359, 0.001287, 0.003334, -0.002332, -0.0004862, -0.0008527]
    # sigmas_histo = [0.006956, 0.004817, 0.01296, -0.009301, 0.006654, 0.006453]
    #
    # means_fit = [1.55624e-03, 6.14139e-04, 3.74024e-03, -8.67148e-04, -2.18261e-04, -4.08056e-04]
    # sigmas_fit = [1.68825e-03, 2.69909e-03, 7.95629e-03, 2.47641e-03, 2.27118e-03, 2.38832e-03]

    n_pfo_survive = ["1", "2", "3", "4", "5", "6"]

    print("debug2")

    n = len(means_histo)

    canvas = ROOT.TCanvas()
    h1 = ROOT.TH1F("h1", "Histo values", n, 0, n)
    h2 = ROOT.TH1F("h2", "Fit values", n, 0, n)
    print("debug3")
    latex = ROOT.TLatex()

    for i in range(n):
        h1.Fill(methods[i], means_histo[i])
        h1.SetBinError(i+1, sigmas_histo[i])

        h2.Fill(methods[i], means_fit[i])
        h2.SetBinError(i+1, sigmas_fit[i])


    h1.SetStats(0)
    h1.LabelsOption("h", "X")
    h1.Draw("E1")
    h1.SetLineWidth(3)
    h1.SetMarkerStyle(20)
    h1.SetMaximum(0.025)
    h1.SetMinimum(-0.025)

    h2.SetStats(0)
    h2.Draw("E1 same")
    h2.SetLineWidth(3)
    h2.SetLineColor(2)
    h2.SetMarkerColor(2)
    h2.SetMarkerStyle(20)

    canvas.BuildLegend()
    # for i in range(n):
    #     text = ROOT.TText((i+1)/n, .9, n_pfo_survive[i])
    #     text.DrawClone()

    canvas.SetGridy()
    canvas.Update()

    input("wait")


def cylinders():
    ROOT.gStyle.SetErrorX(0.00000001);
    methods = ["{} mm".format(i) for i in range(1, 9)]
    methods.append(".....")

    methods.append("All hits")
    means = [-0.001786, -0.001601, -0.001371, -0.0008719, -0.0002143,
    0.0004927, 0.001182, 0.00181, -1, 0.003251]
    sigmas = [0.005582, 0.005656, 0.005856, 0.006463, 0.006984,
    0.007439, 0.007846, 0.008201, 0., 0.013]

    n = len(means)

    canvas = ROOT.TCanvas()
    h1 = ROOT.TH1F("h1", "Histo values; fit methods; #Delta TOF, [ns]", n, 0, n)

    for i in range(n):
        h1.Fill(methods[i], means[i])
        h1.SetBinError(i+1, sigmas[i])


    h1.SetStats(0)
    h1.LabelsOption("h", "X")
    h1.Draw("E1")
    h1.SetLineWidth(3)
    h1.SetMarkerStyle(20)
    h1.SetMaximum(0.025)
    h1.SetMinimum(-0.025)

    canvas.SetGridy()
    canvas.Update()

    input("wait")

def best():
    ROOT.gStyle.SetErrorX(0.00000001);
    methods = ["#splitline{Fit closest to line}{in 10 layers}", "#splitline{Fit 5 mm cyl}{in 30 layers}", "closest hit", "fastest hit"]
    means = [-0.002394, -0.0002143, 0.001272, 0.0003893]
    sigmas = [0.009328, 0.006984, 0.004832, 0.003205]

    n = len(means)

    canvas = ROOT.TCanvas()
    h1 = ROOT.TH1F("h1", "Histo values; fit methods; #Delta TOF, [ns]", n, 0, n)

    for i in range(n):
        h1.Fill(methods[i], means[i])
        h1.SetBinError(i+1, sigmas[i])


    h1.SetStats(0)
    h1.LabelsOption("h", "X")
    h1.Draw("E1")
    h1.SetLineWidth(3)
    h1.SetMarkerStyle(20)
    h1.SetMaximum(0.025)
    h1.SetMinimum(-0.025)

    canvas.SetGridy()
    canvas.Update()

    input("wait")

best()
