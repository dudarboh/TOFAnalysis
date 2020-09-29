from ROOT import TFile, TH1F, TCanvas, TH2F, gStyle, TGraph, THStack
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
import time

# c in mm/ns
c = 299.792458

file = TFile("/nfs/dust/ilc/user/dudarboh/final_files/calib_kaons.root", "READ")
tree = file.ana_tree

def plot_event(ev):

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.grid(False)
    colors = ["b", "g", "r", "c", "m", "y", "k"]
    for idx, event in enumerate(tree):
        if idx != ev:
            continue
        n_pfo = event.nGoodPFOs
        print(n_pfo)
        for pfo in range(n_pfo):
            # Plot tracker hits
            x = []
            y = []
            z = []
            n_hits = event.nHitsTrack[pfo]
            for i in range(n_hits):
                x.append(event.xHit[pfo][i])
                y.append(event.yHit[pfo][i])
                z.append(event.zHit[pfo][i])

            x = np.array(x)
            y = np.array(y)
            z = np.array(z)

            ax.scatter3D(z, x, y, color=colors[pfo%len(colors)])
            # Plot cluster
            x = []
            y = []
            z = []
            n_hits = event.nHitsCluster[pfo]
            for i in range(n_hits):
                x.append(event.xClusterHit[pfo][i])
                y.append(event.yClusterHit[pfo][i])
                z.append(event.zClusterHit[pfo][i])

            x = np.array(x)
            y = np.array(y)
            z = np.array(z)

            ax.scatter3D(z, x, y, color=colors[pfo%len(colors)], s=30, marker="s")


    ax.set_xlabel("Z axis, [mm]")
    ax.set_ylabel("X axis, [mm]")
    ax.set_zlabel("Y axis, [mm]")

    # ECAL
    r = 1805.
    h = 2350.

    x=np.linspace(-r, r, 100)
    z=np.linspace(-h, h, 100)
    Xc, Zc=np.meshgrid(x, z)
    Yc = np.sqrt(r**2-Xc**2)
    # Draw parameters
    rstride = 20
    cstride = 10
    ax.plot_surface(Zc, Xc, Yc, alpha=0.15, rstride=rstride, cstride=cstride, color="cyan")
    ax.plot_surface(Zc, Xc, -Yc, alpha=0.15, rstride=rstride, cstride=cstride, color="cyan")
    # Add endcaps
    p1 = Circle((0, 0), r, alpha=0.15)
    p2 = Circle((0, 0), r, alpha=0.15)
    ax.add_patch(p1)
    ax.add_patch(p2)
    art3d.pathpatch_2d_to_3d(p1, z=h, zdir="x")
    art3d.pathpatch_2d_to_3d(p2, z=-h, zdir="x")

    plt.show()

    input("wait")



def plot_track(ev, track):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.grid(False)
    colors = ["b", "g", "r", "c", "m", "y", "k"]
    for idx, event in enumerate(tree):
        if idx != ev:
            continue
        # plot track
        n_pfo = event.nGoodPFOs
        print(n_pfo)
        pfo = track
        x = []
        y = []
        z = []
        n_hits = event.nHitsTrack[pfo]
        for i in range(n_hits):
            x.append(event.xHit[pfo][i])
            y.append(event.yHit[pfo][i])
            z.append(event.zHit[pfo][i])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        ax.scatter3D(z, x, y, color=colors[pfo%len(colors)])

        x = []
        y = []
        z = []
        # Plot shower
        n_hits = event.nHitsCluster[pfo]
        for i in range(n_hits):
            x.append(event.xClusterHit[pfo][i])
            y.append(event.yClusterHit[pfo][i])
            z.append(event.zClusterHit[pfo][i])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        ax.scatter3D(z, x, y, color=colors[pfo%len(colors)], s=30, marker="s")


        x_calo = event.xRefCalo[pfo]
        y_calo = event.yRefCalo[pfo]
        z_calo = event.zRefCalo[pfo]
        ax.scatter3D(z_calo, x_calo, y_calo, color="red", s=130, marker="*")

        x_last = event.xRefLast[pfo]
        y_last = event.yRefLast[pfo]
        z_last = event.zRefLast[pfo]
        ax.scatter3D(z_last, x_lakst, y_last, color="yellow", s=130, marker="*")

    ax.set_xlabel("Z axis, [mm]")
    ax.set_ylabel("X axis, [mm]")
    ax.set_zlabel("Y axis, [mm]")

    # ECAL
    r = 1805.
    h = 2350.

    x=np.linspace(-r, r, 100)
    z=np.linspace(-h, h, 100)
    Xc, Zc=np.meshgrid(x, z)
    Yc = np.sqrt(r**2-Xc**2)
    # Draw parameters
    rstride = 20
    cstride = 10
    ax.plot_surface(Zc, Xc, Yc, alpha=0.15, rstride=rstride, cstride=cstride, color="cyan")
    ax.plot_surface(Zc, Xc, -Yc, alpha=0.15, rstride=rstride, cstride=cstride, color="cyan")
    # Add endcaps
    p1 = Circle((0, 0), r, alpha=0.15)
    p2 = Circle((0, 0), r, alpha=0.15)
    ax.add_patch(p1)
    ax.add_patch(p2)
    art3d.pathpatch_2d_to_3d(p1, z=h, zdir="x")
    art3d.pathpatch_2d_to_3d(p2, z=-h, zdir="x")

    plt.show()

    input("wait")



def beta_p():
    # canvas = TCanvas()
    gStyle.SetPalette(1)
    nBinsX = 100
    nBinsY = 100
    minBinX = 0.
    maxBinX = 1. # log10(min/max momentum / GeV)
    minBinY = .7
    maxBinY = 1.1; # for TOF beta = v/c values from 0 to 1

    histbinsX = []
    histbinsY = []
    for i in range(nBinsX + 1):
        histbinsX.append(10**(minBinX + 1.*(maxBinX-minBinX)*i/nBinsX))
    for i in range(nBinsY + 1):
        histbinsY.append(minBinY + 1.*(maxBinY-minBinY)*i/nBinsY)

    h1 = TH2F("h1", "title", nBinsX, np.array(histbinsX), nBinsY, np.array(histbinsY))
    h1.SetTitle("Beta vs p plot; p, [GeV]; #beta")
    length = "abs((phiCalo-phi)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo)"
    beta = "{}/(Min$(tHitCluster)*{})".format(length, c)
    tree.Draw("{}:p>>h1".format(beta), "", "COLZ")
    input("wait")


def plot_mass():
    canvas = TCanvas()
    h1 = TH1F("h1", "title", 1000, 0.3, 1.3)
    h1.SetTitle("title; mass, [GeV]; N tracks")

    length = "(abs((phiCalo-phi)/omega)*sqrt(1. + tanL*tanL))"

    beta = "({}/({}*{}))".format(length, tof, c)
    mass = "p*sqrt(1.-{0}*{0})/{0}".format(beta)

    tree.Draw("{}>>h1".format(mass), "nHitsTrack > 200 && p > 1 && p < 10")
    input("wait")



def tof_analysis():
    canvas = TCanvas()
    h1 = TH1F("h1", "TOF fit", 100, 5., 20.)
    h1.SetLineColor(1)
    h2 = TH1F("h2", "TOF avg", 100, 5., 20.)
    h2.SetLineColor(2)
    h3 = TH1F("h3", "TOF closest", 100, 5., 20.)
    h3.SetLineColor(3)

    for idx, event in enumerate(tree):
        if idx % 1000 == 0:
            print(idx, "event")
        # if idx == 1000:
        #     break

        n_pfo = event.nGoodPFOs
        for pfo in range(n_pfo):
            # Make hit list
            hits = []
            for h in range(event.nHitsCluster[pfo]):
                t = event.tHitCluster[pfo][h]
                r = event.dToRefPointHitCluster[pfo][h]
                d = event.dToLineHitCluster[pfo][h]
                layer = event.layerHitCluster[pfo][h]
                hits.append((t, r, d, layer))
            if hits == []:
                continue

            length = event.length[pfo]

            # Closest
            hits_closest = sorted(hits, key=lambda x: x[1])
            tof_closest = hits_closest[0][0] - hits_closest[0][1]/c
            if tof_closest != 0:
                beta_closest = length/(tof_closest*c)
                if 0. < beta_closest < 1.:
                    # h3.Fill(event.p[pfo] * np.sqrt(1. - beta_closest*beta_closest)/beta_closest)
                    h3.Fill(tof_closest)

            # Franks Average
            hits_by_layers = [[] for i in range(10)]
            for t, r, d, layer in hits:
                if layer > 9:
                    continue
                hits_by_layers[int(layer)].append((t, r, d))
            sum = 0.
            arr_length = 0.
            for l in hits_by_layers:
                if l == []:
                    continue
                l.sort(key=lambda x: x[2])
                sum += l[0][0] - l[0][1]/c
                arr_length += 1.

            tof_avg = sum/arr_length if arr_length !=0. else 0.
            if tof_avg != 0:
                beta_avg = length/(tof_avg*c)
                if 0. < beta_avg < 1.:
                    # h2.Fill(event.p[pfo] * np.sqrt(1. - beta_avg*beta_avg)/beta_avg)
                    h2.Fill(tof_avg)

            #Fit
            gr = TGraph()
            n = 0
            for t, r, d, layer in hits:
                # if d < 10. and r < 100:
                gr.SetPoint(n, r, t)
                n += 1
            if n == 0:
                continue

            gr.Fit("pol1", "Q")
            fit = gr.GetFunction("pol1")
            tof_fit = fit.GetParameter(0)
            if 1.95 < event.p[pfo] < 2.05:
                gr.Draw("AP")
                print(event.p[pfo])
                gr.SetMarkerStyle(20)
                gr.SetTitle("Fit of a TOF vs distance to impact point; d, [mm];TOF, [ns]")
                canvas.Update()
            # time.sleep(3)
                raw_input("wait")
            if fit.GetChisquare() > 1.:
                continue
            if tof_fit != 0:
                beta_fit = length/(tof_fit*c)
                if 0. < beta_fit < 1.:
                    # h1.Fill(event.p[pfo] * np.sqrt(1. - beta_fit*beta_fit)/beta_fit)
                    h1.Fill(tof_fit)

    hs = THStack()
    hs.Add(h1)
    hs.Add(h2)
    hs.Add(h3)
    hs.Draw("nostack")
    canvas.BuildLegend()
    canvas.Update()
    input("wait")


def plot_histo():
    canvas = TCanvas()
    h1 = TH1F("h1", "Distance between last tracker hit and impact point in ECAL;d, [mm];N_{PFOs}", 600, 0, 400)
    # var = "(zRefCalo-zRefLast)"

    # length = "abs((phiCalo-phi)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo)"
    # var = "{}/(Sum$(tHitCluster)/Length$(tHitCluster)*{})".format(length, c)

    var = "sqrt((xRefCalo-xRefLast)*(xRefCalo-xRefLast) + (yRefCalo-yRefLast)*(yRefCalo-yRefLast) + (zRefCalo-zRefLast)*(zRefCalo-zRefLast))"

    # var = "abs((phiCalo-phi)*(1./omegaCalo))*sqrt(1. + tanLCalo*tanLCalo)"

    # var = "dToLineHitCluster"

    tree.Draw("{}>>h1".format(var))
    input("wait")



# plot_event(2)
# plot_track(0, 1)
# beta_p()
# plot_mass()
tof_analysis()
# plot_histo()
