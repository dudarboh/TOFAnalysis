from ROOT import TFile, TH1F, TCanvas, TH2F, gStyle, TGraph, THStack, TChain
import time
from algorithms import *

# This all for 3D plots
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt

# c in mm/ns
c = 299.792458

# calib_kaons data
# file = TFile("/nfs/dust/ilc/user/dudarboh/final_files/calib_kaons.root", "READ")
# tree = file.ana_tree

# 2f_Z_hadronic data
tree = TChain("ana_tree")
tree.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/*.root")

# Deprecated. I dont store event info anymore
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


def plot_track(pfo_number):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.grid(False)
    colors = ["b", "g", "r", "c", "m", "y", "k"]
    for idx, pfo in enumerate(tree):
        if idx != pfo_number:
            continue
        # plot track
        x = []
        y = []
        z = []
        n_hits = pfo.nTrHits
        for i in range(n_hits):
            x.append(pfo.xTrHit[i])
            y.append(pfo.yTrHit[i])
            z.append(pfo.zTrHit[i])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        ax.scatter3D(z, x, y, color=colors[pfo%len(colors)])

        x = []
        y = []
        z = []
        # Plot shower
        n_hits = pfo.nCalHits
        for i in range(n_hits):
            x.append(pfo.xCalHit[i])
            y.append(pfo.yCalHit[i])
            z.append(pfo.zCalHit[i])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        ax.scatter3D(z, x, y, color=colors[pfo%len(colors)], s=30, marker="s")


        x_calo = pfo.xRefCalo
        y_calo = pfo.yRefCalo
        z_calo = pfo.zRefCalo
        ax.scatter3D(z_calo, x_calo, y_calo, color="red", s=130, marker="*")

        x_last = pfo.xRefLast
        y_last = pfo.yRefLast
        z_last = pfo.zRefLast
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
    beta = "length/(Min$(tHitCluster)*{})".format(c)
    tree.Draw("{}:p>>h1".format(beta), "", "COLZ")
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


def tof_analysis():
    canvas = TCanvas()
    h1 = TH1F("h1", "TOF fit", 400, 0., 1.)
    h1.SetLineColor(1)
    h2 = TH1F("h2", "TOF avg", 400, 0., 1.)
    h2.SetLineColor(2)
    h3 = TH1F("h3", "TOF closest", 400, 0., 1.)
    h3.SetLineColor(3)
    hs = THStack()
    hs.Add(h1)
    hs.Add(h2)
    hs.Add(h3)

    n_events_total = tree.GetEntries()
    for idx, event in enumerate(tree):
        if idx % 200 == 0:
            print(idx, "event from", n_events_total)
        # if idx == 10000:
        #     break

        n_pfo = event.nGoodPFOs
        for pfo in range(n_pfo):
            p = event.p[pfo]
            l_trk = event.length[pfo]
            hits = []
            for (t, r, d, layer) in zip(event.tHitCluster[pfo], event.dToRefPointHitCluster[pfo], event.dToLineHitCluster[pfo], event.layerHitCluster[pfo]):
                if layer > 9:
                    continue
                hits.append({"t": t, "r" : r, "d" : d, "layer" : layer})
            if hits == []:
                continue

            mass_closest = algo_closest(p, l_trk, hits)
            mass_avg = algo_avg(p, l_trk, hits)
            mass_fit = algo_fit(p, l_trk, hits)
            h1.Fill(mass_fit)
            h2.Fill(mass_avg)
            h3.Fill(mass_closest)
    hs.Draw("nostack")
    canvas.BuildLegend()
    canvas.Update()
    input("wait")





# plot_event(2)
# plot_track(0, 1)
# beta_p()
# plot_mass()
tof_analysis()
# plot_histo()
