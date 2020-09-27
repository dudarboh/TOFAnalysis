from ROOT import TFile, TH1F, TCanvas, TH2F, gStyle, TGraph
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
    beta = "({}/(MinIf$(tHitCluster, tHitCluster>0)*{}))".format(length, c)
    mass = "p*sqrt(1.-{0}*{0})/{0}".format(beta)

    tree.Draw("{}>>h1".format(mass), "nHitsTrack > 200 && p > 1 && p < 10")
    input("wait")



def plot_time():
    canvas = TCanvas()
    for idx, event in enumerate(tree):
        if idx % 1000 == 0:
            print(idx, "event")
        n_pfo = event.nHitsCluster.size()
        for pfo in range(n_pfo):
            gr = TGraph()
            gr.SetMarkerStyle(20)
            n = 0
            times = []
            x0 = event.xRefCalo[pfo]
            y0 = event.yRefCalo[pfo]
            z0 = event.zRefCalo[pfo]
            r0 = np.sqrt(x0*x0+y0*y0+z0*z0)
            for h in range(event.nHitsCluster[pfo]):
                t = event.tHitCluster[pfo][h]
                x = event.xHitCluster[pfo][h]
                y = event.yHitCluster[pfo][h]
                z = event.zHitCluster[pfo][h]
                r = np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
                if event.dToLineHitCluster[pfo][h] > 20.:
                    continue
                times.append((t, r))

            times.sort(key=lambda x: x[1])
            for t, d in times:
                gr.SetPoint(n, d, t)
                n += 1

            gr.Draw("AP")
            canvas.Update()
            time.sleep(1.)


def plot_histo():
    canvas = TCanvas()
    h1 = TH1F("h1", "title", 200, 0, 150)
    h1.SetTitle("title; x; y")
    # var = "(zRefCalo-zRefLast)"

    # length = "abs((phiCalo-phi)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo)"
    # var = "{}/(Sum$(tHitCluster)/Length$(tHitCluster)*{})".format(length, c)

    # var = "sqrt((xRefCalo-xRefLast)*(xRefCalo-xRefLast) + (yRefCalo-yRefLast)*(yRefCalo-yRefLast) + (zRefCalo-zRefLast)*(zRefCalo-zRefLast))"

    # var = "abs((phiCalo-phi)*(1./omegaCalo))*sqrt(1. + tanLCalo*tanLCalo)"

    var = "dToLineHitCluster"

    tree.Draw("{}>>h1".format(var))
    input("wait")



# plot_event(2)
# plot_track(0, 1)
# plot_histo()
# beta_p()
# plot_mass()
plot_time()
