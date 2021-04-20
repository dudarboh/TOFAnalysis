import ROOT
import sys
ROOT.gInterpreter.Declare('#include "extract.hpp"')

ch = ROOT.TChain("TOFAnalysis")

# # 2f_Z_hadronic data
ch.Add("/nfs/dust/ilc/user/dudarboh/final_files/2f_Z_hadronic/result*.root")

def compare_momenta():
    canvas = ROOT.TCanvas()

    df = ROOT.RDataFrame(ch)#.Filter('if(rdfentry_ % 1000000 == 0) {cout<<"Event: "<<rdfentry_<< endl;} return true;')
    df = df.Filter("nECALHits > 0 && abs(xyzCluster.Z()) < 2200.")

    ROOT.gInterpreter.Declare('''
    double diff(const XYZVector& a, const XYZVector& b){
        return a.R() - b.R();
    }

    ''')

    h0 = df.Filter("abs(PDG) != 2212 && abs(PDG) != 321 && abs(PDG) != 211").Define("dMom", "diff(pTrackAtIP, pTrackAtCalo)*1000").Histo1D(("h0", "Backgound; #p_{IP} - p_{calo}, [MeV]; N particles", 1000, -50., 50.), "dMom")
    h1 = df.Filter("abs(PDG) == 211").Define("dMom", "diff(pTrackAtIP, pTrackAtCalo)*1000").Histo1D(("h1", "Pions; #p_{IP} - p_{calo}, [MeV]; N particles", 1000, -50., 50.), "dMom")
    h2 = df.Filter("abs(PDG) == 321").Define("dMom", "diff(pTrackAtIP, pTrackAtCalo)*1000").Histo1D(("h2", "Kaons; #p_{IP} - p_{calo}, [MeV]; N particles", 1000, -50., 50.), "dMom")
    h3 = df.Filter("abs(PDG) == 2212").Define("dMom", "diff(pTrackAtIP, pTrackAtCalo)*1000").Histo1D(("h3", "Protons; #p_{IP} - p_{calo}, [MeV]; N particles", 1000, -50., 50.), "dMom")

    h0.Draw()
    h1.Draw("same")
    h1.SetLineColor(4)
    h2.Draw("same")
    h2.SetLineColor(2)
    h3.Draw("same")
    h3.SetLineColor(3)

    canvas.Update()
    input("wait")

def compare_length():
    canvas = ROOT.TCanvas()

    df = ROOT.RDataFrame(ch)#.Filter('if(rdfentry_ % 1000000 == 0) {cout<<"Event: "<<rdfentry_<< endl;} return true;')
    df = df.Filter("nECALHits > 0 && abs(xyzCluster.Z()) < 2200.")

    ROOT.gInterpreter.Declare('''
    double diff(const XYZVector& a, const XYZVector& b){
        return a.R() - b.R();
    }

    ''')

    h1 = df.Filter("abs(PDG) == 211").Define("dLen", "lengthTrackIP - lengthTrackCalo").Histo1D(("h1", "Pions; l_{IP} - l_{calo}, [mm]; N particles", 1000, -50., 50.), "dLen")
    h2 = df.Filter("abs(PDG) == 321").Define("dLen", "lengthTrackIP - lengthTrackCalo").Histo1D(("h2", "Kaons; l_{IP} - l_{calo}, [mm]; N particles", 1000, -50., 50.), "dLen")
    h3 = df.Filter("abs(PDG) == 2212").Define("dLen", "lengthTrackIP - lengthTrackCalo").Histo1D(("h3", "Protons; l_{IP} - l_{calo}, [mm]; N particles", 1000, -50., 50.), "dLen")

    h1.Draw("")
    h2.Draw("same")
    h2.SetLineColor(2)
    h3.Draw("same")
    h3.SetLineColor(3)

    canvas.Update()
    input("wait")

def check_integrated_length():
    canvas = ROOT.TCanvas()

    df = ROOT.RDataFrame(ch)#.Filter('if(rdfentry_ % 1000000 == 0) {cout<<"Event: "<<rdfentry_<< endl;} return true;')
    df = df.Filter("nECALHits > 0 && abs(xyzCluster.Z()) < 2200.")

    ROOT.gInterpreter.Declare('''
    double diff(const XYZVector& a, const XYZVector& b){
        return a.R() - b.R();
    }

    ''')

    h1 = df.Define("dLen", "lengthTrackIP - lengthTrackIntegral").Histo1D(("h1", "IP; l_{IP} - l_{int}, [mm]; N particles", 1000, -500., 500.), "dLen")
    h2 = df.Define("dLen", "lengthTrackCalo - lengthTrackIntegral").Histo1D(("h2", "Calo; l_{IP} - l_{int}, [mm]; N particles", 1000, -500., 500.), "dLen")

    h1.Draw("")
    h2.Draw("same")
    h2.SetLineColor(2)
    canvas.Update()
    input("wait")


def check_shower():
    canvas = ROOT.TCanvas()
    for idx, event in enumerate(ch):
        if event.nECALHits > 0 and abs(event.xyzCluster.Z()) < 2200.:
            print("PDG", event.PDG)
            x0 = event.xyzTrackAtCalo.X()
            print("x0", x0)
            y0 = event.xyzTrackAtCalo.Y()
            print("y0", y0)
            px = event.pTrackAtCalo.X()
            print("px", px)
            py = event.pTrackAtCalo.Y()
            print("py", py)
            gr_line = ROOT.TF1("f1", "{0} + {1}/{2}*(x - {3})".format(y0, py, px, x0), -3000., 3000.)
            gr_line.SetLineColor(2)
            gr_point = ROOT.TGraph()
            gr_point.SetPoint(0, x0, y0)
            gr_point.SetMarkerStyle(29)
            gr_point.SetMarkerColor(4)
            gr_point.SetMarkerSize(2.5)

            gr_shower = ROOT.TGraph()
            for i, hit in enumerate(event.xyzECALHit):
                gr_shower.SetPoint(i+1, hit.X(), hit.Y())

            gr_shower.Draw("AP")
            gr_shower.SetTitle("Shower vs track extrapolation;x, [mm];y, [mm]")
            gr_shower.SetMarkerStyle(20)
            gr_line.Draw("same")
            gr_point.Draw("Psame")
            gr_shower.GetXaxis().SetRangeUser(-3000., 3000.)
            gr_shower.GetYaxis().SetRangeUser(-3000., 3000.)
            canvas.Update()
            input("wait")
            gr_fit = ROOT.TGraph()
            for i, (pos, t) in enumerate( zip(event.xyzECALHit, event.tECALHit) ):
                gr_fit.SetPoint(i, (pos - event.xyzTrackAtCalo).R(), t)
            gr_fit.Fit("pol1", "Q");
            gr_fit.SetLineStyle(20)
            gr_fit.Draw("AP")
            input("wait")



# check_integrated_length()
check_shower()
