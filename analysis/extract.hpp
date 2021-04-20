#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include <random>
using namespace ROOT::Math;
using namespace ROOT::VecOps;
using std::cout, std::endl, std::stringstream, std::vector, std::string, std::runtime_error;

// in mm/ns
#define SPEED_OF_LIGHT 299.792458
// const double rInner = 1804.8;

// smearing generator and distributions
std::default_random_engine generator;
std::normal_distribution <double> gaus10(0., 0.01);
std::normal_distribution <double> gaus30(0., 0.03);
std::normal_distribution <double> gaus50(0., 0.05);
std::normal_distribution <double> gaus100(0., 0.1);
std::normal_distribution <double> gaus200(0., 0.2);
std::normal_distribution <double> gaus300(0., 0.3);

TCanvas* canvas = new TCanvas("canvas", "title", 1920, 1600);
int pic = 0;
XYZVector getECALPlane(const double& phiCluster){
    int nSides = 8;
    double side = 2 * M_PI / nSides;

    double phiPlane = -M_PI + side;
    for (int i = 0; i < nSides-1; ++i){
        if ( phiPlane - side/2 <= phiCluster && phiCluster < phiPlane + side/2 )
            return XYZVector(Polar3DVector(1804.8, M_PI/2, phiPlane));
        else phiPlane += side;
    }
    // if ( M_PI - side/2 <= phiCluster || phiCluster < -M_PI + side/2 )
    return XYZVector(Polar3DVector(1804.8, M_PI/2, M_PI));
}


XYZVector intersection(const XYZVector& vtx, const XYZVector& p, const XYZVector& plane){
    XYZVector n = plane.Unit();
    double prod = p.Dot(n);
    if (prod == 0.) return XYZVector(); // Plane and line are parallel
    double lambda = (plane - vtx).Dot(n)/prod;
    return vtx + lambda * p;
}


RVec <double> dToImpact(const RVec<XYZVector>& hits, const XYZVector& rImpact){
    auto getDistance = [&](const XYZVector& hit) { return (hit - rImpact).R(); };
    return Map(hits, getDistance);
}


RVec <double> dToLine(const RVec <XYZVector>& hits, const XYZVector& vtx, const XYZVector& p){
    RVec <double> distance;
    for (const auto& hit:hits){
        double d = (hit - vtx).Cross(p.Unit()).R();
        distance.push_back(d);
    }
    return distance;
}

RVec <bool> selectHits(const RVec<double>& dToLine, const RVec<int>& layer_hit, bool only_closest=true, int n_layers=10, double cyl_cut=9999.){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    if(!only_closest){
        for (int i=0; i < nHits; ++i)
            if (dToLine[i] < cyl_cut && layer_hit[i] < n_layers)
                selected[i] = true;
        return selected;
    }

    for (int layer=0; layer<n_layers; ++layer){
        map<int, double> layer_hits;
        for (int i=0; i < nHits; ++i)
            if(dToLine[i] < cyl_cut && layer_hit[i] == layer)
                layer_hits[i] = dToLine[i];

        int min_idx =(*min_element(layer_hits.begin(), layer_hits.end(),
                        [](const auto& l, const auto& r) { return l.second < r.second; }) ).first ;
        selected[min_idx] = true;
    }
    return selected;
}


/////////////////////////TOFS//////////////////////////
/////////////////////////TOFS//////////////////////////
/////////////////////////TOFS//////////////////////////
RVec <double> tofHit(const RVec <double>& tHits, int resolution=0){
    //resolution option in picoseconds
    RVec <double> tofs;
    int nHits = tHits.size();
    double smearing = 0.;
    for (int i = 0; i < nHits; ++i){
        if (resolution == 10) smearing = gaus10(generator);
        else if (resolution == 30) smearing = gaus30(generator);
        else if (resolution == 50) smearing = gaus50(generator);
        else if (resolution == 100) smearing = gaus100(generator);
        else if (resolution == 200) smearing = gaus200(generator);
        else if (resolution == 300) smearing = gaus300(generator);
        else smearing = 0.;
        tofs.push_back( tHits[i] + smearing );
    }
    return tofs;
}


double tofClosest(const RVec<double>& tofHit, const RVec<double>& dToImpact){
    int min_idx = min_element(dToImpact.begin(), dToImpact.end()) - dToImpact.begin();
    return tofHit[min_idx] - dToImpact[min_idx]/SPEED_OF_LIGHT;
}


double tofFastest(const RVec<double>& tofHit, const RVec<double>& dToImpact){
    int min_idx = min_element(tofHit.begin(), tofHit.end()) - tofHit.begin();
    return tofHit[min_idx] - dToImpact[min_idx]/SPEED_OF_LIGHT;
}


double fitFunc(const RVec <double>& tofHit, const RVec<double>& dToImpact, const int opt=0, int limits=0){
    int nHits = tofHit.size();
    if(nHits <= 1) return 0;
    vector<double> x;
    vector<double> x_err;
    vector<double> y;
    vector<double> y_err;

    for(int i=0; i < nHits; ++i){
        x.push_back(dToImpact[i]);
        y.push_back(tofHit[i]);
        x_err.push_back(0.);
        y_err.push_back(0.1);
    }

    TGraphErrors gr(x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
    TF1* fit = new TF1("fit", "pol1");
    if (limits){
        fit->SetParameters(5., 1./(0.95*SPEED_OF_LIGHT));
        fit->SetParLimits(0, 5., 20.);
        fit->SetParLimits(1, 1./(1.01*SPEED_OF_LIGHT), 1./(0.5*SPEED_OF_LIGHT));
    }

    if (limits){
        gr.Fit(fit, "QB");
    }
    else{
        gr.Fit(fit, "Q");
    }

    fit = gr.GetFunction("fit");

    double result = -1;
    switch (opt) {
        case 0: result = fit->GetParameter(0); break;
        case 1: result = 1./fit->GetParameter(1); break;
        case 2: result = fit->GetChisquare(); break;
        case 3: result = fit->GetNDF(); break;
    }
    return result;
}


double tofAvg(const RVec <double>& tofHit, const RVec<double>& dToImpact){
    int nHits = tofHit.size();
    double tofSum = 0.;
    for(int i=0; i < nHits; ++i) tofSum += tofHit[i] - dToImpact[i]/SPEED_OF_LIGHT;
    return tofSum/nHits;
}

int fit_analysis(const RVec <double>& tofHit, const RVec<double>& dToImpact, const RVec <double>& tofHitNoSmearing, int pdg, XYZVector p, int nClusterHits){
    cout<< "Debug1"<<endl;
    int nHits = tofHit.size();
    cout<< "N hits"<<nHits<<endl;
    if(nHits <= 1) return 0;
    vector<double> x;
    vector<double> x_err;
    vector<double> y;
    vector<double> y_err;
    vector<double> y_no_smear;


    for(int i=0; i < nHits; ++i){
        x.push_back(dToImpact[i]);
        y.push_back(tofHit[i]);
        x_err.push_back(0.);
        y_err.push_back(0.1);
        y_no_smear.push_back(tofHitNoSmearing[i]);
    }
    cout<< "Debug2"<<endl;

    TGraphErrors gr(x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
    TGraphErrors gr_no_smear(x.size(), &x[0], &y_no_smear[0], &x_err[0], &y_err[0]);
    TF1* fit = new TF1("fit", "pol1");
    TF1* fit_no_smear = new TF1("fit_no_smear", "pol1");
    TF1* fit_limits = new TF1("fit_limits", "pol1");
    fit_limits->SetParameters(5., 1./(0.95*SPEED_OF_LIGHT));
    fit_limits->SetParLimits(0, 5., 20.);
    fit_limits->SetParLimits(1, 1./(1.01*SPEED_OF_LIGHT), 1./(0.5*SPEED_OF_LIGHT));
    TF1* f_2 = new TF1("f_2", "pol1", 0, 500);
    TF1* f_1 = new TF1("f_1", "pol1", 0, 500);
    TF1* f_22 = new TF1("f_22", "pol1", 0, 500);
    TF1* f_11 = new TF1("f_11", "pol1", 0, 500);

    cout<< "Debug3"<<endl;
    gr.Fit(fit);
    fit = gr.GetFunction("fit");
    gr.Fit(fit_limits, "B+");
    fit_limits = gr.GetFunction("fit_limits");
    gr_no_smear.Fit(fit_no_smear, "B+");
    fit_no_smear = gr_no_smear.GetFunction("fit_no_smear");
    cout<< "Debug4"<<endl;

    gr.SetTitle("ECAL hits, 100 ps");
    gr.SetMarkerStyle(20);
    gr.SetMarkerColor(1);
    gr.SetLineColor(1);
    gr.SetLineWidth(1);
    fit->SetLineColor(1);
    fit->SetLineWidth(3);
    fit->SetTitle("Fit 100 ps");
    fit->SetNpx(1000);

    gr_no_smear.SetTitle("ECAL hits, 0 ps");
    gr_no_smear.SetMarkerStyle(20);
    gr_no_smear.SetMarkerColor(4);
    gr_no_smear.SetLineColor(4);
    gr_no_smear.SetLineWidth(1);
    fit_no_smear->SetLineColor(4);
    fit_no_smear->SetLineWidth(3);
    fit_no_smear->SetTitle("Fit 0 ps");
    fit_no_smear->SetNpx(1000);

    fit_limits->SetLineColor(2);
    fit_limits->SetLineWidth(3);
    fit_limits->SetTitle("Fit 100 ps + param limits");
    fit_limits->SetNpx(1000);
    cout<< "Debug5"<<endl;


    TMultiGraph mg = TMultiGraph();
    mg.Add(&gr, "P");
    mg.Add(&gr_no_smear, "P");
    mg.Draw("A");
    mg.SetTitle("Trying to improve fit; d to impact, [mm]; Time, [ns]");
    mg.GetXaxis()->SetLimits(0., mg.GetXaxis()->GetBinUpEdge(mg.GetXaxis()->GetLast()));

    cout<< "Debug6"<<endl;
    fit->Draw("same");
    fit_no_smear->Draw("same");
    fit_limits->Draw("same");
    cout<< "Debug7"<<endl;
    canvas->BuildLegend(0.1, 0.7, 0.3, 0.9);

    TLatex text;
    text.SetTextSize(0.02);
    text.DrawLatexNDC(0.35, 0.88, (const char*)Form("PDG (%d)" , pdg));
    text.DrawLatexNDC(0.35, 0.86, (const char*)Form("nECALHits (%d)", nClusterHits));
    text.DrawLatexNDC(0.35, 0.84, (const char*)Form("pt (%f.2 GeV)", p.Rho()));
    text.DrawLatexNDC(0.35, 0.82, (const char*)Form("p (%f.2 GeV)", p.R()));
    text.DrawLatexNDC(0.35, 0.8, Form("#theta (%f.2) degree", p.Theta()/(2.*3.1415692)*180.));
    text.DrawLatexNDC(0.12, 0.68, Form("#color[4]{TOF 0ps (%f.2), [ns]}", fit_no_smear->GetParameter(0)));
    text.DrawLatexNDC(0.12, 0.66, Form("TOF 100ps (%f.2), [ns]", fit->GetParameter(0)));
    text.DrawLatexNDC(0.12, 0.64, Form("#color[2]{TOF 100ps par limits (%f.2), [ns]}", fit_limits->GetParameter(0)));
    text.DrawLatexNDC(0.12, 0.62, Form("#color[4]{#beta 0ps (%f.2), [mm/ns]}", 1./fit_no_smear->GetParameter(1)));
    text.DrawLatexNDC(0.12, 0.6, Form("#beta 100ps (%f.2), [mm/ns]", 1./fit->GetParameter(1)));
    text.DrawLatexNDC(0.12, 0.58, Form("#color[2]{#beta 100ps par limits (%f.2), [mm/ns]}", 1./fit_limits->GetParameter(1)));

    canvas->Modified();
    canvas->Update();
    // cin.get();
    canvas->Print(Form("./fit_analysis_pics/%d.png", pic++));
    return true;
}



/////////////////////////EXTRA//////////////////////////
/////////////////////////EXTRA//////////////////////////
/////////////////////////EXTRA//////////////////////////
RVec <double> residual(const RVec <double>& tHits, const RVec<double>& dToImpact,const double& par0,const double& par1){
    RVec <double> residual;
    int nHits = tHits.size();
    for (int i = 0; i < nHits; ++i)
        residual.push_back( tHits[i] - (par0 + par1*dToImpact[i]) );
    return residual;
}


int getAbsMax(const RVec <double>& res){
	return max_element(res.begin(), res.end(), [](const double& a, const double& b) {return abs(a) < abs(b);}) - res.begin();
}
