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


double fitFunc(const RVec <double>& tofHit, const RVec<double>& dToImpact, const int opt=0){
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
        // y_err.push_back(0.0121634); // for ALL hits
        // y_err.push_back(0.00325199); // for 4 mm cylinder hits
        // y_err.push_back(4.72659e-03); // for 10 mm cylinder hits
    }
    TGraphErrors gr(x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
    gr.Fit("pol1", "Q");

    double result = -1;
    switch (opt) {
        case 0: result = gr.GetFunction("pol1")->GetParameter(0); break;
        case 1: result = gr.GetFunction("pol1")->GetParameter(1); break;
        case 2: result = gr.GetFunction("pol1")->GetChisquare(); break;
        case 3: result = gr.GetFunction("pol1")->GetNDF(); break;
    }
    return result;
}


double tofAvg(const RVec <double>& tofHit, const RVec<double>& dToImpact){
    int nHits = tofHit.size();
    double tofSum = 0.;
    for(int i=0; i < nHits; ++i) tofSum += tofHit[i] - dToImpact[i]/SPEED_OF_LIGHT;
    return tofSum/nHits;
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
