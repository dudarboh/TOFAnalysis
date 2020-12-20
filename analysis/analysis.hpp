#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <random>
using namespace ROOT::Math;
using namespace ROOT::VecOps;
using std::cout, std::endl, std::stringstream, std::vector, std::string, std::runtime_error;

// in mm/ns
#define SPEED_OF_LIGHT 299.792458;
// const double rInner = 1804.8;

// 50 ps smearing generator
std::default_random_engine generator;
std::normal_distribution <double> distribution(0.0, 0.05);


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

RVec <double> dToImpact(const RVec<XYZTVector>& hits, const XYZVector& rImpact){
    auto getDistance = [&](const XYZTVector& hit) { return (hit.Vect() - rImpact).R(); };
    return Map(hits, getDistance);
}


double tofClosest(const RVec<XYZTVector>& hits, const RVec<double>& dToImpact){
    int min_idx = min_element(dToImpact.begin(), dToImpact.end()) - dToImpact.begin();
    return hits[min_idx].T() - dToImpact[min_idx]/SPEED_OF_LIGHT;
}


double tofFastest(const RVec<XYZTVector>& hits, const RVec<double>& dToImpact){
    int min_idx = min_element(hits.begin(), hits.end(),
                                    [](const XYZTVector& hit1, const XYZTVector& hit2)
                                    {return hit1.T() < hit2.T();} ) - hits.begin();
    return hits[min_idx].T() - dToImpact[min_idx]/SPEED_OF_LIGHT;
}


RVec <double> dToLine(const RVec <XYZTVector>& hits, const XYZVector& vtx, const XYZVector& p){
    RVec <double> distance;
    for (auto&& hit:hits){
        double d = (hit.Vect() - vtx).Cross(p.Unit()).R();
        distance.push_back(d);
    }
    return distance;
}

RVec <bool> selectHits(const RVec<double>& dToLine, const RVec<int>& layer_hit, int n_layers=10){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    for (int layer=0; layer<n_layers; ++layer){
        map<int, double> selected_hits;
        for (int i=0; i < nHits; ++i)
            if(layer_hit[i] == layer) selected_hits[i] = dToLine[i];

        int min_idx =(*min_element(selected_hits.begin(), selected_hits.end(),
                        [](const auto& l, const auto& r) { return l.second < r.second; }) ).first ;
        selected[min_idx] = true;
    }
    return selected;
}


double tofAvg(const RVec <XYZTVector>& hits, const RVec<double>& dToImpact){
    double tof = 0;
    int nHits = hits.size();
    for (int i=0; i < nHits; ++i) tof += hits[i].T() - dToImpact[i]/SPEED_OF_LIGHT;
    return tof / nHits;
}



double fitFunc(const RVec <XYZTVector>& hits, const RVec<double>& dToImpact, int opt){
    int nHits = hits.size();
    double x[nHits], x_err[nHits], y[nHits], y_err[nHits];
    for(int i=0; i < nHits; ++i){
        x[i] = dToImpact[i];
        y[i] = hits[i].T();
        x_err[i] = 0.;
        y_err[i] = 0.05;
    }
    TGraphErrors gr(nHits, x, y, x_err, y_err);
    gr.Fit("pol1", "Q");
    TF1 fit = *gr.GetFunction("pol1");
    switch (opt) {
        case 0: return fit.GetParameter(0);
        case 1: return fit.GetParameter(1);
        case 2: return fit.GetChisquare();
        case 3: return fit.GetNDF();
    }
    return 0;
}

RVec <double> residual(const RVec <XYZTVector>& hits, const RVec<double>& dToImpact,const double& par0,const double& par1){
    RVec <double> residual;
    int nHits = hits.size();
    for (int i = 0; i < nHits; ++i)
        residual.push_back( hits[i].T() - (par0 + par1*dToImpact[i]) );
    return residual;
}

RVec <double> tofHit(const RVec <XYZTVector>& hits){
    RVec <double> tofs;
    int nHits = hits.size();
    for (int i = 0; i < nHits; ++i)
        tofs.push_back(hits[i].T());
    return tofs;
}


RVec <double> smear_time(const RVec<double> &t_hit){
    RVec <double> result;
        for(auto& t : t_hit) result.push_back(t + distribution(generator));
    return result;
}
