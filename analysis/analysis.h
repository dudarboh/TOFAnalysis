using std::cout, std::endl, std::stringstream, std::vector, std::string, std::runtime_error;
using namespace ROOT::VecOps;

#include "DDRec/Vector3D.h"
using dd4hep::rec::Vector3D;

const double PI = 3.141592653589793238462643383279502884;
const double rInner = 1804.8;
const Vector3D testVec(1804.8, 0., PI/2., Vector3D::spherical);
const double c = 299.792458;

struct Hit{
    double t;
    double r;
    double d;
    int layer;
    Hit (double t_hit, double r_hit, double d_hit, int layer_hit) :
        t(t_hit), r(r_hit), d(d_hit), layer(layer_hit){}
};


// Almost obsolete methods
double line(double x, double x1, double x2, double y1, double y2){
    return y1 + (y2-y1)/(x2-x1) * (x - x1);
}

double track_len(double phi1, double phi2, double omega1, double omega2, double lambda1, double lambda2){
    int nPoints = 10000;
    double h = abs(phi2-phi1)/nPoints;
    // Calculate func values at points
    double func[nPoints+1];
    for (int i = 0; i < nPoints+1; ++i) {
        double phi = phi1 + 1.*i/nPoints*(phi2-phi1);
        double omega = line(phi, phi1, phi2, omega1, omega2);
        double lambda = line(phi, phi1, phi2, lambda1, lambda2);
        func[i] = abs(1./omega) * sqrt(1. + lambda*lambda);
    }

    // Calculate integral by Simpsons formula
    double result = 0.;
    for (int i = 1; i < nPoints; i+=2) result += func[i-1] + 4*func[i] + func[i+1];
    return h/3. * result;
}

RVec <double> rFunc(const RVec <float>& x, const RVec <float>& y, const RVec <float>& z, float x0, float y0, float z0){
    RVec <double> rVec;
    for (size_t i = 0; i < x.size(); ++i) {
        rVec.push_back(sqrt((x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) + (z[i]-z0)*(z[i]-z0)));
    }
    return rVec;
}


double tof_closest(const RVec<double> &t_hit, const RVec<double> &r_hit){
    int min_idx = min_element(r_hit.begin(), r_hit.end()) - r_hit.begin();
    double tof = t_hit[min_idx] - r_hit[min_idx]/c;
    if (tof <= 0.) return 0.;
    return tof;
}

double tof_fastest(const RVec<double> &t_hit, const RVec<double> &r_hit){
    int min_idx = min_element(t_hit.begin(), t_hit.end()) - t_hit.begin();
    double tof = t_hit[min_idx] - r_hit[min_idx]/c;
    if (tof <= 0.) return 0.;
    return tof;
}

// Franks Average
double tof_avg(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit){
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        layer_hits[hit.layer].push_back(hit);
    }
    int n_tofs = 0;
    double sum_tof = 0.;
    for (int layer=0; layer<10; ++layer){
        vector <Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(layer_hits[layer].begin(), layer_hits[layer].end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        double closest_time = closest_hit.t - closest_hit.r/c;
        sum_tof += closest_time;
        ++n_tofs;
    }
    double tof = sum_tof/n_tofs;
    if (tof <= 0.) return 0.;
    return tof;
}


RVec <double> dToLine(const RVec <float>& x, const RVec <float>& y, const RVec <float>& z,
                     const float& xRef, const float& yRef, const float& zRef,
                     const float& phi, const float& theta){
    vector <double> distance;

    float posRef[3] = {xRef, yRef, zRef};
    Vector3D refPoint(posRef);
    Vector3D unitDir(1., phi, theta, Vector3D::spherical);

    for(int i = 0; i < x.size(); ++i){
        float pos[3] = {x[i], y[i], z[i]};
        Vector3D caloHit(pos);
        double d = (caloHit - refPoint).cross(unitDir).r();
        distance.push_back(d);
    }
    return distance;
}


RVec <double> dToRef(const RVec <float>& x, const RVec <float>& y, const RVec <float>& z,
                     const float& xRef, const float& yRef, const float& zRef){
    vector <double> distance;

    float posRef[3] = {xRef, yRef, zRef};
    Vector3D refPoint(posRef);

    for(int i = 0; i < x.size(); ++i){
        float pos[3] = {x[i], y[i], z[i]};
        Vector3D caloHit(pos);
        double d = (caloHit - refPoint).r();
        distance.push_back(d);
    }
    return distance;
}


double tof_fit(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit){
    // Get Franks layer hits
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        layer_hits[hit.layer].push_back(hit);
    }
    // TCanvas canvas_debug;

    vector<Hit> selected_hits;
    for (int layer=0; layer<10; ++layer){
        vector <Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(layer_hits[layer].begin(), layer_hits[layer].end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        selected_hits.push_back(closest_hit);
    }

    const int n_sel_hits = selected_hits.size();
    if (n_sel_hits <= 1) return 0.;

    double r[n_sel_hits];
    double t[n_sel_hits];
    double r_err[n_sel_hits];
    double t_err[n_sel_hits];
    for(int i=0; i < n_sel_hits; ++i){
        r[i] = selected_hits[i].r;
        t[i] = selected_hits[i].t;
        r_err[i] = 0.;
        t_err[i] = 0.05;
    }
    TGraphErrors gr(n_sel_hits, r, t, r_err, t_err);

    gr.Fit("pol1", "Q");
    // gr.Draw("AP");
    // gr.SetMarkerStyle(20);
    // gr.SetTitle("Fit of a TOF vs distance to impact point; d, [mm];TOF, [ns]");
    // canvas_debug.Update();
    // sleep(3);
    TF1 fit = *gr.GetFunction("pol1");
    // if (fit.GetChisquare() > 1.) continue;
    double tof = fit.GetParameter(0);
    if (tof <= 0.) return -999.;
    return tof;
}



// Code to find intersection of a line and a plane
double intersection(double nx, double ny, double nz,
                    double x0, double y0, double z0,
                    double px, double py, double pz,
                    double x1, double y1, double z1){
    //n - normal vector to the plane
    Vector3D n(nx, ny, nz);
    //x0 - any point on a plane
    Vector3D r0(x0, y0, z0);
    //p - directional vector of a line
    Vector3D p(px, py, pz);
    //x1 - any point on a line
    Vector3D r1(x1, y1, z1);

    // If plane and line parallel return 0.
    if (p.dot(n) == 0.) return 0.;

    // Line parameter for a solution
    double t = (r0 - r1).dot(n)/p.dot(n);
    return t;
}


double getPlaneDir(double phi, int dir){
    //Input is cluster direction.
    // This point gives direction AND point of the plane
    double center_phi = 0.;
    for (int i = 0; i < 8; ++i){
        if (center_phi - PI/8. <= phi < center_phi + PI/8.) break;
        else center_phi += PI/4.;
    }
    Vector3D pos(1804.8, center_phi, PI/2., Vector3D::spherical);
    return pos[dir];
}

double getLineDir(double r, double phi, int theta, int dir){
    //Line direction based on cluster direction
    Vector3D pos(r, phi, theta, Vector3D::spherical);
    return pos[dir];
}
