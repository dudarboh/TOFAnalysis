using std::cout, std::endl, std::stringstream, std::vector, std::string, std::runtime_error;
using namespace ROOT::VecOps;

#include <random>
#include "DDRec/Vector3D.h"
using dd4hep::rec::Vector3D;
const double PI = 3.14159265;
const double rInner = 1804.8;
const double c = 299.792458;
std::default_random_engine generator;
std::normal_distribution <double> distribution(0.0, 0.05);

int _image = 0;

struct Hit{
    double t;
    double r;
    double d;
    int layer;
    Hit (double t_hit, double r_hit, double d_hit, int layer_hit) :
        t(t_hit), r(r_hit), d(d_hit), layer(layer_hit){}
};


double dToRef_fastest(const RVec<double> &t_hit, const RVec<double> &r_hit){
    int min_idx = min_element(t_hit.begin(), t_hit.end()) - t_hit.begin();
    double dToRef = r_hit[min_idx];
    if (dToRef <= 0.) return 0.;
    return dToRef;
}

double dToLine_fastest(const RVec<double> &t_hit, const RVec<double> &d_hit){
    int min_idx = min_element(t_hit.begin(), t_hit.end()) - t_hit.begin();
    double dToLine = d_hit[min_idx];
    if (dToLine <= 0.) return 0.;
    return dToLine;
}

double layer_fastest(const RVec<double> &t_hit, const RVec<double> &layer_hit){
    int min_idx = min_element(t_hit.begin(), t_hit.end()) - t_hit.begin();
    double layer = layer_hit[min_idx];
    if (layer <= 0.) return 0.;
    return layer;
}


RVec <double> smear_time(const RVec<double> &t_hit){
    RVec <double> result;
        for(auto& t : t_hit) result.push_back(t + distribution(generator));
    return result;
}




RVec <double> dToLine(const RVec <float>& x, const RVec <float>& y, const RVec <float>& z,
                     const float& xRef, const float& yRef, const float& zRef,
                     const float& px, const float& py, const float& pz){
    RVec <double> distance;

    Vector3D refPoint(xRef, yRef, zRef);
    Vector3D unitDir = Vector3D(px, py, pz).unit();
    for(int i = 0; i < x.size(); ++i){
        Vector3D caloHit(x[i], y[i], z[i]);
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
    int nSides = 8;
    double step = PI/nSides;
    // phi is in range[-pi, pi]. 1 plane will lay between at pi.
    if(phi < (-PI + step) || phi > (PI - step)){
        Vector3D pos(1804.8, PI, PI/2., Vector3D::spherical);
        return pos[dir];

    }
    double center_phi = -PI + 2.*step;
    for (int i = 0; i < nSides-1; ++i){
        if (center_phi - step <= phi && phi < center_phi + step) break;
        else center_phi += 2.*step;
    }
    Vector3D pos(1804.8, center_phi, PI/2., Vector3D::spherical);
    return pos[dir];
}


double getLineDir(double r, double phi, int theta, int dir){
    //Line direction based on cluster direction
    Vector3D pos(r, phi, theta, Vector3D::spherical);
    return pos[dir];
}


int nCaloHits(const RVec<int> &layer_hit, int maxLayer){
    int n_hits = 0;
    for(auto& l : layer_hit){
        if (l < maxLayer) ++n_hits;
    }
    return n_hits;
}


RVec <double> time_hit_fit(const RVec<double> &t_hit, const RVec<double> &d_hit, const RVec<double> &layer_hit, int n_layers=10){
    // Get Franks layer hits
    int n_hits = d_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], 0., d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }

    RVec<double> time_arr;
    for (int layer=0; layer<n_layers; ++layer){
        vector<Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(arr.begin(), arr.end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        time_arr.push_back(closest_hit.t);
    }
    return time_arr;
}


RVec <double> dToLine_hit_fit(const RVec<double> &d_hit, const RVec<double> &layer_hit, int n_layers=10){
    // Get Franks layer hits
    int n_hits = d_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(0., 0., d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }

    RVec<double> dToLine_arr;
    for (int layer=0; layer<n_layers; ++layer){
        vector<Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(arr.begin(), arr.end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        dToLine_arr.push_back(closest_hit.d);
    }
    return dToLine_arr;
}


RVec <double> dToRef_hit_fit(const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<double> &layer_hit, int n_layers=10){
    // Get Franks layer hits
    int n_hits = d_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(0., r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }

    RVec<double> dToRef_arr;
    for (int layer=0; layer<n_layers; ++layer){
        vector<Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(arr.begin(), arr.end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        dToRef_arr.push_back(closest_hit.r);
    }
    return dToRef_arr;
}

RVec <int> layer_hit_fit(const RVec<double> &d_hit, const RVec<double> &layer_hit, int n_layers=10){
    // Get Franks layer hits
    int n_hits = d_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(0., 0., d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }

    RVec<int> layer_arr;
    for (int layer=0; layer<n_layers; ++layer){
        vector<Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(arr.begin(), arr.end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        layer_arr.push_back(closest_hit.layer);
    }
    return layer_arr;
}

int fit_plot(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit, int n_layers=10){
    // Get Franks layer hits
    TCanvas canvas = TCanvas();
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }
    // TCanvas canvas_debug;

    vector<Hit> selected_hits;
    for (int layer=0; layer<n_layers; ++layer){
        vector <Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(layer_hits[layer].begin(), layer_hits[layer].end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        selected_hits.push_back(closest_hit);
    }

    const int n_sel_hits = selected_hits.size();
    if (n_sel_hits <= 1) return true;

    double r[n_sel_hits];
    double t[n_sel_hits];
    double r_err[n_sel_hits];
    double t_err[n_sel_hits];
    for(int i=0; i < n_sel_hits; ++i){
        r[i] = selected_hits[i].r;
        t[i] = selected_hits[i].t;
        r_err[i] = 0.;
        t_err[i] = 0.1;
    }
    TGraphErrors gr(n_sel_hits, r, t, r_err, t_err);

    gr.Fit("pol1", "Q");
    TF1 fit = *gr.GetFunction("pol1");
    gr.SetTitle("Fit with >10ps bias; |#vec{r}_{hit} - #vec{r}_{impact}|, mm; TOF_{hit}, ns");
    gr.Draw("AP");
    canvas.Print(Form("./pics/%i_test.png", _image));
    cout<<_image++<<endl;
    return true;
}


double tof_chi2(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit, int n_layers=10, double d_cut=9999.){
    // Get Franks layer hits
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }
    // TCanvas canvas_debug;

    vector<Hit> selected_hits;
    for (int layer=0; layer<n_layers; ++layer){
        vector <Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(layer_hits[layer].begin(), layer_hits[layer].end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        if (closest_hit.d > d_cut) continue;
        selected_hits.push_back(closest_hit);
    }

    const int n_sel_hits = selected_hits.size();
    if (n_sel_hits <= 1) return -999.;

    double r[n_sel_hits];
    double t[n_sel_hits];
    double r_err[n_sel_hits];
    double t_err[n_sel_hits];
    for(int i=0; i < n_sel_hits; ++i){
        r[i] = selected_hits[i].r;
        t[i] = selected_hits[i].t;
        r_err[i] = 0.;
        t_err[i] = 0.01;
    }
    TGraphErrors gr(n_sel_hits, r, t, r_err, t_err);

    gr.Fit("pol1", "Q");
    TF1 fit = *gr.GetFunction("pol1");
    return fit.GetChisquare()/fit.GetNDF();
}

RVec<double> residual(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit, int n_layers=10, double d_cut=9999.){
    // Get Franks layer hits
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }
    // TCanvas canvas_debug;

    vector<Hit> selected_hits;
    for (int layer=0; layer<n_layers; ++layer){
        vector <Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(layer_hits[layer].begin(), layer_hits[layer].end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
        if (closest_hit.d > d_cut) continue;
        selected_hits.push_back(closest_hit);
    }

    const int n_sel_hits = selected_hits.size();
    if (n_sel_hits <= 1) return {};

    double r[n_sel_hits];
    double t[n_sel_hits];
    double r_err[n_sel_hits];
    double t_err[n_sel_hits];
    for(int i=0; i < n_sel_hits; ++i){
        r[i] = selected_hits[i].r;
        t[i] = selected_hits[i].t;
        r_err[i] = 0.;
        t_err[i] = 0.01;
    }
    TGraphErrors gr(n_sel_hits, r, t, r_err, t_err);

    gr.Fit("pol1", "Q");
    TF1 fit = *gr.GetFunction("pol1");

    RVec<double> residuals;

    for(int i=0; i < n_sel_hits; ++i) residuals.push_back(t[i] - fit.Eval(r[i]));

    return residuals;
}

////////////////////// There are so many things on top, so I start a new section where I try to be more pedantic with my code *_*

RVec <Hit> class_test(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit){
    int n_hits = t_hit.size();
    RVec <Hit> hits;
    for (int i=0; i < n_hits; ++i){
        hits.push_back(Hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]));
    }
    return hits;
}
