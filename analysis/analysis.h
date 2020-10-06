using namespace std;
using namespace ROOT::VecOps;
double c = 299.792458;
struct Hit{
    double t;
    double r;
    double d;
    int layer;
    Hit (double t_hit, double r_hit, double d_hit, int layer_hit) :
        t(t_hit), r(r_hit), d(d_hit), layer(layer_hit){}
};

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
        t_err[i] = 10.;
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
    if (tof <= 0.) return 0.;
    return tof;
}
