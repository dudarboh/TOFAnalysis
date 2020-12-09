#include "analysis.hpp"

double tof_closest(const RVec<double> &t_hit, const RVec<double> &r_hit){
    int min_idx = min_element(r_hit.begin(), r_hit.end()) - r_hit.begin();
    // double tof = t_hit[min_idx];
    double tof = t_hit[min_idx] - r_hit[min_idx]/c;
    if (tof <= 0.) return 0.;
    return tof;
}

double tof_fastest(const RVec<double> &t_hit, const RVec<double> &r_hit){
    int min_idx = min_element(t_hit.begin(), t_hit.end()) - t_hit.begin();
    // double tof = t_hit[min_idx];
    double tof = t_hit[min_idx] - r_hit[min_idx]/c;
    if (tof <= 0.) return 0.;
    return tof;
}

// Franks Average
double tof_avg(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit, int n_layers=10){
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }
    int n_tofs = 0;
    double sum_tof = 0.;
    for (int layer=0; layer<n_layers; ++layer){
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

double tof_fit(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit, int n_layers=10, int func=1, double res_cut=999.){
    // Get Franks layer hits
    int n_hits = t_hit.size();
    map <int, vector <Hit> > layer_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.layer < n_layers) layer_hits[hit.layer].push_back(hit);
    }

    vector<Hit> selected_hits;
    for (int layer=0; layer<n_layers; ++layer){
        vector <Hit>& arr = layer_hits[layer];
        if (arr.size() == 0) continue;
        Hit closest_hit = *min_element(layer_hits[layer].begin(), layer_hits[layer].end(), [](const Hit &hit1, const Hit &hit2){return hit1.d < hit2.d;} );
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

    string fit_func = "pol1";
    if (func == 2) fit_func = "pol2";
    gr.Fit(fit_func.c_str(), "Q");
    TF1 fit = *gr.GetFunction(fit_func.c_str());

    //After the fit check residual for each hit and exclude hits with too high
    // residual and refit
    vector<double> r_new;
    vector<double> t_new;
    vector<double> r_err_new;
    vector<double> t_err_new;

    for(int i=0; i < n_sel_hits; ++i){
        double residual = abs(t[i] - fit.Eval(r[i]));
        if(residual > res_cut) continue;
        r_new.push_back(r[i]);
        t_new.push_back(t[i]);
        r_err_new.push_back(0.);
        t_err_new.push_back(0.01);
    }
    if (r_new.size() == 0) return -999;
    TGraphErrors gr2(r_new.size(), &r_new[0], &t_new[0], &r_err_new[0], &t_err_new[0]);
    gr2.Fit(fit_func.c_str(), "Q");
    TF1 fit2 = *gr2.GetFunction(fit_func.c_str());


    // if (fit.GetChisquare() > 1.) continue;
    double tof = fit2.GetParameter(0);
    if (tof <= 0.) return -999.;
    return tof;
}

double tof_fit_cyl(const RVec<double> &t_hit, const RVec<double> &r_hit, const RVec<double> &d_hit, const RVec<int> &layer_hit, int n_layers=10, int func=1, double res_cut=999.){
    int n_hits = t_hit.size();
    vector<Hit> selected_hits;
    for (int i=0; i < n_hits; ++i){
        Hit hit(t_hit[i], r_hit[i], d_hit[i], layer_hit[i]);
        if (hit.d > 4) continue;
        selected_hits.push_back(hit);
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
        t_err[i] = 0.1;
    }
    TGraphErrors gr(n_sel_hits, r, t, r_err, t_err);

    string fit_func = "pol1";
    if (func == 2) fit_func = "pol2";
    gr.Fit(fit_func.c_str(), "Q");
    TF1 fit = *gr.GetFunction(fit_func.c_str());


    //After the fit check residual for each hit and exclude hits with too high
    // residual and refit
    vector<double> r_new;
    vector<double> t_new;
    vector<double> r_err_new;
    vector<double> t_err_new;

    for(int i=0; i < n_sel_hits; ++i){
        double residual = abs(t[i] - fit.Eval(r[i]));
        if(residual > res_cut) continue;
        r_new.push_back(r[i]);
        t_new.push_back(t[i]);
        r_err_new.push_back(0.);
        t_err_new.push_back(0.01);
    }
    if (r_new.size() == 0) return -999;
    TGraphErrors gr2(r_new.size(), &r_new[0], &t_new[0], &r_err_new[0], &t_err_new[0]);
    gr2.Fit(fit_func.c_str(), "Q");
    TF1 fit2 = *gr2.GetFunction(fit_func.c_str());

    // if (fit.GetChisquare() > 1.) continue;
    double tof = fit2.GetParameter(0);
    if (tof <= 0.) return -999.;
    return tof;
}
