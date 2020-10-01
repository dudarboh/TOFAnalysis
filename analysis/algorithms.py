from ROOT import TGraphErrors, TCanvas
import numpy as np
import time
# c in mm/ns
c = 299.792458


def algo_closest(p, l_trk, hits):
    hit = min(hits, key = lambda x: x["r"])
    tof = hit["t"] - hit["r"]/c
    if tof <= 0.:
        return 0.
    beta = l_trk/(tof*c)
    return p * np.sqrt(1. - beta*beta)/beta

def algo_avg(p, l_trk, hits):
    # Franks Average
    tofs = []
    for layer in range(10):
        layer_hits = [hit for hit in hits if hit["layer"] == layer]
        if layer_hits == []:
            continue
        hit = min(layer_hits, key=lambda x: x["d"])
        tofs.append(hit["t"] - hit["r"]/c)
    n_tofs = len(tofs)
    if n_tofs == 0:
        return 0.
    tof = sum(tofs)/n_tofs
    beta = l_trk/(tof*c)
    return p * np.sqrt(1. - beta*beta)/beta


def algo_fit(p, l_trk, hits):
    # Get Franks layer hits
    # canvas_debug = TCanvas()
    hits_selected = []
    for layer in range(10):
        layer_hits = [hit for hit in hits if hit["layer"] == layer]
        if layer_hits == []:
            continue
        hits_selected.append(min(layer_hits, key=lambda x: x["d"]))
    n_hits = len(hits_selected)
    if n_hits <= 1:
        return 0.

    r = np.array([hit["r"] for hit in hits_selected])
    t = np.array([hit["t"] for hit in hits_selected])
    r_err = np.zeros_like(r)
    t_err = np.full_like(t, 10.)
    gr = TGraphErrors(n_hits, r, t, r_err, t_err)

    gr.Fit("pol1", "Q")
    # gr.Draw("AP")
    # gr.SetMarkerStyle(20)
    # gr.SetTitle("Fit of a TOF vs distance to impact point; d, [mm];TOF, [ns]")
    # canvas_debug.Update()
    # time.sleep(3)
    # raw_input("wait")
    # if fit.GetChisquare() > 1.:
    #     continue
    fit = gr.GetFunction("pol1")
    tof = fit.GetParameter(0)
    if tof == 0:
        return 0.
    beta = l_trk/(tof*c)
    return p * np.sqrt(1. - beta*beta)/beta
