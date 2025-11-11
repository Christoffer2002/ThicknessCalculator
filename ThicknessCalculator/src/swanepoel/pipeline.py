# -*- coding: utf-8 -*-
# src/swanepoel/pipeline.py
import numpy as np
from .io import bandpass
from .extrema import find_extrema, fit_envelopes
from .optics import substrate_refractive_index, film_refractive_index
from .thickness import calculate_initial_thickness, thickness_estimate, summarize_thickness

def run_swanepoel(lam_nm, T,
                  lam_min, lam_max,
                  min_sep_nm, window_size,
                  substrate_model, substrate_coeffs):
    # 1) crop
    lam_b, T_b = bandpass(lam_nm, T, lam_min, lam_max)

    # 2) extrema + envelopes
    max_idx, min_idx = find_extrema(T_b, lam_b, min_sep_nm)
    lam_peaks   = lam_b[max_idx];   T_peaks   = T_b[max_idx]
    lam_valleys = lam_b[min_idx];   T_valleys = T_b[min_idx]
    env = fit_envelopes(lam_b, lam_peaks, T_peaks, lam_valleys, T_valleys, window_size=window_size)

    # 3) substrate index on band + Eq.11 at band
    s_band = substrate_refractive_index(env.lam_band_nm, substrate_model, substrate_coeffs)
    n_band = film_refractive_index(env.TM, env.Tm, s_band)

    # sample n at extrema (linear interp on the band result)
    n_peaks   = np.interp(lam_peaks,   env.lam_band_nm, n_band)
    n_valleys = np.interp(lam_valleys, env.lam_band_nm, n_band)

    # 4) thickness from Eq.23 (use peaks OR valleys; pick peaks if available)
    d1_src_lam = lam_peaks if len(lam_peaks) >= 2 else lam_valleys
    d1_src_n   = n_peaks   if len(lam_peaks) >= 2 else n_valleys
    d1_all = calculate_initial_thickness(d1_src_lam, d1_src_n)

    # 5) (optional) order-refined d2 — you said it’s fine to keep it
    d2_all = thickness_estimate(lam_peaks, n_peaks, lam_valleys, n_valleys, d1_all)

    # 6) summary
    stats = summarize_thickness(d2_all)

    return {
        "lam_band_nm": env.lam_band_nm,
        "TM": env.TM, "Tm": env.Tm,
        "lam_peaks_nm": lam_peaks, "n_peaks": n_peaks,
        "lam_valleys_nm": lam_valleys, "n_valleys": n_valleys,
        "d1_all_m": d1_all,
        "d2_all_m": d2_all,
        "summary": stats,
    }
