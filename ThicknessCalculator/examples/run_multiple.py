# examples/run_multiple.py
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

from swanepoel.io import load_spectrum_csv, bandpass
from swanepoel.extrema import find_extrema, fit_envelopes
from swanepoel.optics import substrate_refractive_index, film_refractive_index
from swanepoel.thickness import calculate_initial_thickness, thickness_estimate, summarize_thickness
from swanepoel.frequency import theoretical_spectrum, fft_compare
from swanepoel.plotting import (
    plot_envelopes, plot_n_band, plot_d_hist,
    plot_fft, plot_measured_vs_theoretical
)

# Settings
folder = "data/GT-Thickness/"
pattern = "Square1_*.csv"   # process all CSVs in folder
lam_min = 600
lam_max = 900
substrate_model = "cauchy"
substrate_coeffs = (1.5690, 0.00531)

files = sorted(glob.glob(os.path.join(folder, pattern)))
print(f"Found {len(files)} files.")

all_results = []   # list to store each sample’s results
lengths = []

for path in files:
    print("\n--- Processing:", os.path.basename(path), "---")

    # === SAME CODE AS run_single ===
    lam, T = load_spectrum_csv(path)
    lam_b, T_b = bandpass(lam, T, lam_min=lam_min, lam_max=lam_max)

    pk_idx, vl_idx = find_extrema(T_b, lam_b, min_sep_nm=2.0)
    lam_pk, T_pk = lam_b[pk_idx], T_b[pk_idx]
    lam_vl, T_vl = lam_b[vl_idx], T_b[vl_idx]
    env = fit_envelopes(lam_b, lam_pk, T_pk, lam_vl, T_vl, window_size=5)

    s = substrate_refractive_index(env.lam_band_nm, substrate_model, substrate_coeffs)
    n_band = film_refractive_index(env.TM, env.Tm, s)

    n_pk = np.interp(lam_pk, env.lam_band_nm, n_band)
    n_vl = np.interp(lam_vl, env.lam_band_nm, n_band)

    d1_pk = calculate_initial_thickness(lam_pk, n_pk)
    d1_vl = calculate_initial_thickness(lam_vl, n_vl)
    d1_all = np.concatenate((d1_pk, d1_vl))

    d2_all = thickness_estimate(lam_pk, n_pk, lam_vl, n_vl, d1_all)
    
    stats = summarize_thickness(d2_all)

    all_results.append((os.path.basename(path), stats))

    print(f"Mean thickness: {stats['mean_m']*1e6:.2f} μm")

    # PLOTTING
    # plot_envelopes(lam, T, env.lam_band_nm, env.TM, env.Tm, lam_pk, lam_vl)
    # plot_n_band(env.lam_band_nm, n_band)
    # plot_d_hist(d2_all)

    # THEORY + FFT
    # T_theory = theoretical_spectrum(env.lam_band_nm, n_band, s, d2_all)
    # T_meas_band = np.interp(env.lam_band_nm, lam, T)

    # plot_measured_vs_theoretical(env.lam_band_nm, T_meas_band, T_theory)
    # f, F_meas, F_theory = fft_compare(env.lam_band_nm, T_meas_band, T_theory)
    # plot_fft(f, F_meas, F_theory)

plt.show()

# ---- Summary table ----
print("\n===== SUMMARY OF ALL SAMPLES =====")
for fname, st in all_results:
    print(f"{fname:40s}  {st['mean_m']*1e6:7.2f} µm ± {st['std_m']*1e9:5.1f} nm")


