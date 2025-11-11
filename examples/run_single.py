# examples/run_single.py
import numpy as np
import matplotlib.pyplot as plt
from swanepoel.io import load_spectrum_csv, bandpass
from swanepoel.extrema import find_extrema, fit_envelopes
from swanepoel.optics import substrate_refractive_index, film_refractive_index
from swanepoel.thickness import calculate_initial_thickness, thickness_estimate, summarize_thickness
from swanepoel.plotting import plot_envelopes, plot_n_band, plot_d_hist
from swanepoel.frequency import theoretical_spectrum, fft_compare
from swanepoel.plotting import plot_fft, plot_measured_vs_theoretical

#Run from ThicknessCalculator folder - "pip install -e ."
# load data
lam, T = load_spectrum_csv("data/GT-Thickness/Square1_SpotA_Rep1.csv")

# crop signal
lam_b, T_b = bandpass(lam, T, lam_min=600, lam_max=900)

# extrema + envelopes
pk_idx, vl_idx = find_extrema(T_b, lam_b, min_sep_nm=2.0)
lam_pk, T_pk = lam_b[pk_idx], T_b[pk_idx]
lam_vl, T_vl = lam_b[vl_idx], T_b[vl_idx]
env = fit_envelopes(lam_b, lam_pk, T_pk, lam_vl, T_vl, window_size=5)

# substrate + film n(λ)
s = substrate_refractive_index(env.lam_band_nm, "cauchy", (1.5690, 0.00531))
n_band = film_refractive_index(env.TM, env.Tm, s)

# interpolate n at extrema
n_pk = np.interp(lam_pk, env.lam_band_nm, n_band)
n_vl = np.interp(lam_vl, env.lam_band_nm, n_band)

# thickness (Eq.23) → scalar d1 → refine d2
d1_pk = calculate_initial_thickness(lam_pk, n_pk)
d1_vl = calculate_initial_thickness(lam_vl, n_vl)
d1_all = np.concatenate((d1_pk, d1_vl))

d2_all = thickness_estimate(lam_pk, n_pk, lam_vl, n_vl, d1_all)

# summary
stats = summarize_thickness(d2_all)
print(f"Mean thickness: {stats['mean_m']*1e6:.1f} μm")
print(f"Std: {stats['std_m']*1e9:.1f} nm")
print(f"95% CI: [{stats['CI95'][0]*1e9:.1f}, {stats['CI95'][1]*1e9:.1f}] nm")

# Plots
plot_envelopes(lam, T, env.lam_band_nm, env.TM, env.Tm, lam_pk, lam_vl)
plot_n_band(env.lam_band_nm, n_band)
plot_d_hist(d2_all)


# Theoretical spectrum
T_theory = theoretical_spectrum(env.lam_band_nm, n_band, s, d2_all)

# FFT comparison
f, F_meas, F_theory = fft_compare(env.lam_band_nm, T_b, T_theory)
plot_fft(f, F_meas, F_theory)

# Plot in time domaine
T_theory = theoretical_spectrum(env.lam_band_nm, n_band, s, d2_all)
T_meas_band = np.interp(env.lam_band_nm, lam, T)
plot_measured_vs_theoretical(env.lam_band_nm, T_meas_band, T_theory)

plt.show()