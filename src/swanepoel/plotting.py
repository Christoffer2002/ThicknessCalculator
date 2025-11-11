# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

def plot_envelopes(lam_nm, T, lam_band_nm, TM, Tm, lam_peaks_nm, lam_valleys_nm):
    fig, ax = plt.subplots()
    ax.plot(lam_nm, T, label="T (raw)", alpha=0.6)
    ax.plot(lam_band_nm, TM, label="Upper envelope TM")
    ax.plot(lam_band_nm, Tm, label="Lower envelope Tm")
    if len(lam_peaks_nm):
        ax.scatter(lam_peaks_nm, np.interp(lam_peaks_nm, lam_band_nm, TM), marker="^", label="peaks")
    if len(lam_valleys_nm):
        ax.scatter(lam_valleys_nm, np.interp(lam_valleys_nm, lam_band_nm, Tm), marker="v", label="valleys")
    ax.set_xlabel("Wavelength [nm]"); ax.set_ylabel("Transmittance")
    ax.legend(); ax.set_title("Swanepoel envelopes")
    return fig

def plot_n_band(lam_band_nm, n_band):
    fig, ax = plt.subplots()
    ax.plot(lam_band_nm, n_band, label="n(λ) film")
    ax.set_xlabel("Wavelength [nm]"); ax.set_ylabel("Refractive index n")
    ax.legend(); ax.set_title("Film refractive index (Eq. 11)")
    return fig

def plot_d_hist(d2_all_m):
    fig, ax = plt.subplots()
    ax.hist(d2_all_m * 1e6, bins="auto")
    ax.set_xlabel("Thickness [μm]"); ax.set_ylabel("Count")
    ax.set_title("Thickness distribution (d2)")
    return fig


def plot_fft(f, F1, F2, labels=("Measured","Theory")):
    fig, ax = plt.subplots()
    ax.plot(f, F1, label=labels[0])
    ax.plot(f, F2, label=labels[1])
    ax.set_xlabel("Spatial frequency [1/nm]")
    ax.set_ylabel("Amplitude")
    ax.grid(alpha=0.3)
    ax.legend()
    return fig

def plot_measured_vs_theoretical(lam_nm, T_meas, T_theory, title="Measured vs Theoretical"):
    """
    Plot measured and theoretical transmission spectra on the same axes.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(lam_nm, T_meas,   label="Measured", lw=1.3, color="black")
    ax.plot(lam_nm, T_theory, label="Theoretical", lw=1.5, color="tab:red")

    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Transmittance")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    return fig

