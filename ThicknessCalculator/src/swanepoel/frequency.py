# -*- coding: utf-8 -*-
"""
Reconstruct the theoretical signal using d2. 

"""
import numpy as np

def theoretical_spectrum(lam_nm, n, s, d2_all_m):
    """
    Compute theoretical T(λ) using Swanepoel Appendix A1, k=0.

    Input:
        lam_nm (array) : Wavelength (nm).
        n (array) : n(λ) - refractive index dependent on lambda.
    Return:
        T (array) : Theoretical transmittance (fraction in [0,1]).
    """
    lam_m = lam_nm * 1e-9
    d_m = float(np.nanmean(d2_all_m))

    #k = 0.0
    x = 1.0  # transparent film
    phi = 4.0 * np.pi * n * d_m / lam_m

    A = 16.0 * s * (n**2)
    B = ((n + 1)**2) * ((n + 1)*(n + s**2))
    C = ( (n**2 - 1)*(n**2 - s**2) ) * 2*np.cos(phi)
    D = ((n - 1)**2) * ((n - 1)*(n - s**2))

    T = (A * x) / (B - C * x + D * x**2)
    return np.clip(T, 0.0, 1.0)


def fft_compare(lam_nm, T1, T2):
    """
    Resample onto uniform λ-grid and return FFT magnitudes.

    Input:
        lam_nm (array) : Wavelength (nm).
        T1 (array) : Recorded transmittance (fraction in [0,1]).
        T2 (array) : Theoretical transmittance (fraction in [0,1]).
    Return:
        f (array) : frequency axis
        F1 (array) : Recorded frequency components
        F2 (array) : Theoretical frequency components
    """ 
    step = float(np.median(np.diff(lam_nm))) # Step size (in nm)

    lam_u = np.arange(lam_nm.min(), lam_nm.max()+step/2, step)
    #T1_u = np.interp(lam_u, lam_nm, T1)
    #T2_u = np.interp(lam_u, lam_nm, T2)

    Fs = 1.0 / step # Sampling frequency - samples pr nm
    # Frequency axis
    f = np.fft.rfftfreq(len(lam_nm), d=1/Fs)

    #Hanning window to account for spectral leakage
    w = np.hanning(len(lam_u))
    F1 = np.abs(np.fft.rfft((T1 - T1.mean()) * w))
    F2 = np.abs(np.fft.rfft((T2 - T2.mean()) * w))
    # No Hanning window
    # F1 = np.abs(np.fft.rfft((T1 - T1.mean())))
    # F2 = np.abs(np.fft.rfft((T2 - T2.mean()))) 


    return f, F1, F2

