# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 08:53:57 2025

@author: s224492
"""
import numpy as np
from scipy.ndimage import uniform_filter1d
from .models import Envelopes



def find_extrema(signal, lam, min_sep_nm):
    """
    Return sorted indices of local maxima and minima in T(λ).
    
    Input:
        signal (array) : Transmission values T(λ) ranging from [0:1]
        lam (array) : Wavelengths
        min_sep_nm (float) : Minimum seperation length of each extrema in nm
    Return:
        max_idx (array) : array of index values where peaks occur
        min_idx (array) : array of index values where valleys occur
    """

    mid = np.arange(1, len(signal)-1)  # skip ends

    # raw extrema by neighbor test
    is_max = (signal[mid] > signal[mid-1]) & (signal[mid] > signal[mid+1])
    is_min = (signal[mid] < signal[mid-1]) & (signal[mid] < signal[mid+1])

    max_idx = mid[is_max]
    min_idx = mid[is_min]

    # wavelength-based thinning (greedy)
    def thin_by_lambda(idxs):
        if idxs.size == 0: return idxs
        kept = [idxs[0]]
        for i in idxs[1:]:
            if lam[i] - lam[kept[-1]] >= min_sep_nm:
                kept.append(i)
        return np.array(kept, dtype=int)

    max_idx = thin_by_lambda(np.sort(max_idx))
    min_idx = thin_by_lambda(np.sort(min_idx))

    return max_idx, min_idx


def fit_envelopes(
    lam_band_nm: np.ndarray,
    lam_peaks_nm: np.ndarray,
    T_peaks: np.ndarray,
    lam_valleys_nm: np.ndarray,
    T_valleys: np.ndarray,
    *,
    fit: str = "poly3",
    window_size: int = 5,
) -> Envelopes:
    """
    Within "min_lam" and "max_lam" detect and smooth extrema (using MA-filter).
    Fit 3'rd order polynomial to peaks and valleys creating upper
    and lower envelopes (TM & Tm).
    
    Input:
        T (array) : Transmittance (fraction in [0,1])
        lam (array) : Wavelength (nm)
        min_lam (float) : User-defined minimum wavelength for polynomial fit
        max_lam (float): User-defined maximum wavelength for polynomial fit
    Returns:
        lam_band (array) : Wavelength (nm) cut  to [min_lam:max_lam]
        TM_band (array) : Upper Swanepoel envelope
        Tm_band (array) : Lower SwanePoel envelope
        lam_peaks (array) : Wavelength (nm) where peaks are located
        lam_valleys (array) : Wavelength (nm) where valleys are located
    """
    # Smooth peaks and valleys
    T_peaks = uniform_filter1d(T_peaks, size=window_size, mode="nearest")
    T_valleys = uniform_filter1d(T_valleys, size=window_size, mode="nearest")
    
    # calculate envelopes
    # lower envelope
    if lam_valleys_nm.size >= 4:
        # Tm envelope - 3'rd order fit
        Tm_coef = np.polyfit(lam_valleys_nm, T_valleys, deg=3)
        poly_valley = np.poly1d(Tm_coef)
        Tm_band = poly_valley(lam_band_nm)
    else:
        Tm_band = np.interp(lam_band_nm, lam_valleys_nm, T_valleys)
    
    # upper envelope
    if lam_peaks_nm.size >= 4:
        TM_coef = np.polyfit(lam_peaks_nm, T_peaks, deg=3)
        poly_peak = np.poly1d(TM_coef)
        TM_band = poly_peak(lam_band_nm)
    else:
        TM_band = np.interp(lam_band_nm, lam_peaks_nm, T_peaks)
    
    return Envelopes(
    lam_band_nm=lam_band_nm,
    TM=TM_band,
    Tm=Tm_band,
    lam_peaks_nm=lam_peaks_nm,
    lam_valleys_nm=lam_valleys_nm,
)