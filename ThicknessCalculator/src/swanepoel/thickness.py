# -*- coding: utf-8 -*-
"""
Thickness calculations (Swanepoel Eq. 23 and Eq. 3).
"""
import numpy as np
from scipy.stats import norm


def calculate_initial_thickness(lam_nm, n_vals):
    """
    Calculate film thickness using Swanepoel Eq. (23).
    
    Input:
        lam_nm (array): array of wavelengths at extrema [nm]
        n_vals : array of refractive indices at those wavelengths
    Returns:
        d_all : array of thickness values [m] computed from adjacent pairs
    """
    lam_m = lam_nm * 1e-9  # convert nm → m

    d1 = []
    for i in range(len(lam_m)-1):
            num = lam_m[i]*lam_m[i+1]
            den = (2*(lam_m[i]*n_vals[i+1]-lam_m[i+1]*n_vals[i]))
            if den != 0:
                
                d1.append( abs(num / den) )
            else:
                print("Error: Division with 0, when using Swanepoel Eq. (23)")
    return np.array(d1)

def thickness_estimate(lam_peaks_nm, n_peaks,
                       lam_valleys_nm, n_valleys,
                       d1, trials=3):
    """
    Calculate film thickness d2 using Swanepoel Eq. (3).
    
    Input:
        lam_peaks_nm (array): array of wavelengths at peaks [nm]
        n_peaks (array): Refractive indices at those wavelengths
        lam_valleys_nm (array): array of wavelengths at valleys [nm]
        n_valleys (array): Refractive indices at those wavelengths
    Returns:
        d_all : array of thickness values [m] computed from adjacent pairs
    """

    d1 = float(np.mean(d1))

    # Combine peaks (integer m) & valleys (half-integer m)
    lam_all_nm = np.concatenate([lam_peaks_nm, lam_valleys_nm])
    n_all   = np.concatenate([n_peaks,     n_valleys])
    type_flag = np.concatenate([np.ones_like(lam_peaks_nm),
                                np.zeros_like(lam_valleys_nm)])

    idx = np.argsort(lam_all_nm)
    lam_all_nm = lam_all_nm[idx]
    n_all = n_all[idx]
    type_flag = type_flag[idx]

    lam_all_m = lam_all_nm * 1e-9

    # Raw m estimate using d1
    m_raw = 2 * n_all * d1 / lam_all_m

    # First extremum’s raw m
    m0 = m_raw[0]

    # Define starting m trials, e.g., nearest integer +/- 0,1,2
    start_candidates = np.arange(np.floor(m0)-trials,
                                 np.floor(m0)+trials+1)

    best_d2 = None
    best_var = np.inf

    # Try several possible starting m-values
    for m_start in start_candidates:

        m_trial = np.zeros_like(m_raw)

        # Assign the first m:
        if type_flag[0] == 1:      # peak → integer
            m_trial[0] = m_start
        else:                      # valley → half-integer
            m_trial[0] = m_start + 0.5

        # ---- Assign subsequent m-values by following m_raw downwards ----
        for i in range(1, len(m_raw)):
            # expected m from d1 ratio:
            m_expected = m_raw[i] - m_raw[0] + m_trial[0]

            if type_flag[i] == 1:   # peak: nearest integer
                m_trial[i] = np.round(m_expected)
            else:                   # valley: nearest half-integer
                m_trial[i] = np.round(m_expected - 0.5) + 0.5

        # Compute d_i from eq. (3)
        d_i = m_trial * lam_all_m / (2 * n_all)

        #Check variance
        var = np.var(d_i)

        # Keep the best m-set
        if var < best_var:
            best_var = var
            best_d2 = d_i

    return best_d2

def summarize_thickness(d_m):
    """
    Compute summary statistics (mean and std) for thickness array.
    
    Input
        d_m : array of thickness values [m]
    Returns
        dict with:
        "mean_m" : float
        "std_m" : float
        "CI95" : float
    """
    #d_m = np.asarray(d_m, float)

    #Get mean and std
    mean = np.mean(d_m)
    std = np.std(d_m)
    
    #Calculate 95% confidence interval
    alpha = 0.05
    z = norm.ppf(1 - alpha/2)
    N = len(d_m)
    SEM = std / np.sqrt(N)
    ci95 = (mean - z*SEM, mean + z*SEM)

    return {
        "mean_m": mean,
        "std_m": std,
        "CI95": ci95}