# -*- coding: utf-8 -*-
"""
Bandpass filter to cut signal to the length being examined.
"""
import numpy as np
import pandas as pd

def bandpass(lam_nm, T, lam_min, lam_max):
    m = (lam_nm >= lam_min) & (lam_nm <= lam_max)
    return lam_nm[m].astype(float), T[m].astype(float)

def load_spectrum_csv(path):
    """
    Load transmission spectrum file.

    Input:
        path (string) : String containing path to file
    Return:
        lam (array) : Wavelength (nm)
        T (array) : Transmittance (fraction in [0,1])
    """
    
    # Autodetect delimiter
    df = pd.read_csv(path, sep=None, engine="python")
    # try:
    #     df = pd.read_csv(path, sep=None, engine="python")
    # except Exception:
    #     df = pd.read_csv(path, sep=r"[;\t,]+", engine="python", header=None)

    x = df.iloc[:, 0].astype(str)
    y = df.iloc[:, 1].astype(str)

    # Clean up: remove % and normalize decimal commas to dots
    x = x.str.replace("%", "", regex=False).str.replace(",", ".", regex=False)
    y = y.str.replace("%", "", regex=False).str.replace(",", ".", regex=False)

    # Wavelengths and transmittance
    lam = pd.to_numeric(x, errors="coerce").to_numpy()
    T   = pd.to_numeric(y, errors="coerce").to_numpy()
    
    # Convert % to fraction
    T = T / 100.0

    # Filter out potential NaN's
    m = np.isfinite(lam) & np.isfinite(T)
    lam, T = lam[m], T[m]
    idx = np.argsort(lam)
    lam, T = lam[idx], T[idx]

    return lam, T