# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 08:53:57 2025

@author: s224492
"""

from dataclasses import dataclass
import numpy as np

@dataclass
class Envelopes:
    lam_band_nm: np.ndarray
    TM: np.ndarray
    Tm: np.ndarray
    lam_peaks_nm: np.ndarray
    lam_valleys_nm: np.ndarray