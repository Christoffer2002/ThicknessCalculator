# -*- coding: utf-8 -*-
"""
Refractive-index models + Swanepoel Eq. 11.
"""
import numpy as np

def substrate_refractive_index(lam_nm, model, coeffs):
    """
    Evaluate substrate refractive index s(λ) at wavelengths lam_nm (nm) using
    3-term Sellmeier with λ in µm and C in µm², 2 term Cauchy model or constant
    refractive index.
    
    Input:
        lam_nm (array) : Wavelength (nm)
        coeffs (array) : sellmeier or cauchy coeffecients [B1,C1,B2,C2,B3,C3], with C in µm²
        model (string) : "sellmeier" "cauchy_alt" or "cauchy".
    Returns:
        s(λ) (array) : Substrate refractive index
    """
    
    lam_um = np.asarray(lam_nm, float) * 1e-3 #Convert to μm
    L2 = lam_um**2
    model = model.lower()

    if model == "sellmeier":
        B1, C1, B2, C2, B3, C3 = np.asarray(coeffs, float)
        
        # sellmeier: n^2(λ) = 1 + Σ (B_i λ^2)/(λ^2 - C_i)
        n2 = 1.0 + (B1*L2)/(L2 - C1) + (B2*L2)/(L2 - C2) + (B3*L2)/(L2 - C3) 
        return np.sqrt(n2)

    if model == "cauchy":  # classic
        #A, B, C = np.asarray(coeffs, float)
        A, B = np.asarray(coeffs, float)
        return A + B*(lam_um**-2) # cauchy: n(λ) = A + B/λ^2

    if model == "cauchy_alt":
        A, B, C = np.asarray(coeffs, float)
        
        #Alternative cauchy form: n(λ) = A + B*λ^2 + C/λ^2
        return A + B*L2 + C*(lam_um**-2) 

    if model in ("const", "constant"):
        n0 = float(np.asarray(coeffs).reshape(-1)[0])
        return np.full_like(lam_um, n0) # Constant: n(λ) = c

    raise ValueError(f'Unknown model "{model}"')


def film_refractive_index(TM, Tm, s):
    """
    Utilizes Swanepoel Eq. 11 to calculate refractive index from envelopes and
    substrate index s.
    
    Input:
        TM (array) : Upper envelope
        Tm (array) : Lower envelope
        s (array) : Substrate refrative index at wavelengths matching TM and Tm
    Returns:
        n (array) : n(λ) - refractive index dependent on lambda
    """

    N = 2*s*(TM-Tm)/(TM*Tm)+(s**2+1)/2
    n = np.sqrt(N+np.sqrt(N**2-s**2))

    return n

