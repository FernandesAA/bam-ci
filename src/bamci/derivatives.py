"""Divergencia e vorticidade por derivadas espectrais.

Alternativa fiel ao legado (GetConv/GetVort): em vez de diferencas finitas na
grade, as derivadas de u e v sao calculadas no espaco espectral:

    div = 1/(a cos fi) [ d u/d lambda + cos fi d v/d fi - v sen fi ]
    vort = 1/(a cos fi) [ d v/d lambda - cos fi d u/d fi + u sen fi ]

com d/d lambda = i m (nao muda os harmonicos) e d/d fi = cos fi d/d mu obtido
pela recorrencia de Legendre. Validado com rotacao solida (u = u0 cos fi),
que deve dar vort = 2 u0 sen fi / a e div = 0.
"""

from __future__ import annotations

import numpy as np

from .levels import gaussian_latitudes
from .spectral import analysis, gauss_weights, legendre_norm

A_EARTH = 6.37e6


def _scale_table(mend: int) -> np.ndarray:
    from math import lgamma
    s = np.zeros((mend + 1, mend + 1))
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            logn = 0.5 * (np.log(2.0 if m > 0 else 1.0) + np.log(2 * n + 1)
                          + lgamma(n - m + 1) - lgamma(n + m + 1))
            s[m, n] = np.exp(logn) * ((1.0 / np.sqrt(2.0)) if m == 0 else 0.5)
    return s


def _dpln(mend: int, mu) -> np.ndarray:
    """d pln / d mu."""
    mu = np.asarray(mu, dtype=np.float64).ravel()
    pln = legendre_norm(mend, mu)
    scale = _scale_table(mend)
    out = np.zeros_like(pln)
    one_minus = np.maximum(1.0 - mu ** 2, 1e-12)
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            p = pln[m, n] / scale[m, n]
            # (1-mu^2) dP_n^m/dmu = (n+m) P_{n-1}^m - n mu P_n^m
            term = -n * mu * p
            if n > m:      # P_{n-1}^m = 0 quando n == m
                term = term + (n + m) * (pln[m, n - 1] / scale[m, n - 1])
            out[m, n] = scale[m, n] * term / one_minus
    return out


def _synth_fourier(coef, mend: int, nlat: int, nlon: int, weight=None) -> np.ndarray:
    """Sintese a partir de coeficientes (m+1, n+1), com peso opcional por (m,n)."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    pln = legendre_norm(mend, mu)
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m in range(mend + 1):
        c = coef[m, m:].copy()
        if weight is not None:
            c = c * weight[m, m:]
        trig = np.cos(m * lam) if m else np.ones(nlon)
        sint = np.sin(m * lam) if m else np.zeros(nlon)
        factor = 1.0 if m == 0 else 2.0
        y = pln[m, m:, :].T @ c
        out += factor * (np.real(y)[:, None] * trig[None, :] + np.imag(y)[:, None] * sint[None, :])
    return out


def _synth_dmu(coef, mend: int, nlat: int, nlon: int) -> np.ndarray:
    """Sintese de d(.)/d mu."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    dpln = _dpln(mend, mu)
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m in range(mend + 1):
        c = coef[m, m:]
        trig = np.cos(m * lam) if m else np.ones(nlon)
        sint = np.sin(m * lam) if m else np.zeros(nlon)
        factor = 1.0 if m == 0 else 2.0
        y = dpln[m, m:, :].T @ c
        out += factor * (np.real(y)[:, None] * trig[None, :] + np.imag(y)[:, None] * sint[None, :])
    return out


def spectral_div_vort(u, v, nlat: int, nlon: int, mend: int) -> tuple:
    """Divergencia e vorticidade (1/s) por derivadas espectrais."""
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    lat = np.radians(gaussian_latitudes(nlat))
    coslat = np.cos(lat)[:, None]
    sinlat = np.sin(lat)[:, None]

    cu = analysis(u, mend)
    cv = analysis(v, mend)

    # du/dlambda = i m u  (no espaco de coeficientes) -> sintese direta
    lam_w = np.zeros((mend + 1, mend + 1), dtype=complex)
    for m in range(mend + 1):
        lam_w[m, :] = 1j * m
    du_dlam = _synth_fourier(cu, mend, nlat, nlon, weight=lam_w)
    dv_dlam = _synth_fourier(cv, mend, nlat, nlon, weight=lam_w)

    du_dmu = _synth_dmu(cu, mend, nlat, nlon)
    dv_dmu = _synth_dmu(cv, mend, nlat, nlon)

    du_df = coslat * du_dmu
    dv_df = coslat * dv_dmu

    div = (du_dlam + coslat * dv_df - sinlat * v) / (A_EARTH * coslat)
    vort = (dv_dlam - coslat * du_df + sinlat * u) / (A_EARTH * coslat)
    return div, vort
