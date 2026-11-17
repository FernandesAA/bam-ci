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
from .spectral import analysis, legendre_columns

A_EARTH = 6.37e6


def dpln_columns(mend: int, mu):
    """Gera (m, d pln / d mu) por ordem zonal, sem tabela de escala.

    Usa a recorrencia em termos dos pln ja normalizados:

        (1 - mu^2) dpln_n^m/dmu = -n mu pln_n^m + (n+m) r pln_{n-1}^m

    com `r = s_{n-1}/s_n = sqrt((2n-1)(n+m) / ((2n+1)(n-m)))` o quociente das
    constantes de normalizacao. A forma antiga dividia pela escala
    `exp(logn)` (que underflowa para 0 em n grande) e estourava em T213+.
    """
    mu = np.asarray(mu, dtype=np.float64).ravel()
    one_minus = np.maximum(1.0 - mu ** 2, 1e-12)
    for m, pln_m in legendre_columns(mend, mu):
        ns = np.arange(m, mend + 1, dtype=np.float64)
        out = np.empty_like(pln_m)
        # n == m: P_{n-1}^m = 0
        out[0] = (-ns[0] * mu * pln_m[0]) / one_minus
        for k in range(1, ns.size):
            n = ns[k]
            # s_n / s_{n-1} (inverso do quociente das normalizacoes)
            ratio = np.sqrt(((2 * n + 1) * (n - m)) / ((2 * n - 1) * (n + m)))
            out[k] = (-n * mu * pln_m[k] + (n + m) * ratio * pln_m[k - 1]) / one_minus
        yield m, out


def _dpln(mend: int, mu) -> np.ndarray:
    """d pln / d mu (versao materializada; use :func:`dpln_columns`)."""
    mu = np.asarray(mu, dtype=np.float64).ravel()
    out = np.zeros((int(mend) + 1, int(mend) + 1, mu.size))
    for m, col in dpln_columns(mend, mu):
        out[m, m:, :] = col
    return out


def _accumulate(out, m, c, y, lam, nlon):
    """Soma a contribuicao da ordem zonal m (nucleo comum das sinteses)."""
    trig = np.cos(m * lam) if m else np.ones(nlon)
    sint = np.sin(m * lam) if m else np.zeros(nlon)
    factor = 1.0 if m == 0 else 2.0
    out += factor * (np.real(y)[:, None] * trig[None, :] + np.imag(y)[:, None] * sint[None, :])


def _synth_fourier(coef, mend: int, nlat: int, nlon: int, weight=None) -> np.ndarray:
    """Sintese a partir de coeficientes (m+1, n+1), com peso opcional por (m,n)."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m, pln_m in legendre_columns(mend, mu):
        c = coef[m, m:].copy()
        if weight is not None:
            c = c * weight[m, m:]
        _accumulate(out, m, c, pln_m.T @ c, lam, nlon)
    return out


def _synth_dmu(coef, mend: int, nlat: int, nlon: int) -> np.ndarray:
    """Sintese de d(.)/d mu."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m, dpln_m in dpln_columns(mend, mu):
        c = coef[m, m:]
        _accumulate(out, m, c, dpln_m.T @ c, lam, nlon)
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
