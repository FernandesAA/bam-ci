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


def _fourier_to_grid(y, nlon):
    """(..., m+1) coeficientes complexos por latitude -> (..., nlon) na grade.

    Mesma convencao de bamci.spectral.synthesis:

        f = A(0) + sum_{m>0} 2 [Re(Y_m) cos m lam + Im(Y_m) sin m lam]

    A soma em m vira dois produtos matriz (nmax x nlon) em vez de acumular m a
    m: mesma matematica, ~10x mais rapido com 42 niveis.
    """
    y = np.asarray(y, dtype=np.complex128)
    nmax = y.shape[-1]
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    m = np.arange(nmax, dtype=np.float64)
    factor = np.where(m > 0, 2.0, 1.0)[:, None]
    cos_m = factor * np.cos(m[:, None] * lam[None, :])
    sin_m = factor * np.sin(m[:, None] * lam[None, :])
    flat = y.reshape(-1, nmax)
    out = flat.real @ cos_m + flat.imag @ sin_m
    return out.reshape(y.shape[:-1] + (nlon,))


def _coef_slice(coef, m, weight=None):
    """Fatia (..., m:, m:) dos coeficientes e aplica o peso por (m,n)."""
    c = coef[..., m, m:].copy()
    if weight is not None:
        c = c * weight[m, m:]
    return c


def _synth_fourier(coef, mend: int, nlat: int, nlon: int, weight=None,
                   columns=None) -> np.ndarray:
    """Sintese a partir de coeficientes (m+1, n+1) ou (nlev, m+1, n+1).

    Com uma dimensao extra de niveis a sintese e' vetorizada (mesmo resultado,
    uma unica passada de legendre_columns). 'columns' permite reaproveitar
    (m, pln_m) ja calculados por legendre_columns (ex.: varios campos na mesma
    grade) - nao e' modificado.
    """
    coef = np.asarray(coef, dtype=np.complex128)
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    y = np.zeros(coef.shape[:-2] + (nlat, mend + 1), dtype=np.complex128)
    for m, pln_m in (columns if columns is not None else legendre_columns(mend, mu)):
        c = _coef_slice(coef, m, weight)
        y[..., :, m] = c @ pln_m
    return _fourier_to_grid(y, nlon)


def _synth_dmu(coef, mend: int, nlat: int, nlon: int, columns=None) -> np.ndarray:
    """Sintese de d(.)/d mu (aceita uma dimensao extra de niveis)."""
    coef = np.asarray(coef, dtype=np.complex128)
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    y = np.zeros(coef.shape[:-2] + (nlat, mend + 1), dtype=np.complex128)
    for m, dpln_m in (columns if columns is not None else dpln_columns(mend, mu)):
        c = _coef_slice(coef, m)
        y[..., :, m] = c @ dpln_m
    return _fourier_to_grid(y, nlon)


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


def wind_from_div_vort(div_sh, vort_sh, mend: int, nlat: int, nlon: int,
                       *, columns=None, dcolumns=None) -> tuple:
    """u, v em grade a partir dos coeficientes espectrais de div/vort.

    Inverso de spectral_div_vort na mesma base: resolve o Laplaciano
    (psi = -a^2 zeta / n(n+1), chi = -a^2 delta / n(n+1)) e usa

        u = -(1/a) dpsi/dphi + (1/(a cos fi)) dchi/dlam
        v =  (1/(a cos fi)) dpsi/dlam + (1/a) dchi/dphi

    O modo (0,0) nao tem vento associado e fica zerado. Usado por
    bamci.mb09 (vento do IC/historico e' sempre reconstruido de div/vort) e
    por bamci.bamhist (escrita do .ctl a partir do .fct).
    """
    div_sh = np.asarray(div_sh, dtype=np.complex128)
    vort_sh = np.asarray(vort_sh, dtype=np.complex128)
    n = np.arange(mend + 1, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(n > 0, -A_EARTH ** 2 / (n * (n + 1.0)), 0.0)
    fac = fac.reshape((1,) * (div_sh.ndim - 1) + (mend + 1,))
    psi = vort_sh * fac
    chi = div_sh * fac

    lat = np.radians(gaussian_latitudes(nlat))
    coslat = np.cos(lat)[:, None]
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    wlam = np.zeros((mend + 1, mend + 1), dtype=np.complex128)
    for m in range(mend + 1):
        wlam[m, :] = 1j * m

    dpsi_dlam = _synth_fourier(psi, mend, nlat, nlon, weight=wlam, columns=columns)
    dchi_dlam = _synth_fourier(chi, mend, nlat, nlon, weight=wlam, columns=columns)
    dpsi_dmu = _synth_dmu(psi, mend, nlat, nlon, columns=dcolumns)
    dchi_dmu = _synth_dmu(chi, mend, nlat, nlon, columns=dcolumns)

    u = -(coslat / A_EARTH) * dpsi_dmu + dchi_dlam / (A_EARTH * coslat)
    v = dpsi_dlam / (A_EARTH * coslat) + (coslat / A_EARTH) * dchi_dmu
    return u, v
