"""Transformada espectral do BAM (grade <-> harmonicos esfericos).

Porte de ICnAGCM/Decomposition.f90 + Legendre.f90 + Fourier.f90:

* FFT: grade -> onda usa kernel e^{+im lambda} com 1/n (np.fft.ifft); a
  reconstrucao usa f = A(0) + 2 Re[sum_{m>0} A(m) e^{-im lambda}] (+ modo de
  Nyquist).  Um cosseno unitario da Re A(1) = 0.5;
* SymAsy: fp = f_norte + f_sul, fm = f_norte - f_sul (sem fator 1/2);
* Legendre normalizado: pln = Pbar/sqrt(2) para m=0 e Pbar/2 para m>0, com
  Pbar o harmonico real 4pi-normalizado; assim pln(0,0)=1/sqrt(2) e
  sum_meia 2*w*pln^2 = 1 (Cte = 1/sqrt(2), logo a(0,0) = sqrt(2)*media);
* coeficiente: c(m,n) = sum_meia w_j pln(m,n,j) * (fp se (n-m) par senao fm);
* ordem do arquivo (TransSpherHarm, ICase=1): m externo, n interno.
"""

from __future__ import annotations

from math import lgamma

import numpy as np

from .levels import gaussian_latitudes

__all__ = ["legendre_norm", "gauss_weights", "analysis", "synthesis", "file_order"]


def gauss_weights(jmax: int) -> np.ndarray:
    """Pesos de Gauss (soma 2), na ordem norte -> sul."""
    _, w = np.polynomial.legendre.leggauss(int(jmax))
    return np.asarray(w[::-1], dtype=np.float64)


def legendre_norm(mend: int, mu) -> np.ndarray:
    """pln[m, n, j] normalizado (n>=m) para as latitudes mu (sin(lat))."""
    mu = np.asarray(mu, dtype=np.float64).ravel()
    ny = mu.size
    m_max = int(mend)
    p = np.zeros((m_max + 1, m_max + 1, ny))
    p[0, 0] = 1.0
    somx2 = np.sqrt(np.maximum(1.0 - mu ** 2, 0.0))
    for m in range(1, m_max + 1):
        p[m, m] = (2 * m - 1) * somx2 * p[m - 1, m - 1]
    for m in range(m_max + 1):
        if m < m_max:
            p[m + 1, m] = mu * (2 * m + 1) * p[m, m]
        for n in range(m + 2, m_max + 1):
            p[n, m] = ((2 * n - 1) * mu * p[n - 1, m] - (n + m - 1) * p[n - 2, m]) / (n - m)

    pln = np.zeros_like(p)
    for m in range(m_max + 1):
        for n in range(m, m_max + 1):
            logn = 0.5 * (np.log(2.0 if m > 0 else 1.0) + np.log(2 * n + 1)
                          + lgamma(n - m + 1) - lgamma(n + m + 1))
            scale = (1.0 / np.sqrt(2.0)) if m == 0 else 0.5
            pln[m, n] = np.exp(logn) * scale * p[n, m]
    return pln


def _fourier(field, mend):
    """Coeficientes de Fourier (complexos) dos primeiros mend+1 modos zonais."""
    f = np.asarray(field, dtype=np.float64)
    if f.ndim != 2:
        raise ValueError("field deve ser (nlat, nlon)")
    return np.fft.ifft(f, axis=1)[:, :mend + 1]


def analysis(field, mend: int) -> np.ndarray:
    """Decompoe field(nlat, nlon) em coeficientes (m+1, n+1) complexos (cos + i sen)."""
    f = np.asarray(field, dtype=np.float64)
    nlat, nlon = f.shape
    if nlat % 2:
        raise ValueError("nlat deve ser par (grade gaussiana)")
    half = nlat // 2
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))[:half]
    w = gauss_weights(nlat)[:half]
    pln = legendre_norm(mend, mu)
    A = _fourier(f, mend)

    coef = np.zeros((mend + 1, mend + 1), dtype=np.complex128)
    for m in range(mend + 1):
        cn = A[:half, m]
        cs = A[::-1][:half, m]
        fp = cn + cs
        fm = cn - cs
        ns = np.arange(m, mend + 1)
        parity = (ns - m) % 2 == 0
        s = np.where(parity[:, None], fp[None, :], fm[None, :])      # (n, lat)
        coef[m, m:] = (pln[m, m:, :] * (w * s)).sum(axis=1)
    return coef


def synthesis(coef, mend: int, nlat: int, nlon: int) -> np.ndarray:
    """Reconstroi a grade (nlat, nlon) a partir dos coeficientes (m+1, n+1)."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    pln = legendre_norm(mend, mu)
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m in range(mend + 1):
        trig = np.cos(m * lam) if m else np.ones(nlon)
        sint = np.sin(m * lam) if m else np.zeros(nlon)
        factor = 1.0 if m == 0 else 2.0
        c = coef[m, m:]
        y = pln[m, m:, :].T @ c          # (nlat,) complexo
        out += factor * (np.real(y)[:, None] * trig[None, :] + np.imag(y)[:, None] * sint[None, :])
    return out


def file_order(coef) -> np.ndarray:
    """Reordena (m+1, n+1) para a ordem do arquivo: m externo, n interno (pares re/im).

    TransSpherHarm com ICase=1: l percorre mm=1..Mend1, nn=1..Mend2-mm e copia
    qCoef(la0(mm,nn)) para a posicao sequencial 2l-1, 2l.
    """
    c = np.asarray(coef)
    mend = c.shape[0] - 1
    out = np.zeros(2 * ((mend + 1) * (mend + 2) // 2))
    pos = 0
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            out[pos] = c[m, n].real
            out[pos + 1] = c[m, n].imag
            pos += 2
    return out


def from_file_order(vector, mend: int) -> np.ndarray:
    v = np.asarray(vector, dtype=np.float64).ravel()
    coef = np.zeros((mend + 1, mend + 1), dtype=np.complex128)
    pos = 0
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            coef[m, n] = complex(v[pos], v[pos + 1])
            pos += 2
    return coef

def hoskins_weights(mend: int, per_cut: float = 0.12, iterations: int = 10) -> np.ndarray:
    """Pesos do filtro de Hoskins usado em SmoothCoef (SpectralICn.f90).

    w(n) = exp(cx*n*(n+1))**Iter, com
    rm = 2M(2M+1), rn = (M-1)M, red = per_cut**(-rm^2/(rn^2*Iter)),
    cx = -log(red)/rm^2.
    """
    M = float(mend)
    rm = 2.0 * M * (2.0 * M + 1.0)
    rn = (M - 1.0) * M
    red = per_cut ** (-(rm * rm) / (rn * rn) / iterations)
    cx = -np.log(red) / (rm * rm)
    n = np.arange(mend + 1, dtype=float)
    rmn = n * (n + 1.0)
    return np.exp(cx * rmn * rmn) ** iterations


def smooth_coefficients(coef, mend: int, per_cut: float = 0.12, iterations: int = 10) -> np.ndarray:
    """Aplica o filtro de Hoskins aos coeficientes (m+1, n+1)."""
    w = hoskins_weights(mend, per_cut, iterations)
    c = np.asarray(coef).copy()
    c *= w[None, :]
    return c
