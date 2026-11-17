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

import numpy as np

from .levels import gaussian_latitudes

__all__ = ["legendre_norm", "gauss_weights", "analysis", "synthesis", "file_order"]


def gauss_weights(jmax: int) -> np.ndarray:
    """Pesos de Gauss (soma 2), na ordem norte -> sul."""
    _, w = np.polynomial.legendre.leggauss(int(jmax))
    return np.asarray(w[::-1], dtype=np.float64)


def _norm_scale(m: int) -> float:
    """Fator de normalizacao por ordem zonal (1/sqrt2 em m=0, 1/2 caso geral)."""
    return (1.0 / np.sqrt(2.0)) if m == 0 else 0.5


def legendre_columns(mend: int, mu):
    """Gera as colunas `pln[m, m:]` de uma ordem zonal por vez.

    Mesma recorrencia e normalizacao de :func:`legendre_norm`, mas sem
    materializar o cubo `(m+1, n+1, nlat)`: em T666 isso seriam ~3,6 GB
    (667 x 667 x 500 x 8 B) por campo. Aqui o pico e uma coluna
    `(mend+1-m, nlat)`, o que deixa T213/T254/T666 viaveis.

    Yields
    ------
    (m, pln_m)
        `pln_m` tem forma `(mend+1-m, nlat)` e cobre n = m .. mend.
    """
    mu = np.asarray(mu, dtype=np.float64).ravel()
    ny = mu.size
    n_max = int(mend)
    somx2 = np.sqrt(np.maximum(1.0 - mu ** 2, 0.0))

    # A recorrencia crua de P_n^m cresce como O(10^n): em T126 ja anda perto do
    # limite do double e em T213 estoura (NaN). Mantemos a mesma recorrencia,
    # mas reescalonamos a coluna quando ela passa de RESCALE_LIMIT e acumulamos
    # o log do fator em `log_scale`; a normalizacao final (que e minúscula,
    # exp(logn) ~ 1e-300 para n grande) e somada em log, sem overflow.
    RESCALE_LIMIT = 1e150
    pmm = np.ones(ny)                      # P(0,0), ja dividido por exp(log_diag)
    log_diag = 0.0
    for m in range(n_max + 1):
        if m > 0:
            pmm = (2 * m - 1) * somx2 * pmm    # P(m,m) pela diagonal
            peak = float(np.max(np.abs(pmm)))
            if peak > RESCALE_LIMIT:
                pmm /= peak
                log_diag += np.log(peak)
        ns = np.arange(m, n_max + 1)
        col = np.empty((ns.size, ny))
        col[0] = pmm
        log_scale = log_diag
        if ns.size > 1:
            col[1] = mu * (2 * m + 1) * pmm
            for n in range(m + 2, n_max + 1):
                k = n - m
                col[k] = ((2 * n - 1) * mu * col[k - 1] - (n + m - 1) * col[k - 2]) / (n - m)
                peak = float(np.max(np.abs(col[k])))
                if peak > RESCALE_LIMIT:
                    col[:k + 1] /= peak
                    log_scale += np.log(peak)
        logn = 0.5 * (np.log(2.0 if m > 0 else 1.0) + np.log(2 * ns + 1)
                      + _lgamma_vec(ns - m + 1) - _lgamma_vec(ns + m + 1))
        yield m, (np.exp(logn + log_scale)[:, None] * _norm_scale(m)) * col


def _lgamma_vec(values: np.ndarray) -> np.ndarray:
    """lgamma vetorizado (os argumentos sao inteiros pequenos)."""
    from math import lgamma as _lg
    return np.fromiter((_lg(int(v)) for v in values), dtype=np.float64, count=len(values))


def legendre_norm(mend: int, mu) -> np.ndarray:
    """pln[m, n, j] normalizado (n>=m) para as latitudes mu (sin(lat)).

    Materializa o cubo completo (memoria O(mend^2 * nlat)): em T666 sao
    ~3,6 GB. Prefira :func:`legendre_columns` em codigo novo; esta funcao
    existe para diagnostico e testes e usa o mesmo caminho numerico.
    """
    mu = np.asarray(mu, dtype=np.float64).ravel()
    pln = np.zeros((int(mend) + 1, int(mend) + 1, mu.size))
    for m, col in legendre_columns(mend, mu):
        pln[m, m:, :] = col
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
    A = _fourier(f, mend)

    coef = np.zeros((mend + 1, mend + 1), dtype=np.complex128)
    for m, pln_m in legendre_columns(mend, mu):
        cn = A[:half, m]
        cs = A[::-1][:half, m]
        fp = cn + cs
        fm = cn - cs
        ns = np.arange(m, mend + 1)
        parity = (ns - m) % 2 == 0
        s = np.where(parity[:, None], fp[None, :], fm[None, :])      # (n, lat)
        coef[m, m:] = (pln_m * (w * s)).sum(axis=1)
    return coef


def synthesis(coef, mend: int, nlat: int, nlon: int) -> np.ndarray:
    """Reconstroi a grade (nlat, nlon) a partir dos coeficientes (m+1, n+1)."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m, pln_m in legendre_columns(mend, mu):
        trig = np.cos(m * lam) if m else np.ones(nlon)
        sint = np.sin(m * lam) if m else np.zeros(nlon)
        factor = 1.0 if m == 0 else 2.0
        c = coef[m, m:]
        y = pln_m.T @ c                  # (nlat,) complexo
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
