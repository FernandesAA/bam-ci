"""Interpolacao horizontal das analises para a grade gaussiana do modelo.

O legado fazia isso no Fortran (AreaInterpolation.f90 quando a grade de entrada e
mais fina que a de saida, LinearInterpolation.f90 no caso contrario). Aqui comeca
a versao Python: bilinear com longitude periodica, latitudes norte->sul.
"""

from __future__ import annotations

import numpy as np

__all__ = ["interpolate_bilinear", "regular_grid"]


def regular_grid(nlat: int, nlon: int, lat0: float = 90.0, lat1: float = -90.0,
                 lon0: float = 0.0, lon1: float = 360.0) -> tuple:
    """Latitudes/longitudes de uma grade regular (norte->sul, 0..360)."""
    lats = np.linspace(lat0, lat1, nlat)
    lons = np.linspace(lon0, lon1, nlon, endpoint=False)
    return lats, lons


def _axis_weights(src: np.ndarray, dst: np.ndarray) -> tuple:
    """Indices (na ordem original) e pesos lineares para src -> dst."""
    src = np.asarray(src, dtype=np.float64)
    dst = np.asarray(dst, dtype=np.float64)
    order = np.argsort(src)
    ordered = src[order]
    if ordered.size < 2:
        raise ValueError("eixo de origem precisa de ao menos 2 pontos")
    idx = np.clip(np.searchsorted(ordered, dst) - 1, 0, ordered.size - 2)
    w = (dst - ordered[idx]) / (ordered[idx + 1] - ordered[idx])
    return order[idx], order[idx + 1], np.clip(w, 0.0, 1.0)


def interpolate_bilinear(field: np.ndarray, src_lat: np.ndarray, src_lon: np.ndarray,
                         dst_lat: np.ndarray, dst_lon: np.ndarray) -> np.ndarray:
    """Interpola field(ny, nx) de (src_lat, src_lon) para (dst_lat, dst_lon).

    * latitudes em qualquer ordem;
    * longitudes periodicas (0..360), sem pontos duplicados nas pontas;
    * fora do dominio em latitude repete o valor da borda (extrapolacao constante).
    """
    f = np.asarray(field, dtype=np.float64)
    src_lat = np.asarray(src_lat, dtype=np.float64)
    src_lon = np.asarray(src_lon, dtype=np.float64)
    if f.shape != (src_lat.size, src_lon.size):
        raise ValueError("field %s incompativel com lat(%d) x lon(%d)"
                         % (f.shape, src_lat.size, src_lon.size))

    j0, j1, wy = _axis_weights(src_lat, dst_lat)

    lon = np.mod(src_lon, 360.0)
    order = np.argsort(lon)
    ordered = lon[order]
    n = ordered.size
    step = 360.0 / n
    frac = np.mod(np.asarray(dst_lon, dtype=np.float64) - ordered[0], 360.0) / step
    k0 = np.floor(frac).astype(np.int64) % n
    wl = frac - np.floor(frac)
    i0 = order[k0]
    i1 = order[(k0 + 1) % n]

    f0 = f[j0]                                    # (nlat_dst, nlon_src)
    f1 = f[j1]
    vy0 = f0[:, i0] * (1.0 - wl)[None, :] + f0[:, i1] * wl[None, :]
    vy1 = f1[:, i0] * (1.0 - wl)[None, :] + f1[:, i1] * wl[None, :]
    return vy0 * (1.0 - wy)[:, None] + vy1 * wy[:, None]
