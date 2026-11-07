"""Grades e niveis verticais do BAM.

* grade gaussiana: latitudes de Gauss e resolucao FFT (equivalente a GetImaxJmax);
* HybridLevels.Lxxx: coeficientes (a, b) em Pa nas interfaces, p = a + b*ps;
* DeltaSigma.Lxxx: espessuras das camadas em coordenada sigma;
* PressureLevels.*: niveis de pressao das analises (hPa).
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np


def gaussian_latitudes(jmax: int) -> np.ndarray:
    """Latitudes de Gauss (graus, norte -> sul) para jmax linhas."""
    nodes, _ = np.polynomial.legendre.leggauss(int(jmax))
    lat = np.degrees(np.arcsin(nodes))[::-1]   # norte -> sul
    return np.asarray(lat, dtype=np.float64)


def fft_friendly_size(n: int, limit: int = 40000) -> int:
    """Menor inteiro 2^a*3^b*5^c >= n (mesma tabela do GetImaxJmax)."""
    best = None
    value = 1
    while value <= limit:
        v3 = value
        while v3 <= limit:
            v5 = v3
            while v5 <= limit:
                if v5 >= n and (best is None or v5 < best):
                    best = v5
                v5 *= 5
            v3 *= 3
        value *= 2
    if best is None:
        raise ValueError("sem tamanho FFT amigavel >= %d" % n)
    return best


def grid_shape(mend: int) -> Tuple[int, int]:
    """(ImaxOut, JmaxOut) do modelo para truncacao mend (GetImaxJmax)."""
    imax = fft_friendly_size(3 * mend + 1)
    jmax = imax // 2
    if jmax % 2:
        jmax += 1
    return imax, jmax


def read_hybrid(path: str | Path) -> Tuple[np.ndarray, np.ndarray]:
    """Le HybridLevels.Lxxx: (a, b) nas K+1 interfaces, a convertido de Pa para hPa."""
    a, b = [], []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        a.append(float(parts[0]))
        b.append(float(parts[1]))
    a_arr = np.asarray(a, dtype=np.float64) / 100.0   # Pa -> hPa (como o Fortran)
    b_arr = np.asarray(b, dtype=np.float64)
    return a_arr, b_arr


def read_pressure_levels(path: str | Path) -> np.ndarray:
    """Le PressureLevels.* no formato (8F7.1), em hPa."""
    text = Path(path).read_text()
    values = [float(t) for t in text.split()]
    return np.asarray(values, dtype=np.float64)


def read_delta_sigma(path: str | Path) -> np.ndarray:
    """Le DeltaSigma.Lxxx no formato (5F9.6) e normaliza a soma para 1."""
    text = Path(path).read_text()
    values = np.asarray([float(t) for t in text.split()], dtype=np.float64)
    total = values.sum()
    if not np.isclose(total, 1.0, atol=1e-5):
        values = values / total
    return values


def sigma_interfaces(delta_sigma: np.ndarray) -> np.ndarray:
    """Interfaces sigma a partir das espessuras: 0, d1, d1+d2, ..., 1."""
    return np.concatenate([[0.0], np.cumsum(np.asarray(delta_sigma, dtype=np.float64))])


def sigma_layers(delta_sigma: np.ndarray) -> np.ndarray:
    """Centros das camadas sigma (media das interfaces)."""
    inter = sigma_interfaces(delta_sigma)
    return 0.5 * (inter[:-1] + inter[1:])
