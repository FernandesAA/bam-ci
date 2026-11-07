"""Formato do arquivo de condicao inicial espectral do BAM.

Especificacao (doc/Reanalysis2CPTECAGCM.pdf e Sources/ICnAGCM/SpectralICn.f90):

    registro 1: ForecastDay(i4), TimeOfDay(r4), DateInp(4xi4), DateCur(4xi4),
                a(K+1) [Pa] e b(K+1)              (caso hibrido)
                ou sigma_interfaces(K+1), sigma_camadas(K)   (caso sigma)
    registro 2: Topo(M coeffs)      [m]
    registro 3: LnPs(M coeffs)      [ln(ps/10), ps em hPa]
    registros 4 .. K+3: Tvir(k)     [K]
    registros K+4 .. 3K+3: Divg(k), Vort(k) alternados por nivel [1/s]
    registros 3K+4 .. 4K+3: SpHu(k) [kg/kg]

Todos os registros sao unformatted sequencial com marcadores de 4 bytes,
big endian, coerentes com fortran_io.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np

from .fortran_io import encode_record, read_records

CTE = 1.0 / np.sqrt(2.0)


def n_coeffs(mend: int) -> int:
    """Numero de coeficientes reais por nivel: (M+1)(M+2)."""
    return (int(mend) + 1) * (int(mend) + 2)


def coeff_bytes(mend: int) -> int:
    """Bytes de um registro de coeficientes (coeffs + 2 marcadores)."""
    return 4 * (n_coeffs(mend) + 2)


def header_words(kmax: int, style: str = "sigma") -> int:
    """Palavras do registro 1 (sem marcadores).

    style="sigma":   ForecastDay, TimeOfDay, 4+4 datas, SigmaI(K+1), SigmaL(K)  -> 11+2K
    style="hybrid":  ForecastDay, TimeOfDay, 4+4 datas, a(K+1), b(K+1)          -> 12+2K
    """
    if style == "sigma":
        return 11 + 2 * int(kmax)
    if style == "hybrid":
        return 12 + 2 * int(kmax)
    raise ValueError("style deve ser 'sigma' ou 'hybrid'")


def header_bytes(kmax: int, style: str = "sigma") -> int:
    return 4 * (header_words(kmax, style) + 2)


def ic_file_size(mend: int, kmax: int, style: str = "sigma") -> int:
    """Tamanho total do arquivo de IC.

    Com style="sigma" reproduz o arquivo de producao medido na jaci
    (GANL...TQ126L042 = 11.055.828 bytes para T126L42).
    """
    return header_bytes(kmax, style) + (2 + 4 * kmax) * coeff_bytes(mend)


def pair_index(mend: int) -> np.ndarray:
    """Matriz (M+1, M+1) com o indice 0-based do par (m, n); -1 fora do triangulo.

    Ordem do arquivo (TransSpherHarm, ICase=1): m externo, n interno,
    p = m(M+1) - m(m-1)/2 + (n-m).  Real em 2p, imaginario em 2p+1.
    """
    m_max = int(mend)
    idx = np.full((m_max + 1, m_max + 1), -1, dtype=np.int64)
    for m in range(m_max + 1):
        base = m * (m_max + 1) - m * (m - 1) // 2
        for n in range(m, m_max + 1):
            idx[m, n] = base + (n - m)
    return idx


def pack_coeffs(mend: int, coeffs: np.ndarray) -> np.ndarray:
    """Converte pares (cos, sen) em vetor achatado na ordem do modelo."""
    return np.asarray(coeffs, dtype=np.float64)


@dataclass
class ICHeader:
    """Cabecalho do arquivo de IC."""

    forecast_day: int = 0
    time_of_day: float = 0.0
    date_inp: tuple = (0, 1, 1, 2000)      # hora, mes, dia, ano
    date_cur: tuple = (0, 1, 1, 2000)
    a_hybr: Optional[np.ndarray] = None    # (K+1) em Pa
    b_hybr: Optional[np.ndarray] = None    # (K+1)
    sigma_inter: Optional[np.ndarray] = None
    sigma_layers: Optional[np.ndarray] = None

    @property
    def hybrid(self) -> bool:
        return self.a_hybr is not None and self.b_hybr is not None

    def encode(self) -> bytes:
        parts = [
            struct.pack(">i", int(self.forecast_day)),
            struct.pack(">f", float(self.time_of_day)),
            struct.pack(">4i", *[int(v) for v in self.date_inp]),
            struct.pack(">4i", *[int(v) for v in self.date_cur]),
        ]
        if self.hybrid:
            parts.append(np.asarray(self.a_hybr, dtype=">f4").tobytes())
            parts.append(np.asarray(self.b_hybr, dtype=">f4").tobytes())
        else:
            if self.sigma_inter is None or self.sigma_layers is None:
                raise ValueError("cabecalho sigma exige sigma_inter e sigma_layers")
            parts.append(np.asarray(self.sigma_inter, dtype=">f4").tobytes())
            parts.append(np.asarray(self.sigma_layers, dtype=">f4").tobytes())
        return b"".join(parts)


def write_ic(path: str | Path, header: ICHeader, topo: np.ndarray, lnps: np.ndarray,
             tvir: np.ndarray, divg: np.ndarray, vort: np.ndarray,
             sphu: np.ndarray) -> Path:
    """Escreve o arquivo de IC no formato do BAM.

    topo, lnps: (ncoeff); tvir, divg, vort, sphu: (ncoeff, kmax).
    """
    topo = np.asarray(topo, dtype=">f4").ravel()
    lnps = np.asarray(lnps, dtype=">f4").ravel()
    tvir = np.asarray(tvir, dtype=">f4")
    divg = np.asarray(divg, dtype=">f4")
    vort = np.asarray(vort, dtype=">f4")
    sphu = np.asarray(sphu, dtype=">f4")
    if tvir.ndim != 2:
        raise ValueError("tvir deve ser (ncoeff, kmax)")
    kmax = tvir.shape[1]
    for name, arr in (("divg", divg), ("vort", vort), ("sphu", sphu)):
        if arr.shape != tvir.shape:
            raise ValueError("%s tem forma %s, esperada %s" % (name, arr.shape, tvir.shape))

    p = Path(path)
    with p.open("wb") as fh:
        fh.write(encode_record(header.encode()))
        fh.write(encode_record(topo.tobytes()))
        fh.write(encode_record(lnps.tobytes()))
        for k in range(kmax):
            fh.write(encode_record(tvir[:, k].tobytes()))
        for k in range(kmax):
            fh.write(encode_record(divg[:, k].tobytes()))
            fh.write(encode_record(vort[:, k].tobytes()))
        for k in range(kmax):
            fh.write(encode_record(sphu[:, k].tobytes()))
    return p


def read_ic(path: str | Path) -> Dict[str, object]:
    """Le o arquivo de IC (equivalente ao CheckICn)."""
    records = list(read_records(path))
    if len(records) < 3:
        raise ValueError("%s: registros insuficientes" % path)
    head = records[0]
    kmax = (len(records) - 3) // 4
    forecast_day, time_of_day = struct.unpack(">if", head[:8])
    date_inp = struct.unpack(">4i", head[8:24])
    date_cur = struct.unpack(">4i", head[24:40])
    rest = np.frombuffer(head[40:], dtype=">f4").astype(np.float64)
    if rest.size == 2 * (kmax + 1):
        a = rest[:kmax + 1]
        b = rest[kmax + 1:]
    elif rest.size == 2 * kmax + 1:
        a = rest[:kmax + 1]
        b = rest[kmax + 1:]
    else:
        raise ValueError("registro 1 com %d palavras incompativel com kmax=%d" % (rest.size, kmax))
    topo = np.frombuffer(records[1], dtype=">f4").astype(np.float64)
    lnps = np.frombuffer(records[2], dtype=">f4").astype(np.float64)
    ncoeff = topo.size
    tvir = np.empty((ncoeff, kmax))
    divg = np.empty((ncoeff, kmax))
    vort = np.empty((ncoeff, kmax))
    sphu = np.empty((ncoeff, kmax))
    pos = 3
    for k in range(kmax):
        tvir[:, k] = np.frombuffer(records[pos], dtype=">f4")
        pos += 1
    for k in range(kmax):
        divg[:, k] = np.frombuffer(records[pos], dtype=">f4")
        vort[:, k] = np.frombuffer(records[pos + 1], dtype=">f4")
        pos += 2
    for k in range(kmax):
        sphu[:, k] = np.frombuffer(records[pos], dtype=">f4")
        pos += 1
    return dict(forecast_day=forecast_day, time_of_day=time_of_day, date_inp=date_inp,
                date_cur=date_cur, a_hybr=a, b_hybr=b, topo=topo, lnps=lnps,
                tvir=tvir, divg=divg, vort=vort, sphu=sphu, kmax=kmax, ncoeff=ncoeff)


def check_ic(path: str | Path, mend: int | None = None) -> List[str]:
    """Checagens objetivas que o CheckICn do Fortran nao faz."""
    data = read_ic(path)
    problems: List[str] = []
    kmax = int(data["kmax"])
    ncoeff = int(data["ncoeff"])
    topo = np.asarray(data["topo"])
    if mend is not None and ncoeff != n_coeffs(mend):
        problems.append("ncoeff=%d, esperado %d para T%d" % (ncoeff, n_coeffs(mend), mend))
    for name in ("topo", "lnps", "tvir", "divg", "vort", "sphu"):
        arr = np.asarray(data[name], dtype=float)
        if not np.all(np.isfinite(arr)):
            problems.append("%s tem valores nao finitos" % name)
        if np.nanmax(np.abs(arr)) > 1e30:
            problems.append("%s com magnitude absurda (%.3g)" % (name, np.nanmax(np.abs(arr))))
    # segundo coeficiente de cada par: seno de m=0 tem de ser nulo
    pairs = ncoeff // 2
    for name in ("topo", "lnps", "tvir", "divg", "vort", "sphu"):
        arr = np.asarray(data[name], dtype=float)
        if arr.ndim == 1:
            arr = arr[:, None]
        problems.extend(_check_zonal(name, arr, mend))
    return problems


def _check_zonal(name: str, arr: np.ndarray, mend: int | None) -> List[str]:
    if mend is None:
        return []
    idx = pair_index(mend)
    m_max = mend
    bad = []
    col = idx[0, 0]
    imag_idx = 2 * col + 1
    if imag_idx < arr.shape[0] and np.max(np.abs(arr[imag_idx])) > 1e-3:
        bad.append("%s: seno de m=0 nao nulo (max %.3g)" % (name, np.max(np.abs(arr[imag_idx]))))
    return bad[:3]
