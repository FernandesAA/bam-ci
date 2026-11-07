"""Comparacao entre dois arquivos de condicao inicial.

Uso tipico: comparar a CI gerada pelo bam-ci com a CI de producao do BAM.
Como os campos sao espectrais, a comparacao e feita em duas frentes:

* nos coeficientes (diferenca absoluta e relativa por campo);
* na grade, sintetizando o campo fisico (topografia em m, pressao de superficie
  em hPa, temperatura virtual em K, umidade em kg/kg) e medindo o RMS da
  diferenca em unidades fisicas.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np

from .ic_format import read_ic
from .levels import grid_shape
from .spectral import from_file_order, synthesis

FIELDS = ("topo", "lnps", "tvir", "divg", "vort", "sphu")


def _to_physical(name: str, coef) -> np.ndarray:
    """Campo fisico na grade (usa o primeiro nivel para campos 3D)."""
    arr = np.asarray(coef)
    if arr.ndim > 1:
        arr = arr[:, 0]
    return arr


def compare_ic(path_a: str | Path, path_b: str | Path, mend: int = 126,
               nlat: int | None = None, nlon: int | None = None,
               synth: bool = True) -> Dict[str, object]:
    """Compara duas CIs e devolve um dicionario com as metricas."""
    a = read_ic(path_a)
    b = read_ic(path_b)
    out: Dict[str, object] = {"kmax": (a["kmax"], b["kmax"]),
                              "ncoeff": (a["ncoeff"], b["ncoeff"]),
                              "date_inp": (a["date_inp"], b["date_inp"])}
    same_grid = a["ncoeff"] == b["ncoeff"] and a["kmax"] == b["kmax"]
    if not same_grid:
        out["compatible"] = False
        return out
    out["compatible"] = True
    out["sigma_a"] = np.asarray(a["a_hybr"])
    out["sigma_b"] = np.asarray(b["a_hybr"])
    out["sigma_maxdiff"] = float(np.max(np.abs(np.asarray(a["a_hybr"]) - np.asarray(b["a_hybr"]))))

    if nlat is None or nlon is None:
        nlon, nlat = grid_shape(mend)
    metrics = {}
    for name in FIELDS:
        va = _to_physical(name, a[name])
        vb = _to_physical(name, b[name])
        diff = va - vb
        base_rms = float(np.sqrt(np.mean(va ** 2))) or 1.0
        item = {"coef_rms": float(np.sqrt(np.mean(diff ** 2))),
                "coef_max": float(np.max(np.abs(diff))),
                "coef_rel": float(np.sqrt(np.mean(diff ** 2)) / base_rms)}
        if synth:
            ga = synthesis(from_file_order(va, mend), mend, nlat, nlon)
            gb = synthesis(from_file_order(vb, mend), mend, nlat, nlon)
            d = ga - gb
            if name == "lnps":
                ga = 10.0 * np.exp(ga)
                gb = 10.0 * np.exp(gb)
                d = ga - gb
            item.update({"grid_rms": float(np.sqrt(np.mean(d ** 2))),
                         "grid_max": float(np.max(np.abs(d))),
                         "field_rms": float(np.sqrt(np.mean(ga ** 2)))})
        metrics[name] = item
    out["fields"] = metrics
    return out


def format_report(res: Dict[str, object]) -> str:
    lines = []
    if not res.get("compatible", False):
        return "arquivos incompativeis: %s" % {k: res[k] for k in ("kmax", "ncoeff")}
    lines.append("kmax=%s ncoeff=%s data=%s" % (res["kmax"], res["ncoeff"], res["date_inp"]))
    lines.append("sigma: diferenca maxima = %.3e" % res["sigma_maxdiff"])
    lines.append("%-5s %12s %12s %10s %12s %12s" % ("campo", "rms_coef", "max_coef", "rel", "rms_grade", "max_grade"))
    for name, m in res["fields"].items():
        lines.append("%-5s %12.4e %12.4e %10.3e %12.4e %12.4e" % (
            name, m["coef_rms"], m["coef_max"], m["coef_rel"],
            m.get("grid_rms", float("nan")), m.get("grid_max", float("nan"))))
    return "\n".join(lines)
