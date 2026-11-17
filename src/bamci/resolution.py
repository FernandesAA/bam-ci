"""Resolucao do modelo: truncacao espectral (mend) + numero de camadas (kmax).

Centraliza o que antes era fixo em `126`/`42`:

* `parse_resolution("TQ0213L64")` -> `("Q", 213, 64)`;
* `format_resolution("Q", 213, 64)` -> `"TQ0213L64"`;
* `grid_shape(mend)` -> grade gaussiana (ImaxOut, JmaxOut) — reexportado de
  :mod:`bamci.levels` para que a resolucao seja um unico ponto de verdade.

Formatos aceitos (case-insensitive, espacos ignorados)::

    TQ0213L64   TQ213L64   T213L64   Q213L64   213L64
    TQ0213      T213       Q213      213        (kmax omitido)

O `T` inicial e o `Q`/`L` (quadratica/linear) sao opcionais; sem a letra de
grade vale o default do :class:`bamci.config.Config` (`"Q"`).
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

from .levels import grid_shape, fft_friendly_size

__all__ = [
    "GRIDS",
    "parse_resolution",
    "format_resolution",
    "resolve",
    "describe",
    "grid_shape",
    "fft_friendly_size",
]

#: grades suportadas: Q = quadratica (TS), L = linear (TL)
GRIDS = ("Q", "L")

_PATTERN = re.compile(
    r"""
    ^\s*
    [Tt]?                        # prefixo do modelo, opcional
    (?P<grid>[QLql])?            # grade, opcional
    (?P<mend>\d{1,5})            # truncacao
    (?:[Ll](?P<kmax>\d{1,4}))?    # camadas, opcional
    \s*$
    """,
    re.VERBOSE,
)


def parse_resolution(text: str) -> Tuple[str, int, Optional[int]]:
    """Interpreta `"TQ0213L64"` -> `("Q", 213, 64)`.

    Levanta `ValueError` com a lista de formatos aceitos quando o texto nao
    casa. `kmax` volta `None` quando o texto nao traz o `L`.
    """
    if not isinstance(text, str):
        raise ValueError("resolucao deve ser texto (ex.: TQ0213L64), recebido %r" % (text,))
    match = _PATTERN.match(text)
    if match is None:
        raise ValueError(
            "resolucao invalida: %r (formatos: TQ0213L64, T213L42, Q213, 213L64)" % (text,)
        )
    grid = (match.group("grid") or "").upper()
    mend = int(match.group("mend"))
    kmax = match.group("kmax")
    kmax = int(kmax) if kmax is not None else None
    if mend <= 0:
        raise ValueError("truncacao deve ser positiva em %r" % (text,))
    if kmax is not None and kmax <= 0:
        raise ValueError("numero de camadas deve ser positivo em %r" % (text,))
    return grid, mend, kmax


def format_resolution(grid: str, mend: int, kmax: Optional[int] = None) -> str:
    """Monta o rotulo canonico, ex.: `("Q", 213, 64)` -> `"TQ0213L64"`."""
    text = "T%s%04d" % (grid.upper(), int(mend))
    if kmax is not None:
        text += "L%03d" % int(kmax)
    return text


def resolve(
    *,
    resolution: Optional[str] = None,
    grid: Optional[str] = None,
    mend: Optional[int] = None,
    kmax: Optional[int] = None,
    default_grid: str = "Q",
) -> Tuple[str, Optional[int], Optional[int]]:
    """Combina `-r TQ0213L64` com `--mend/--kmax/--grid` avulsos.

    Os argumentos explicitos vencem a resolucao textual (permite
    `-r TQ0213 --kmax 91`). Devolve `(grid, mend, kmax)`; `mend` so volta
    `None` se nem a resolucao nem `--mend` o informarem.
    """
    grid_r = mend_r = kmax_r = None
    if resolution:
        grid_r, mend_r, kmax_r = parse_resolution(resolution)

    out_grid = (grid or grid_r or default_grid).upper()
    if out_grid not in GRIDS:
        raise ValueError("grade invalida: %r (use %s)" % (out_grid, " ou ".join(GRIDS)))

    out_mend = mend if mend is not None else mend_r
    out_kmax = kmax if kmax is not None else kmax_r
    if out_mend is not None and out_mend <= 0:
        raise ValueError("--mend deve ser positivo")
    if out_kmax is not None and out_kmax <= 0:
        raise ValueError("--kmax deve ser positivo")
    return out_grid, out_mend, out_kmax


def describe(grid: str, mend: int, kmax: Optional[int] = None) -> dict:
    """Resumo legivel de uma resolucao (usado por `bamci res`)."""
    from .ic_format import ic_file_size

    nlon, nlat = grid_shape(mend)
    info = {
        "label": format_resolution(grid, mend, kmax),
        "grid": grid.upper(),
        "mend": int(mend),
        "nlon": nlon,
        "nlat": nlat,
        "ncoeff": (int(mend) + 1) * (int(mend) + 2),
    }
    if kmax is not None:
        info["kmax"] = int(kmax)
        info["bytes"] = ic_file_size(int(mend), int(kmax))
    return info
