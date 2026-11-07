"""Leitor de descritores GrADS (.ctl) do BAM 1.2 subsazonal.

As analises que alimentam o BAM chegam hoje como campos GrADS: um descritor
.ctl e um binario (.dat) com float32 em ordem x-fastest. Este modulo traduz o
descritor para estruturas Python, sem depender do executavel grads.

Formato: https://cola.gmu.edu/grads/gadoc/descriptorfile.html
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterator, Tuple

import numpy as np


class CtlError(ValueError):
    """Descritor GrADS invalido."""


# ---------------------------------------------------------------------------
# eixos
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Axis:
    """Eixo de grade (x, y ou z).

    kind == "linear": values[i] = first + i * step.
    kind == "levels": valores explicitos em levels.
    """

    n: int
    kind: str
    first: float | None = None
    step: float | None = None
    levels: Tuple[float, ...] = ()

    def values(self) -> np.ndarray:
        if self.kind == "linear":
            if self.first is None or self.step is None:
                raise CtlError("eixo linear sem inicio/passo")
            return self.first + self.step * np.arange(self.n, dtype=np.float64)
        return np.asarray(self.levels, dtype=np.float64)

    def __len__(self) -> int:
        return self.n


_STEP_RE = re.compile(r"^(\d+)(mn|hr|dy|mo|yr)$")
_MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
}


def _add_months(base: datetime, months: int) -> datetime:
    month0 = base.month - 1 + months
    year = base.year + month0 // 12
    month = month0 % 12 + 1
    return base.replace(year=year, month=month)


@dataclass(frozen=True)
class TimeAxis:
    """Eixo temporal: inicio + incremento (constante ou mensal/anual)."""

    n: int
    kind: str = "linear"
    start: datetime | None = None
    step: str | None = None

    def times(self) -> Iterator[datetime]:
        if self.start is None:
            raise CtlError("tdef sem data inicial")
        step = self.step or "1dy"
        m = _STEP_RE.match(step.lower())
        if not m:
            raise CtlError("incremento de tempo nao suportado: " + repr(step))
        amount, unit = int(m.group(1)), m.group(2)
        current = self.start
        for _ in range(self.n):
            yield current
            if unit == "mo":
                current = _add_months(current, amount)
            elif unit == "yr":
                current = _add_months(current, 12 * amount)
            else:
                delta = {
                    "mn": timedelta(minutes=amount),
                    "hr": timedelta(hours=amount),
                    "dy": timedelta(days=amount),
                }[unit]
                current = current + delta

    def __len__(self) -> int:
        return self.n


@dataclass(frozen=True)
class Var:
    """Variavel do descritor; levels vazio = usa todos os niveis do zdef."""

    name: str
    levels: Tuple[float, ...] = ()
    description: str = ""


# ---------------------------------------------------------------------------
# descritor
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Ctl:
    """Descritor GrADS completo."""

    dset: str
    title: str
    undef: float
    options: frozenset[str]
    xdef: Axis
    ydef: Axis
    zdef: Axis
    tdef: TimeAxis
    variables: Tuple[Var, ...]
    path: Path | None = None

    # -- conveniencias -----------------------------------------------------
    @property
    def nx(self) -> int:
        return self.xdef.n

    @property
    def ny(self) -> int:
        return self.ydef.n

    @property
    def nz(self) -> int:
        return self.zdef.n

    @property
    def nt(self) -> int:
        return self.tdef.n

    def variable(self, name: str) -> Var:
        for v in self.variables:
            if v.name == name:
                return v
        raise KeyError(
            "variavel " + repr(name) + " ausente; disponiveis: "
            + repr([v.name for v in self.variables])
        )

    def grid(self) -> dict:
        """Coordenadas em dicionario pronto para xarray."""
        lats = self.ydef.values()
        if "yrev" in self.options:
            # o dado bruto vem de sul para norte; o eixo logico e norte->sul.
            lats = np.asarray(lats[::-1], dtype=np.float64)
        return {
            "lon": self.xdef.values(),
            "lat": lats,
            "level": self.zdef.values(),
            "time": np.asarray([np.datetime64(t) for t in self.tdef.times()]),
        }

    def __str__(self) -> str:  # pragma: no cover - so apresentacao
        lines = [
            "dset   = " + self.dset,
            "title  = " + self.title,
            "undef  = " + format(self.undef, "g"),
            "options= " + (" ".join(sorted(self.options)) or "-"),
            "grid   = %d x %d x %d x %d" % (self.nx, self.ny, self.nz, self.nt),
            "vars:",
        ]
        for v in self.variables:
            lev = " ".join(format(x, "g") for x in v.levels) or "(todos)"
            lines.append("  %-10s %-24s %s" % (v.name, lev, v.description))
        return "\n".join(lines)


def _parse_axis(tokens, lineno: int) -> Axis:
    n = int(tokens[1])
    kind = tokens[2].lower()
    if kind == "linear":
        return Axis(n=n, kind="linear", first=float(tokens[3]), step=float(tokens[4]))
    if kind == "levels":
        levels = tuple(float(t) for t in tokens[3:])
        if len(levels) != n:
            raise CtlError(
                "linha %d: zdef declara %d niveis mas lista %d" % (lineno, n, len(levels))
            )
        return Axis(n=n, kind="levels", levels=levels)
    raise CtlError("linha %d: eixo %r nao suportado" % (lineno, kind))


def _parse_time(text: str) -> datetime:
    m = re.match(r"^(\d{1,2})([zZ])(\d{2})([a-zA-Z]{3})(\d{4})$", text)
    if not m:
        raise CtlError("data GrADS invalida: " + repr(text))
    hour, _, day, mon, year = m.groups()
    month = _MONTHS.get(mon.lower())
    if month is None:
        raise CtlError("mes invalido em " + repr(text))
    return datetime(int(year), month, int(day), int(hour))


def _parse_tdef(tokens, lineno: int) -> TimeAxis:
    n = int(tokens[1])
    kind = tokens[2].lower()
    if kind != "linear":
        raise CtlError("linha %d: tdef %r nao suportado (use linear)" % (lineno, kind))
    return TimeAxis(n=n, kind="linear", start=_parse_time(tokens[3]), step=tokens[4].lower())


def parse_ctl(text: str, path: Path | None = None) -> Ctl:
    """Interpreta o texto de um .ctl e devolve Ctl."""
    dset = title = ""
    undef = -9.99e8
    options: set[str] = set()
    xdef = ydef = zdef = tdef = None
    variables: list[Var] = []
    nvars_expected = None

    lines = text.splitlines()
    in_vars = False
    for lineno, raw in enumerate(lines, start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        tokens = line.split()
        key = tokens[0].lower()
        if in_vars:
            if key == "endvars":
                in_vars = False
                continue
            if len(tokens) < 2:
                raise CtlError("linha %d: variavel sem nome" % lineno)
            name = tokens[0]
            levels: Tuple[float, ...] = ()
            description = ""
            # formato GrADS: nome <n_niveis> <nivel1> ... <unidades> <descricao>
            if len(tokens) >= 2:
                try:
                    count = int(tokens[1])
                except ValueError:
                    raise CtlError("linha %d: numero de niveis invalido em %r" % (lineno, line))
                if count < 0 or len(tokens) < 2 + count:
                    raise CtlError("linha %d: niveis insuficientes em %r" % (lineno, line))
                levels = tuple(float(t) for t in tokens[2:2 + count])
                rest = tokens[2 + count:]
                units = rest[0] if rest else ""
                description = " ".join(rest[1:])
                if units:
                    description = (units + " " + description).strip()
            variables.append(Var(name=name, levels=levels, description=description))
            continue
        if key == "dset":
            dset = line[len("dset"):].strip()
        elif key == "title":
            title = line[len("title"):].strip()
        elif key == "undef":
            undef = float(tokens[1])
        elif key == "options":
            options.update(t.lower() for t in tokens[1:])
        elif key == "xdef":
            xdef = _parse_axis(tokens, lineno)
        elif key == "ydef":
            ydef = _parse_axis(tokens, lineno)
        elif key == "zdef":
            zdef = _parse_axis(tokens, lineno)
        elif key == "tdef":
            tdef = _parse_tdef(tokens, lineno)
        elif key == "vars":
            nvars_expected = int(tokens[1])
            in_vars = True
        elif key == "cdate":
            continue
        else:
            raise CtlError("linha %d: chave desconhecida %r" % (lineno, key))

    missing = [n for n, v in (("xdef", xdef), ("ydef", ydef), ("zdef", zdef), ("tdef", tdef)) if v is None]
    if missing:
        raise CtlError("descritor sem " + ", ".join(missing))
    if nvars_expected is not None and len(variables) != nvars_expected:
        raise CtlError(
            "vars declara %d variaveis, encontradas %d" % (nvars_expected, len(variables))
        )
    if not dset:
        raise CtlError("descritor sem dset")
    return Ctl(
        dset=dset, title=title, undef=undef, options=frozenset(options),
        xdef=xdef, ydef=ydef, zdef=zdef, tdef=tdef, variables=tuple(variables), path=path,
    )


def read_ctl(path: str | Path) -> Ctl:
    p = Path(path)
    return parse_ctl(p.read_text(encoding="utf-8", errors="replace"), path=p)
