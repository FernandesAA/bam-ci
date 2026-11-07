"""Leitura de binarios GrADS (.dat) descritos por um .ctl.

Layout classico (dado sequencial, float32, x mais rapido), por tempo:
    para cada t: para cada variavel: para cada nivel: para cada y: para cada x

Com a opcao "template" o .ctl aponta para um arquivo por tempo
(ex.: dset ^analise.%y4%m2%d2.dat). O pacote cobre os dois casos.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from .ctl import Ctl, CtlError, read_ctl


def _expand_template(pattern: str, when, options=frozenset()) -> str:
    """Troca %y4 %y2 %m2 %d2 %h2 %n2 pela data do tempo (1-based)."""
    text = pattern
    if "365_day_calendar" in options:
        # o calendario de 365 dias do GrADS ignora bissexto; a diferenca so aparece
        # em 29/fev, que nao existe nesse calendario.
        pass
    table = {
        "%y4": "%04d" % when.year,
        "%y2": "%02d" % (when.year % 100),
        "%m2": "%02d" % when.month,
        "%d2": "%02d" % when.day,
        "%h2": "%02d" % when.hour,
        "%n2": "%02d" % when.month,
    }
    for key, value in table.items():
        text = text.replace(key, value)
    return text


class GradsFile:
    """Acesso aleatorio a um binario GrADS descrito por um descritor .ctl."""

    def __init__(self, ctl: Ctl | str | Path, *, byteorder: str | None = None,
                 undef_to_nan: bool = True):
        if not isinstance(ctl, Ctl):
            ctl = read_ctl(ctl)
        self.ctl = ctl
        if byteorder is None:
            if "little_endian" in ctl.options and "big_endian" in ctl.options:
                raise CtlError("options declara little_endian e big_endian ao mesmo tempo")
            byteorder = "little" if "little_endian" in ctl.options else "big"
        if byteorder not in ("little", "big"):
            raise ValueError("byteorder deve ser 'little' ou 'big'")
        self.byteorder = byteorder
        self.dtype = np.dtype("<f4" if byteorder == "little" else ">f4")
        self.undef_to_nan = bool(undef_to_nan)
        self.base = ctl.path.parent if ctl.path is not None else Path(".")
        self._layout = self._build_layout()

    # -- layout ------------------------------------------------------------
    def _build_layout(self) -> Dict[str, Tuple[int, int]]:
        """Offset (em bytes) e numero de niveis de cada variavel dentro do bloco de tempo."""
        ctl = self.ctl
        layout: Dict[str, Tuple[int, int]] = {}
        offset = 0
        cell = 4 * ctl.nx * ctl.ny
        for var in ctl.variables:
            nlev = len(var.levels) if var.levels else ctl.nz
            layout[var.name] = (offset, nlev)
            offset += nlev * cell
        return layout

    @property
    def record_bytes(self) -> int:
        return 4 * self.ctl.nx * self.ctl.ny

    @property
    def time_block_bytes(self) -> int:
        return sum(nlev * self.record_bytes for _, nlev in self._layout.values())

    # -- arquivos ----------------------------------------------------------
    def _file_for(self, t: int) -> Path:
        ctl = self.ctl
        if not (1 <= t <= ctl.nt):
            raise IndexError("t=%d fora de 1..%d" % (t, ctl.nt))
        dset = ctl.dset
        if "template" in ctl.options:
            when = list(ctl.tdef.times())[t - 1]
            dset = _expand_template(dset, when, ctl.options)
        if dset.startswith("^"):
            dset = dset[1:].strip()
            return self.base / dset
        return Path(dset)

    def files(self) -> list[Path]:
        return [self._file_for(t) for t in range(1, self.ctl.nt + 1)]

    # -- leitura -----------------------------------------------------------
    def _read_raw(self, path: Path, offset: int, count: int) -> np.ndarray:
        if not path.exists():
            raise FileNotFoundError(
                "binario GrADS nao encontrado: %s (dset=%r)" % (path, self.ctl.dset)
            )
        need = offset + count * 4
        size = path.stat().st_size
        if size < need:
            raise ValueError(
                "arquivo %s tem %d bytes; leitura pede %d (offset %d, %d valores)"
                % (path, size, need, offset, count)
            )
        with path.open("rb") as fh:
            fh.seek(offset)
            data = fh.read(count * 4)
        arr = np.frombuffer(data, dtype=self.dtype, count=count).astype(np.float32, copy=False)
        return np.asarray(arr, dtype=np.float64)

    def _mask(self, values: np.ndarray) -> np.ndarray:
        if not self.undef_to_nan:
            return values
        undef = float(self.ctl.undef)
        out = values.astype(np.float64, copy=True)
        mask = np.abs(out) >= 1.0e30
        if undef != 0.0:
            mask |= np.isclose(out, undef, rtol=1e-6, atol=0.0)
        out[mask] = np.nan
        return out

    def read(self, name: str, t: int = 1, z: int = 1) -> np.ndarray:
        """Le um campo 2D (y, x) da variavel name, tempo t e nivel z (1-based)."""
        ctl = self.ctl
        if name not in self._layout:
            ctl.variable(name)  # levanta KeyError com a lista de variaveis
        st_off, nlev = self._layout[name]
        if not (1 <= z <= nlev):
            raise IndexError("z=%d fora de 1..%d para %s" % (z, nlev, name))
        path = self._file_for(t)
        # com template cada arquivo tem apenas um tempo: offset nao inclui o bloco t
        base = 0 if "template" in ctl.options else (t - 1) * self.time_block_bytes
        offset = base + st_off + (z - 1) * self.record_bytes
        raw = self._read_raw(path, offset, ctl.nx * ctl.ny)
        field = raw.reshape(ctl.ny, ctl.nx)
        if "yrev" in ctl.options:
            field = field[::-1, :]
        return self._mask(field)

    def read_volume(self, name: str, t: int = 1) -> np.ndarray:
        """Le todos os niveis (nz, ny, nx) de uma variavel num tempo."""
        _, nlev = self._layout[name]
        return np.stack([self.read(name, t=t, z=k) for k in range(1, nlev + 1)], axis=0)

    def read_time(self, t: int = 1) -> Dict[str, np.ndarray]:
        """Le todas as variaveis de um tempo: {nome: (nz, ny, nx)}."""
        return {v.name: self.read_volume(v.name, t=t) for v in self.ctl.variables}

    def as_dataset(self, t: int = 1):
        """Devolve xarray.Dataset (requer xarray instalado)."""
        import xarray as xr

        grid = self.ctl.grid()
        data = {}
        for var in self.ctl.variables:
            vol = self.read_volume(var.name, t=t)
            if var.levels:
                grid_z = np.asarray(var.levels, dtype=np.float64)
            else:
                grid_z = grid["level"]
            data[var.name] = (("level", "lat", "lon"), vol, {"units": var.description})
            _ = grid_z
        return xr.Dataset(data, coords={"lat": grid["lat"], "lon": grid["lon"],
                                        "level": grid["level"]})
