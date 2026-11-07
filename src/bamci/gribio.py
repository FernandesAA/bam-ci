"""Leitura das analises em GRIB (ERA5, FNL, GFS) e montagem dos campos.

O legado fazia isso com GrADS: um .gs por fonte, escrevendo um .dat float32
big endian com topo, psnm, T, u, v, q. Aqui o GRIB e lido direto (eccodes) e os
campos ja saem normalizados nas unidades que o Fortran esperava.

Fontes suportadas:
  era5 : shortName ECMWF (t, q, u, v, z, sp, msl, sst, ci)
  nl   : shortName NCEP (t, r, u, v, hgt, pres, prmsl, sst, icec)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence

import numpy as np

EPS = 0.622
EO = 6.1078      # hPa
TO = 273.16      # K
A_TETENS = 17.2693882
B_TETENS = 35.86
GRAV = 9.80665


@dataclass
class RawField:
    """Uma mensagem GRIB ja decodificada."""

    short_name: str
    level_type: str
    level: float
    values: np.ndarray            # (ny, nx) na ordem do arquivo
    latitudes: np.ndarray
    longitudes: np.ndarray
    units: str = ""
    metadata: Dict[str, object] = field(default_factory=dict)

    def __post_init__(self):
        self.values = np.asarray(self.values, dtype=np.float64)
        self.latitudes = np.asarray(self.latitudes, dtype=np.float64)
        self.longitudes = np.asarray(self.longitudes, dtype=np.float64)


# ---------------------------------------------------------------------------
# selecao por fonte
# ---------------------------------------------------------------------------
SOURCE_MAP: Dict[str, Dict[str, List[str]]] = {
    # campos 3D em niveis de pressao
    "era5": {"t": ["t"], "q": ["q"], "u": ["u"], "v": ["v"], "topo": ["z"],
             "sp": ["sp"], "msl": ["msl"], "sst": ["sst"], "ci": ["ci"]},
    "nl": {"t": ["t"], "rh": ["r"], "u": ["u"], "v": ["v"], "topo": ["orog", "hgt"],
           "sp": ["sp", "pres"], "msl": ["prmsl"], "sst": ["t"], "ci": ["ci", "icec"]},
}


def _matches(fld: RawField, names: Sequence[str], level_type: Optional[str] = None) -> bool:
    if fld.short_name not in names:
        return False
    if level_type is not None and fld.level_type != level_type:
        return False
    return True


def select_levels(fields: Iterable[RawField], names: Sequence[str], level_type: str) -> Dict[float, np.ndarray]:
    out: Dict[float, np.ndarray] = {}
    for fld in fields:
        if _matches(fld, names, level_type):
            out[float(fld.level)] = fld.values
    return out


def select_surface(fields: Iterable[RawField], names: Sequence[str],
                   level_types=("surface",)) -> Optional[np.ndarray]:
    for fld in fields:
        if fld.short_name in names and fld.level_type in level_types:
            return fld.values
    return None


def specific_humidity_from_rh(rh: np.ndarray, temperature: np.ndarray, pressure_hpa: float) -> np.ndarray:
    """Tetens + conversao UR -> umidade especifica (identico ao .gs do FNL)."""
    es = EO * np.exp(A_TETENS * (temperature - TO) / (temperature - B_TETENS))
    e = (np.asarray(rh, dtype=np.float64) / 100.0) * es
    return EPS * e / (pressure_hpa - (1.0 - EPS) * e)


def to_hpa(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=np.float64) / 100.0


def geopotential_to_height(z: np.ndarray) -> np.ndarray:
    return np.asarray(z, dtype=np.float64) / GRAV


# ---------------------------------------------------------------------------
# leitura com eccodes
# ---------------------------------------------------------------------------
def _require_eccodes():
    try:
        import eccodes  # noqa: F401
    except ImportError as exc:  # pragma: no cover - depende do ambiente
        raise RuntimeError(
            "eccodes nao instalado: use 'module load eccodes' e pip install eccodes, "
            "ou converta o GRIB com wgrib2"
        ) from exc
    return eccodes


def read_grib(path: str | Path, short_names: Optional[Sequence[str]] = None) -> List[RawField]:
    """Le todas as mensagens de um arquivo GRIB (filtra por shortName, se dado)."""
    ecc = _require_eccodes()
    fields: List[RawField] = []
    with open(path, "rb") as fh:
        while True:
            gid = ecc.codes_grib_new_from_file(fh)
            if gid is None:
                break
            try:
                name = ecc.codes_get(gid, "shortName")
                if short_names is not None and name not in short_names:
                    continue
                values = ecc.codes_get_values(gid)
                ny = ecc.codes_get(gid, "Nj")
                nx = ecc.codes_get(gid, "Ni")
                lats = np.asarray(ecc.codes_get_array(gid, "latitudes"), dtype=float)
                lons = np.asarray(ecc.codes_get_array(gid, "longitudes"), dtype=float)
                if lats.size == ny * nx:
                    lats = lats.reshape(ny, nx)[:, 0]
                if lons.size == ny * nx:
                    lons = lons.reshape(ny, nx)[0, :]
                try:
                    level = float(ecc.codes_get(gid, "level"))
                except Exception:
                    level = 0.0
                try:
                    units = ecc.codes_get(gid, "units")
                except Exception:
                    units = ""
                fields.append(RawField(
                    short_name=name,
                    level_type=ecc.codes_get(gid, "typeOfLevel"),
                    level=level,
                    values=np.asarray(values).reshape(ny, nx),
                    latitudes=lats,
                    longitudes=lons,
                    units=units,
                ))
            finally:
                ecc.codes_release(gid)
    return fields


@dataclass
class Analysis:
    """Campos de uma analise na grade de origem (antes da interpolacao)."""

    pressure_hpa: np.ndarray                    # (K,) crescente? ordem do arquivo
    temperature: np.ndarray                     # (K, ny, nx)
    humidity: np.ndarray                        # (Kq, ny, nx)
    u_wind: np.ndarray                          # (K, ny, nx)
    v_wind: np.ndarray                          # (K, ny, nx)
    surface_pressure: np.ndarray                # (ny, nx) hPa
    sea_level_pressure: Optional[np.ndarray]    # (ny, nx) hPa
    topography: Optional[np.ndarray]            # (ny, nx) m
    sst: Optional[np.ndarray] = None            # (ny, nx) K
    ice: Optional[np.ndarray] = None            # (ny, nx) 0..1
    latitudes: Optional[np.ndarray] = None
    longitudes: Optional[np.ndarray] = None


def _axis(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values


def build_analysis(fields: Sequence[RawField], source: str, levels_hpa: Sequence[float],
                   level_type: str = "isobaricInhPa") -> Analysis:
    """Monta a Analysis a partir das mensagens GRIB de uma fonte."""
    if source not in SOURCE_MAP:
        raise ValueError("fonte desconhecida: %r" % source)
    spec = SOURCE_MAP[source]
    levels = [float(p) for p in levels_hpa]

    def by_level(key: str) -> Dict[float, np.ndarray]:
        return select_levels(fields, spec[key], level_type)

    t_map = by_level("t")
    u_map = by_level("u")
    v_map = by_level("v")

    def stack(mapping: Dict[float, np.ndarray], pressure_key: str) -> np.ndarray:
        missing = [p for p in levels if not any(abs(p - q) < 1e-6 for q in mapping)]
        if missing:
            raise KeyError("niveis ausentes para %s: %s" % (pressure_key, missing))
        return np.stack([mapping[min(mapping, key=lambda q: abs(q - p))] for p in levels], axis=0)

    temperature = stack(t_map, "t")
    u_wind = stack(u_map, "u")
    v_wind = stack(v_map, "v")

    if "q" in spec:
        humidity = stack(by_level("q"), "q")
    else:
        rh = stack(by_level("rh"), "rh")
        humidity = np.stack([
            specific_humidity_from_rh(rh[k], temperature[k], levels[k]) for k in range(len(levels))
        ], axis=0)

    sp = select_surface(fields, spec["sp"], ("surface", "isobaricInhPa"))
    if sp is None:
        raise KeyError("pressao de superficie ausente")
    msl = select_surface(fields, spec["msl"], ("meanSea", "surface"))

    topo = select_surface(fields, spec.get("topo", []), ("surface",))
    if topo is not None and source == "era5":
        topo = geopotential_to_height(topo)
    if topo is not None:
        topo = np.maximum(topo, 0.0)

    sst = select_surface(fields, spec["sst"], ("surface",))
    ice = select_surface(fields, spec["ci"], ("surface",))

    lat = next((f.latitudes for f in fields if _matches(f, spec["t"], level_type)), None)
    lon = next((f.longitudes for f in fields if _matches(f, spec["t"], level_type)), None)
    return Analysis(
        pressure_hpa=np.asarray(levels, dtype=np.float64),
        temperature=temperature, humidity=humidity, u_wind=u_wind, v_wind=v_wind,
        surface_pressure=to_hpa(sp),
        sea_level_pressure=None if msl is None else to_hpa(msl),
        topography=topo, sst=sst, ice=ice,
        latitudes=None if lat is None else np.asarray(lat, dtype=np.float64),
        longitudes=None if lon is None else np.asarray(lon, dtype=np.float64),
    )
