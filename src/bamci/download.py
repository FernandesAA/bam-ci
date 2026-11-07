"""Download das analises ERA5, FNL e GFS.

Decisao do autor (2026-11-07): **ERA5 e o default**, na grade nativa (0.25 graus,
37 niveis), **sem regrid**; FNL e GFS ficam como alternativas.

Fontes:
  era5 : CDS (cdsapi) -> GRIB, grade nativa
  fnl  : NOMADS (grib filter CGI)
  gfs  : NOMADS (arquivo direto pgrb2.0p25)

Nenhuma funcao aqui depende de rede para ser testada: a construcao de URLs e dos
pedidos e separada da transferencia.
"""

from __future__ import annotations

import os
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence

NOMADS = "https://nomads.ncep.noaa.gov"
CDS_URL = "https://cds.climate.copernicus.eu/api"

#: niveis de pressao do ERA5 (37) e do FNL/GFS (26), em hPa
ERA5_LEVELS = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600,
               550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 70,
               50, 30, 20, 10, 7, 5, 3, 2, 1]
NCEP_LEVELS = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450,
               400, 350, 300, 250, 200, 150, 125, 70, 50, 30, 20, 10]

ERA5_PL_VARIABLES = ["temperature", "specific_humidity", "u_component_of_wind",
                     "v_component_of_wind", "geopotential"]
ERA5_SL_VARIABLES = ["surface_pressure", "mean_sea_level_pressure",
                     "geopotential", "sea_surface_temperature", "sea_ice_cover"]


def parse_date(text: str) -> datetime:
    """Aceita YYYYMMDDHH e devolve datetime."""
    text = str(text).strip()
    if len(text) != 10 or not text.isdigit():
        raise ValueError("data deve ser YYYYMMDDHH: %r" % text)
    return datetime(int(text[:4]), int(text[4:6]), int(text[6:8]), int(text[8:10]))


# ---------------------------------------------------------------------------
# ERA5 (CDS)
# ---------------------------------------------------------------------------
def era5_pressure_request(date: datetime, levels: Sequence[int] = ERA5_LEVELS,
                          variables: Sequence[str] = ERA5_PL_VARIABLES,
                          flavor: str = "beta") -> Dict[str, object]:
    """Pedido CDS para os niveis de pressao, grade nativa (sem regrid)."""
    req: Dict[str, object] = {
        "product_type": ["reanalysis"],
        "variable": list(variables),
        "year": "%04d" % date.year,
        "month": "%02d" % date.month,
        "day": "%02d" % date.day,
        "time": ["%02d:00" % date.hour],
        "pressure_level": [str(int(p)) for p in levels],
    }
    req.update(_era5_format(flavor))
    return req


def era5_surface_request(date: datetime, variables: Sequence[str] = ERA5_SL_VARIABLES,
                         flavor: str = "beta") -> Dict[str, object]:
    """Pedido CDS para os campos de superficie, grade nativa."""
    req: Dict[str, object] = {
        "product_type": ["reanalysis"],
        "variable": list(variables),
        "year": "%04d" % date.year,
        "month": "%02d" % date.month,
        "day": "%02d" % date.day,
        "time": ["%02d:00" % date.hour],
    }
    req.update(_era5_format(flavor))
    return req


def _era5_format(flavor: str) -> Dict[str, str]:
    if flavor == "beta":
        return {"data_format": "grib", "download_format": "unarchived"}
    if flavor == "legacy":
        return {"format": "grib"}
    raise ValueError("flavor do CDS deve ser 'beta' ou 'legacy'")


def era5_paths(out_dir: Path, date: datetime) -> Dict[str, Path]:
    stamp = date.strftime("%Y%m%d%H")
    return {"pressure": Path(out_dir) / ("era5_pl_%s.grib" % stamp),
            "surface": Path(out_dir) / ("era5_sl_%s.grib" % stamp)}


def download_era5(out_dir: str | Path, date: datetime, levels: Sequence[int] = ERA5_LEVELS,
                  flavor: str = "beta", client=None) -> Dict[str, Path]:
    """Baixa ERA5 do CDS (precisa de ~/.cdsapirc). Sem regrid."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = era5_paths(out, date)
    if client is None:
        try:
            import cdsapi
        except ImportError as exc:  # pragma: no cover - ambiente
            raise RuntimeError("cdsapi nao instalado: pip install cdsapi e configure ~/.cdsapirc") from exc
        # na jaci a credencial vive em /oper/scripts/ioper/.cdsapirc
        fallback = "/oper/scripts/ioper/.cdsapirc"
        if not os.environ.get("CDSAPI_RC") and not (Path.home() / ".cdsapirc").exists()                 and Path(fallback).exists():
            os.environ["CDSAPI_RC"] = fallback
        client = cdsapi.Client(url=os.environ.get("CDS_URL", CDS_URL))
    client.retrieve("reanalysis-era5-pressure-levels", era5_pressure_request(date, levels, flavor=flavor),
                    str(paths["pressure"]))
    client.retrieve("reanalysis-era5-single-levels", era5_surface_request(date, flavor=flavor),
                    str(paths["surface"]))
    return paths


# ---------------------------------------------------------------------------
# FNL / GFS (NOMADS)
# ---------------------------------------------------------------------------
def fnl_url(date: datetime, base: str = NOMADS, mode: str = "filter") -> str:
    """URL do FNL no NOMADS (grib filter por padrao)."""
    stamp = date.strftime("%Y%m%d_%H_00")
    if mode == "direct":
        return "%s/pub/data/nccf/com/fnl/prod/fnl.%s/fnl_%s.grib2" % (
            base, date.strftime("%Y%m%d"), stamp)
    params = {
        "file": "fnl_%s.grib2" % stamp,
        "dir": "/fnl.%s" % date.strftime("%Y%m%d"),
        "all_var": "on",
        "all_lev": "on",
    }
    return "%s/cgi-bin/filter_fnl.pl?%s" % (base, urllib.parse.urlencode(params))


def gfs_url(date: datetime, base: str = NOMADS, resolution: str = "0p25") -> str:
    """URL de um arquivo GFS pgrb2 (analise f000)."""
    return "%s/pub/data/nccf/com/gfs/prod/gfs.%s/%02d/atmos/gfs.t%02dz.pgrb2.%s.f000" % (
        base, date.strftime("%Y%m%d"), date.hour, date.hour, resolution)


def ncep_path(out_dir: Path, source: str, date: datetime) -> Path:
    return Path(out_dir) / ("%s_%s.grib2" % (source, date.strftime("%Y%m%d%H")))


def download_url(url: str, dest: str | Path, timeout: float = 300.0) -> Path:
    """Baixa uma URL para o arquivo destino (sem validar conteudo)."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=timeout) as resp, tmp.open("wb") as fh:
        if getattr(resp, "status", 200) != 200:
            raise RuntimeError("HTTP %s em %s" % (resp.status, url))
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.replace(dest)
    return dest


def download_ncep(out_dir: str | Path, date: datetime, source: str = "gfs",
                  base: str = NOMADS, fnl_mode: str = "filter") -> Path:
    """Baixa FNL ou GFS do NOMADS."""
    if source == "fnl":
        url = fnl_url(date, base=base, mode=fnl_mode)
    elif source == "gfs":
        url = gfs_url(date, base=base)
    else:
        raise ValueError("source deve ser 'fnl' ou 'gfs'")
    return download_url(url, ncep_path(Path(out_dir), source, date))
