"""SST semanal (OISST v2.1) na grade do BAM.

Formato do arquivo que o modelo le (model/datain/SSTWeeklyYYYYMMDD.GZZZZZ):
UNFORMATTED com ACCESS='DIRECT' e recl = Imax*Jmax*4, sem marcadores, dois
registros de 192x384 float32 big endian:

    registro 1: mascara oceano/terra (+1 oceano, -1 terra)
    registro 2: SST em K, com -999 sobre terra

A grade e a gaussiana do modelo, com longitudes decrescentes (a mesma convencao
da condicao inicial).
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np

from .levels import gaussian_latitudes, grid_shape

UNDEF = -999.0


def read_oisst(path: str | Path, end_date: str, days: int = 7) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Media semanal da OISST v2.1 (dias terminando em end_date).

    Devolve (sst_k, ice, lat, lon) com NaN onde nao ha dado.
    """
    import xarray as xr

    ds = xr.open_dataset(path)
    time_name = "time"
    times = ds[time_name].values
    end = np.datetime64(end_date)
    sel = ds.sel(time=slice(end - np.timedelta64(days - 1, "D"), end))
    if sel.sizes.get(time_name, 0) == 0:
        raise ValueError("sem dados OISST no periodo terminando em %s" % end_date)
    sst = sel["sst"].mean(dim=time_name) if "sst" in sel else sel["sst"].mean(dim=time_name)
    ice = sel["ice"].mean(dim=time_name) if "ice" in sel else None
    lat = ds["lat"].values
    lon = ds["lon"].values
    ds.close()
    sst_k = np.asarray(sst.values, dtype=np.float64)
    if sst_k.max() < 200:      # veio em Celsius
        sst_k = sst_k + 273.15
    ice_arr = np.zeros_like(sst_k) if ice is None else np.nan_to_num(np.asarray(ice.values, dtype=np.float64))
    return sst_k, ice_arr, np.asarray(lat, dtype=np.float64), np.asarray(lon, dtype=np.float64)


def to_model_grid(field: np.ndarray, lat: np.ndarray, lon: np.ndarray,
                  nlat: int = 192, nlon: int = 384, fill=UNDEF,
                  reverse_lon: bool = False) -> np.ndarray:
    """Interpola para a grade gaussiana do modelo.

    Medido contra o SSTWeekly de producao: a grade da SST e lat norte->sul com
    longitudes CRESCENTES (0->360) — ao contrario dos campos espectrais da CI.
    """
    from .horizontal import interpolate_bilinear

    dst_lat = gaussian_latitudes(nlat)
    dst_lon = np.arange(nlon) * (360.0 / nlon)
    valid = np.isfinite(field)
    data = np.where(valid, field, 0.0)
    wgt = valid.astype(float)
    num = interpolate_bilinear(data, lat, lon, dst_lat, dst_lon)
    den = interpolate_bilinear(wgt, lat, lon, dst_lat, dst_lon)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[den < 0.5] = fill
    return out[..., ::-1] if reverse_lon else out


def write_sstweekly(path: str | Path, ice: np.ndarray, sst: np.ndarray) -> Path:
    """Grava o arquivo de SST semanal no formato do modelo."""
    ice = np.ascontiguousarray(np.asarray(ice, dtype=">f4"))
    sst = np.ascontiguousarray(np.asarray(sst, dtype=">f4"))
    if ice.shape != sst.shape or ice.ndim != 2:
        raise ValueError("ice e sst devem ter a mesma forma 2D")
    p = Path(path)
    with p.open("wb") as fh:
        fh.write(ice.tobytes())
        fh.write(sst.tobytes())
    return p


def build_weekly(oisst_path, end_date, out_path, nlat: int = 192, nlon: int = 384,
                 days: int = 7) -> Path:
    """Pipeline completo: OISST -> grade do modelo -> arquivo do modelo."""
    sst, ice, lat, lon = read_oisst(oisst_path, end_date, days=days)
    sst_g = to_model_grid(sst, lat, lon, nlat=nlat, nlon=nlon, fill=UNDEF)
    ice_g = to_model_grid(ice, lat, lon, nlat=nlat, nlon=nlon, fill=0.0)
    ice_g = np.clip(np.nan_to_num(ice_g), 0.0, 1.0)
    sst_g = np.where(np.isfinite(sst_g) & (sst_g > -900), sst_g, UNDEF)
    return write_sstweekly(out_path, ice_g, sst_g)

# ---------------------------------------------------------------------------
# RTG SST do NCEP (0.083 graus) — alternativa operacional a OISST
# ---------------------------------------------------------------------------
def read_rtg(paths) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Le um ou mais arquivos rtgssthr_grb_0.083.grib2 e devolve a media.

    Cada arquivo tem uma unica mensagem: T em superficie (K), grade 2160x4320.
    Sem fracao de gelo propria — ela e derivada da temperatura (regra do BAM:
    SST <= 271.35 K => gelo).
    """
    from .gribio import read_grib

    if isinstance(paths, (str, Path)):
        paths = [paths]
    fields = []
    lat = lon = None
    for p in paths:
        msgs = [m for m in read_grib(p) if m.level_type == "surface"]
        if not msgs:
            raise ValueError("sem campo de superficie em %s" % p)
        m = msgs[0]
        val = np.asarray(m.values, dtype=np.float64)
        if np.nanmax(val) < 200.0:      # veio em Celsius
            val = val + 273.15
        fields.append(np.where(val > 1e4, np.nan, val))
        lat, lon = m.latitudes, m.longitudes
    sst = np.nanmean(np.stack(fields, axis=0), axis=0)
    ice = ice_from_sst(sst)
    return sst, ice, np.asarray(lat), np.asarray(lon)


def ice_from_sst(sst: np.ndarray, freeze: float = 271.35, width: float = 2.5) -> np.ndarray:
    """Fracao de gelo a partir da SST (regra declarada do BAM)."""
    x = (np.asarray(sst, dtype=np.float64) - freeze) / width
    return np.clip(1.0 - x, 0.0, 1.0)


#: mascara de terra do modelo BAM (arquivo texto: 192 linhas de 384 caracteres, 1=terra)
DEFAULT_MASK = "/p/projetos/ioper/fix/BAM_hibrido/pre/datain/ModelLandSeaMask.G00192"


def load_land_mask(path: str | Path = DEFAULT_MASK, nlat: int = 192, nlon: int = 384) -> np.ndarray:
    """Mascara de terra do modelo (1 = terra, 0 = oceano).

    O arquivo e texto: uma linha por latitude, um caractere por longitude,
    384 caracteres + newline (73.920 bytes para 192x384). A contagem de '1'
    (24.469) coincide exatamente com os pontos -999 do SSTWeekly de producao.
    """
    text = Path(path).read_text()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if len(lines) != nlat:
        raise ValueError("mascara com %d linhas, esperado %d" % (len(lines), nlat))
    mask = np.array([[ch == "1" for ch in ln] for ln in lines], dtype=bool)
    if mask.shape != (nlat, nlon):
        raise ValueError("mascara com forma %s, esperado (%d, %d)" % (mask.shape, nlat, nlon))
    return mask


def build_weekly_rtg(paths, out_path, nlat: int = 192, nlon: int = 384,
                     topo_path: str | None = DEFAULT_MASK, days: int = 1) -> Path:
    """Pipeline com RTG: GRIB2 do NCEP -> grade do modelo -> arquivo do modelo.

    O RTG e dado PONTUAL (analise diaria): o arquivo do modelo usa o campo do dia
    da data, nao uma media. days > 1 so para diagnostico.
    """
    paths = list(paths)
    if days and days > 0:
        paths = paths[-days:]
    sst, ice, lat, lon = read_rtg(paths)
    sst_g = to_model_grid(sst, lat, lon, nlat=nlat, nlon=nlon, fill=UNDEF)
    ice_g = to_model_grid(ice, lat, lon, nlat=nlat, nlon=nlon, fill=0.0)
    ice_g = np.clip(np.nan_to_num(ice_g), 0.0, 1.0)
    if topo_path:
        land = load_land_mask(topo_path, nlat, nlon)
        sst_g = np.where(land, UNDEF, sst_g)
        # registro 1 e a mascara oceano/terra do modelo: +1 oceano, -1 terra
        ice_g = np.where(land, -1.0, 1.0)
    sst_g = np.where(np.isfinite(sst_g) & (sst_g > -900), sst_g, UNDEF)
    return write_sstweekly(out_path, ice_g, sst_g)
