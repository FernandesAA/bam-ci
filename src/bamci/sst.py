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

import os
from pathlib import Path
from typing import Tuple

import numpy as np

from .levels import gaussian_latitudes, grid_shape

UNDEF = -999.0


def normalize_date(text: str) -> str:
    """Aceita YYYYMMDD, YYYYMMDDHH, YYYY-MM-DD e YYYY-MM-DDTHH -> `YYYYMMDD`.

    O `--date` do `bamci run` tem 10 digitos (com hora); o SST e diario, entao a
    hora e simplesmente descartada em vez de virar erro de formato.
    """
    raw = str(text).strip()
    if raw.endswith("Z"):
        raw = raw[:-1]
    if "T" in raw:
        raw = raw.split("T", 1)[0]
    digits = raw.replace("-", "").replace("/", "").replace(":", "").replace(" ", "")
    if not digits.isdigit():
        raise ValueError("data invalida: %r (use YYYYMMDD ou YYYYMMDDHH)" % (text,))
    if len(digits) == 10:
        digits = digits[:8]
    elif len(digits) == 12:
        digits = digits[:8]
    if len(digits) != 8:
        raise ValueError("data invalida: %r (use YYYYMMDD ou YYYYMMDDHH)" % (text,))
    from datetime import datetime

    try:
        datetime.strptime(digits, "%Y%m%d")     # datetime valida mes/dia de verdade
    except ValueError:
        raise ValueError("data inexistente: %r" % (text,)) from None
    return digits


def open_dataset(path: str | Path):
    """Abre o NetCDF com mensagem acionavel quando falta o motor de leitura."""
    try:
        import xarray as xr
    except ImportError as exc:  # pragma: no cover - depende do ambiente
        raise RuntimeError(
            "xarray nao instalado; instale os extras de leitura:\n"
            "  pip install 'bam-ci[netcdf]'"
        ) from exc

    try:
        return xr.open_dataset(path)
    except (ValueError, ImportError) as exc:
        text = str(exc).lower()
        if "backend" in text or "dependenc" in text or "engine" in text:
            raise RuntimeError(
                "nao ha motor NetCDF instalado para ler\n  %s\n"
                "Instale um deles e repita:\n"
                "  pip install netCDF4      (recomendado)\n"
                "  pip install h5netcdf     (alternativa, precisa de h5py)\n"
                "ou: pip install 'bam-ci[netcdf]'"
                % (path,)
            ) from exc
        raise


def _date64(text: str) -> np.datetime64:
    """`"20261112"` -> `np.datetime64('2026-11-12')`.

    Cuidado: `np.datetime64("20261112", "D")` NAO le como data — o numpy
    interpreta os 8 digitos como **ano** e devolve `20261112-01-01`, o que faz
    toda comparacao de janela falhar. Por isso a conversao passa por datetime.
    """
    from datetime import datetime

    return np.datetime64(datetime.strptime(str(text), "%Y%m%d").date())


def _as_datetimes(values) -> np.ndarray | None:
    """Converte o eixo de tempo para datetime64[D]; None se nao for possivel."""
    arr = np.asarray(values)
    if np.issubdtype(arr.dtype, np.datetime64):
        return arr.astype("datetime64[D]")
    try:
        return np.array([np.datetime64(str(v)[:10]) for v in arr], dtype="datetime64[D]")
    except (ValueError, TypeError):
        return None


def _select_window(ds, time_name: str, end_date: str, days: int):
    """Seleciona os tempos da janela SEM depender de indice monotonico.

    `ds.sel(time=slice(...))` exige indice monotonico (e cita o rotulo na
    mensagem de erro): arquivos de media semanal baixados fora de ordem
    quebravam com "non-monotonic index". Aqui o recorte e posicional
    (`isel` sobre uma mascara booleana), que funciona em qualquer ordem.
    """
    times = _as_datetimes(ds[time_name].values)
    if times is None:
        return ds.isel({time_name: slice(0, 0)})
    end = _date64(end_date)
    lo = end - np.timedelta64(max(int(days) - 1, 0), "D")
    mask = (times >= lo) & (times <= end)
    return ds.isel({time_name: np.flatnonzero(mask)})


def _nearest_time(values, end_date: str):
    """(indice, rotulo) do tempo mais proximo da data pedida."""
    times = _as_datetimes(values)
    if times is None:
        raise ValueError(
            "nao foi possivel interpretar o eixo de tempo: %r" % (np.asarray(values)[:3],)
        )
    end = _date64(end_date)
    idx = int(np.argmin(np.abs(times - end)))
    return idx, str(times[idx])


def read_oisst(path: str | Path, end_date: str, days: int = 7) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Media semanal da OISST v2.1 (dias terminando em end_date).

    Devolve (sst_k, ice, lat, lon) com NaN onde nao ha dado. Se o arquivo ja e
    uma media semanal (um unico tempo), usa esse tempo em vez de falhar.
    """
    import warnings

    ds = open_dataset(path)
    time_name = "time"
    if time_name not in ds.dims:
        raise ValueError("arquivo sem dimensao de tempo: %s" % path)
    end_date = normalize_date(end_date)
    sel = _select_window(ds, time_name, end_date, days)
    if sel.sizes.get(time_name, 0) == 0:
        idx, label = _nearest_time(np.asarray(ds[time_name].values), end_date)
        warnings.warn(
            "sem dado na janela de %d dias ate %s; usando o tempo mais proximo "
            "do arquivo (%s), que pode ser uma media semanal ja calculada"
            % (days, end_date, label),
            RuntimeWarning,
            stacklevel=2,
        )
        sel = ds.isel({time_name: [idx]})
    sst = sel["sst"].mean(dim=time_name)
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


def build_weekly(oisst_path, end_date, out_path, nlat: int | None = None,
                 nlon: int | None = None, days: int = 7, mend: int = 126) -> Path:
    """Pipeline completo: OISST -> grade do modelo -> arquivo do modelo.

    A grade vem de `mend` (T126 -> 384x192) quando `nlat`/`nlon` nao sao dados.
    """
    if nlat is None or nlon is None:
        nlat, nlon = model_grid(mend)
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


#: diretorio das mascaras de terra do modelo BAM (arquivo texto, 1=terra)
MASK_DIR = "/p/projetos/ioper/fix/BAM_hibrido/pre/datain"

#: mascara padrao historica (TQ126 -> grade 384x192)
DEFAULT_MASK = MASK_DIR + "/ModelLandSeaMask.G00192"


#: nome do arquivo que o modelo le em model/datain
SST_TEMPLATE = "SSTWeekly{date}.G{nlat:05d}"


def resolve_sst_output(out: str | Path, date: str | None, nlat: int = 192) -> Path:
    """Aceita `--out` como arquivo ou como diretorio.

    Se `out` for um diretorio existente (ou terminar em `/`), monta
    `SSTWeekly<AAAAMMDD>.G00<nlat>` dentro dele. Se for um caminho de arquivo,
    respeita o nome dado.
    """
    raw = str(out)
    path = Path(out)
    looks_dir = path.is_dir() or raw.endswith("/") or raw.endswith(os.sep)
    if not looks_dir:
        return path
    if date is None:
        raise ValueError(
            "--out foi um diretorio e sem --date nao ha como montar o nome "
            "(SSTWeekly<AAAAMMDD>.G00%d)" % nlat
        )
    return path / SST_TEMPLATE.format(date=normalize_date(date), nlat=int(nlat))


def model_grid(mend: int = 126) -> tuple:
    """(nlat, nlon) da grade gaussiana do modelo para a truncacao `mend`."""
    nlon, nlat = grid_shape(mend)
    return int(nlat), int(nlon)


def default_mask(mend: int = 126) -> str:
    """Caminho da mascara de terra para uma truncacao (G00<latitudes>).

    O numero do arquivo e a quantidade de latitudes da grade gaussiana:
    T126 -> 192 (G00192), T213 -> 320 (G00320), T254 -> 384, T666 -> 1000.
    """
    nlat, _ = model_grid(mend)
    return "%s/ModelLandSeaMask.G%05d" % (MASK_DIR, nlat)


def load_land_mask(path: str | Path | None = None, nlat: int = 192,
                   nlon: int = 384) -> np.ndarray:
    """Mascara de terra do modelo (1 = terra, 0 = oceano).

    `path=None` usa :func:`default_mask` para a grade pedida.
    """
    if path is None:
        path = default_mask_for_shape(nlat)
    return _load_land_mask_file(Path(path), nlat, nlon)


def default_mask_for_shape(nlat: int) -> str:
    """Mascara cujo numero casa com a quantidade de latitudes."""
    return "%s/ModelLandSeaMask.G%05d" % (MASK_DIR, int(nlat))


def _load_land_mask_file(path: Path, nlat: int, nlon: int) -> np.ndarray:
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


def build_weekly_rtg(paths, out_path, nlat: int | None = None, nlon: int | None = None,
                     topo_path: str | None = None, days: int = 1,
                     mend: int = 126) -> Path:
    """Pipeline com RTG: GRIB2 do NCEP -> grade do modelo -> arquivo do modelo.

    O RTG e dado PONTUAL (analise diaria): o arquivo do modelo usa o campo do dia
    da data, nao uma media. days > 1 so para diagnostico.

    `mend` define a grade (T126 -> 384x192) e, quando `topo_path` nao e dado, a
    mascara padrao dessa grade (:func:`default_mask`).
    """
    if nlat is None or nlon is None:
        nlat, nlon = model_grid(mend)
    if topo_path is None:
        topo_path = default_mask(mend)
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
