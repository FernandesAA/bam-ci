"""Pipeline completo: analise GRIB -> condicao inicial espectral do BAM.

Etapas (equivalentes ao GetICnAGCM/runPre do CPTEC):
  1. ler GRIB (ERA5/FNL/GFS) e montar os campos;
  2. interpolar horizontalmente para a grade gaussiana TQ126 (384x192);
  3. recalcular a pressao de superficie com a topografia (NewPs);
  4. interpolar verticalmente para os niveis sigma do modelo (L42);
  5. calcular divergencia e vorticidade (diferencas finitas na grade);
  6. decompor em harmonicos esfericos;
  7. escrever o arquivo UNFORMATTED big-endian.

Limitacao conhecida: divergencia/vorticidade usam diferencas finitas (a versao
Fortran usa derivadas espectrais); o resto segue as formulas do legado.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np

from .config import Config
from .derivatives import spectral_div_vort
from .gribio import Analysis, build_analysis, read_grib
from .horizontal import interpolate_bilinear
from .ic_format import ICHeader, write_ic
from .levels import gaussian_latitudes, grid_shape, read_delta_sigma
from .spectral import analysis as spectral_analysis
from .spectral import file_order, smooth_coefficients
from .vertical import new_ps, new_sigma, vert_sigma_inter

EARTH_R = 6.37e6


def divergence_vorticity(u, v, lat, lon) -> tuple:
    """Divergencia e vorticidade por diferencas finitas centradas (grade regular)."""
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    lat_r = np.radians(np.asarray(lat, dtype=np.float64))
    dlon = np.radians(float(lon[1] - lon[0]))
    dlat = np.radians(abs(float(lat[1] - lat[0])))
    coslat = np.cos(lat_r)[:, None]

    def ddx(f):
        return (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2.0 * dlon)

    def ddy(f):
        g = np.gradient(f, axis=0) / dlat
        return g

    div = ddx(u) / (EARTH_R * coslat) + ddy(v) / EARTH_R - v * np.tan(lat_r)[:, None] / EARTH_R
    vort = ddx(v) / (EARTH_R * coslat) - ddy(u) / EARTH_R + u * np.tan(lat_r)[:, None] / EARTH_R
    return div, vort


def _to_model_grid(field, lat_in, lon_in, lat_out, lon_out):
    return interpolate_bilinear(field, lat_in, lon_in, lat_out, lon_out)


@dataclass
class PipelineResult:
    path: Path
    cfg: Config
    fields: Dict[str, object]


def run(cfg: Config, analysis: Optional[Analysis] = None, grib_path: Optional[str] = None,
        delta_sigma_path: Optional[str] = None, p_levels: Optional[Sequence[float]] = None) -> PipelineResult:
    """Executa o pipeline e grava o arquivo de IC."""
    if analysis is None:
        if grib_path is None:
            raise ValueError("informe analysis ou grib_path")
        fields = read_grib(grib_path)
        levels = list(p_levels) if p_levels is not None else sorted(
            {f.level for f in fields if f.level_type == "isobaricInhPa"}, reverse=True)
        analysis = build_analysis(fields, "nl", levels)

    nlon, nlat = grid_shape(cfg.mend_out)
    lat_out = gaussian_latitudes(nlat)
    lon_out = np.arange(nlon) * (360.0 / nlon)
    lat_in = analysis.latitudes
    lon_in = analysis.longitudes
    if lat_in is None or lon_in is None:
        raise ValueError("analysis sem latitudes/longitudes")

    shift = int(round(nlon * (cfg.lon_shift / 360.0))) % nlon

    def roll(field):
        out = field
        if cfg.lon_reverse:
            out = out[..., ::-1]
        if shift:
            out = np.roll(out, shift, axis=-1)
        return out

    def hz(field):
        return np.stack([roll(_to_model_grid(field[k], lat_in, lon_in, lat_out, lon_out))
                         for k in range(field.shape[0])], axis=0)

    t = hz(analysis.temperature)
    q = hz(analysis.humidity)
    u = hz(analysis.u_wind)
    v = hz(analysis.v_wind)
    tvir = t * (1.0 + cfg.ctvr * q)

    if analysis.topography is None:
        raise ValueError("analysis sem topografia")
    topo = roll(_to_model_grid(analysis.topography, lat_in, lon_in, lat_out, lon_out))

    ps_in = analysis.sea_level_pressure if analysis.sea_level_pressure is not None else analysis.surface_pressure
    ps_in = roll(_to_model_grid(ps_in, lat_in, lon_in, lat_out, lon_out))
    ps_out = new_ps(ps_in, analysis.pressure_hpa, tvir, topo)

    if delta_sigma_path is None:
        raise ValueError("informe o arquivo DeltaSigma.Lxxx")
    delta = read_delta_sigma(delta_sigma_path)
    sig_inter, sig_layer = new_sigma(delta)
    p_out = sig_layer[:, None, None] * ps_out[None, :, :]

    remap = vert_sigma_inter({"u": u, "v": v, "t": t, "q": q}, analysis.pressure_hpa, p_out)
    tvir_sig = remap["t"] * (1.0 + cfg.ctvr * remap["q"])

    kmax_in = u.shape[0]
    div = np.empty_like(u)
    vort = np.empty_like(u)
    for k in range(kmax_in):
        div[k], vort[k] = spectral_div_vort(u[k], v[k], nlat, nlon, cfg.mend_out)
    div_sig = vert_sigma_inter({"u": div, "v": div, "t": div, "q": div},
                               analysis.pressure_hpa, p_out)["t"]
    vort_sig = vert_sigma_inter({"u": vort, "v": vort, "t": vort, "q": vort},
                                analysis.pressure_hpa, p_out)["t"]

    lnps = np.log(ps_out / 10.0)
    coef_topo = spectral_analysis(topo, cfg.mend_out)
    coef_lnps = spectral_analysis(lnps, cfg.mend_out)
    if cfg.smooth_fields:
        coef_topo = smooth_coefficients(coef_topo, cfg.mend_out, cfg.smooth_percent, cfg.iter_smooth)
        coef_lnps = smooth_coefficients(coef_lnps, cfg.mend_out, cfg.smooth_percent, cfg.iter_smooth)

    kmax = int(tvir_sig.shape[0])
    coef_tvir = np.zeros((kmax, (cfg.mend_out + 1) * (cfg.mend_out + 2)))
    coef_div = np.zeros_like(coef_tvir)
    coef_vort = np.zeros_like(coef_tvir)
    coef_q = np.zeros_like(coef_tvir)
    def smooth(c):
        return smooth_coefficients(c, cfg.mend_out, cfg.smooth_percent, cfg.iter_smooth)             if cfg.smooth_fields else c

    for k in range(kmax):
        coef_tvir[k] = file_order(smooth(spectral_analysis(tvir_sig[k], cfg.mend_out)))
        coef_div[k] = file_order(smooth(spectral_analysis(div_sig[k], cfg.mend_out)))
        coef_vort[k] = file_order(smooth(spectral_analysis(vort_sig[k], cfg.mend_out)))
        coef_q[k] = file_order(smooth(spectral_analysis(remap["q"][k], cfg.mend_out)))

    header = ICHeader(forecast_day=0, time_of_day=0.0,
                      date_inp=_date_tuple(cfg.date), date_cur=_date_tuple(cfg.date),
                      sigma_inter=sig_inter, sigma_layers=sig_layer)
    path = write_ic(cfg.output_path(), header,
                    file_order(coef_topo), file_order(coef_lnps),
                    coef_tvir.T, coef_div.T, coef_vort.T, coef_q.T)
    return PipelineResult(path=path, cfg=cfg, fields={"ps": ps_out, "topo": topo,
                                                      "sigma": sig_layer})


def _date_tuple(text: str):
    return (int(text[8:10]), int(text[4:6]), int(text[6:8]), int(text[0:4]))
