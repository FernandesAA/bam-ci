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
from .gribio import Analysis, build_analysis, read_grib, source_for_case
from .horizontal import interpolate_bilinear
from .ic_format import ICHeader, write_ic
from .levels import gaussian_latitudes, grid_shape, read_delta_sigma
from .spectral import analysis as spectral_analysis
from .spectral import file_order, smooth_coefficients, synthesis
from .vertical import clip_physical, new_ps, new_sigma, vert_sigma_inter

EARTH_R = 6.37e6

#: Tolerancia de positividade de q na grade sintetizada (kg/kg).
Q_POSITIVITY_TOL = 1e-8


def q_positive_coefficients(q_grid, mend, nlat, nlon, max_iter=3,
                            tol=Q_POSITIVITY_TOL, smooth=None):
    """Coeficientes de q com positividade na GRADE sintetizada (anti-ringing).

    Para cada nivel, parte dos coeficientes (opcionalmente suavizados por
    'smooth') e itera ate 'max_iter' vezes:

        c = analysis(max(synthesis(c), 0))

    O criterio de parada e' o minimo da GRADE sintetizada >= -tol (nao o
    coeficiente). NAO reaplica o filtro de Hoskins a cada iteracao: o filtro
    entra so no chute inicial, senao a convergencia fica muito lenta.

    A iteracao de recorte sozinha converge devagar para -1e-8 (ringing de
    alta ordem); por isso, se ainda restar min < -tol depois das iteracoes,
    aplica-se um DESLOCAMENTO CONSTANTE minimo delta = tol - min. Um campo
    constante e' exatamente representavel na base espectral (so o coeficiente
    (0,0)), logo a grade final fica com min == tol sem custo espectral. O
    preco e' somar delta ao q daquele nivel (medido em info["lift_max"] /
    ["lift_total"]); com N=3 o delta tipico e' ordens de grandeza menor que o
    minimo original.

    Parametros: q_grid (K, nlat, nlon) em kg/kg; mend/nlat/nlon da grade.
    Retorna (coef (K, (mend+1)(mend+2)) na ordem do arquivo, info) com
    info = {iterations (max por nivel), iterations_total, min (grade final),
            levels_below_tol, levels_lifted, lift_max, lift_total, tol}.
    Com max_iter=0 nao sintetiza (info["min"] fica None) e devolve o
    comportamento antigo, sem garantia de positividade.
    """
    q = np.asarray(q_grid, dtype=np.float64)
    if q.ndim != 3:
        raise ValueError("q_grid deve ser (kmax, nlat, nlon); recebido %s" % (q.shape,))
    kmax = int(q.shape[0])
    ncoeff = (int(mend) + 1) * (int(mend) + 2)
    coef = np.empty((kmax, ncoeff), dtype=np.float64)
    if max_iter <= 0:
        for k in range(kmax):
            c = spectral_analysis(q[k], mend)
            if smooth is not None:
                c = smooth(c)
            coef[k] = file_order(c)
        return coef, {"iterations": 0, "iterations_total": 0, "min": None,
                      "levels_below_tol": 0, "tol": float(tol)}

    total = 0
    max_used = 0
    min_final = np.inf
    levels_below = 0
    levels_lifted = 0
    lift_max = 0.0
    lift_total = 0.0
    for k in range(kmax):
        c = spectral_analysis(q[k], mend)
        if smooth is not None:
            c = smooth(c)
        qg = synthesis(c, mend, nlat, nlon)
        used = 0
        while float(np.min(qg)) < -tol and used < max_iter:
            c = spectral_analysis(np.maximum(qg, 0.0), mend)
            qg = synthesis(c, mend, nlat, nlon)
            used += 1
        m = float(np.min(qg))
        if m < -tol:
            # garantia final: campo constante e' exato na base espectral
            delta = tol - m
            c = spectral_analysis(qg + delta, mend)
            m = tol
            levels_lifted += 1
            lift_max = max(lift_max, delta)
            lift_total += delta
        if m < -tol:
            levels_below += 1
        coef[k] = file_order(c)
        total += used
        max_used = max(max_used, used)
        min_final = min(min_final, m)
    return coef, {"iterations": max_used, "iterations_total": total,
                  "min": float(min_final), "levels_below_tol": levels_below,
                  "levels_lifted": levels_lifted, "lift_max": float(lift_max),
                  "lift_total": float(lift_total), "tol": float(tol)}


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
        delta_sigma_path: Optional[str] = None, p_levels: Optional[Sequence[float]] = None,
        source: Optional[str] = None) -> PipelineResult:
    """Executa o pipeline e grava o arquivo de IC.

    `source` e o leitor GRIB ("era5" ou "nl"); quando omitido vem do caso
    (`Config.case`), via :func:`bamci.gribio.source_for_case`.
    """
    if analysis is None:
        if grib_path is None:
            raise ValueError("informe analysis ou grib_path")
        reader = source or source_for_case(cfg.case)
        fields = read_grib(grib_path)
        levels = list(p_levels) if p_levels is not None else sorted(
            {f.level for f in fields if f.level_type == "isobaricInhPa"}, reverse=True)
        if not levels:
            raise ValueError(
                "o GRIB nao tem campos em niveis de pressao (isobaricInhPa); "
                "confira se o download incluiu os campos 3D"
            )
        analysis = build_analysis(fields, reader, levels)

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
    # sanidade fisica antes de NewPs/Tv: q negativo/absurdo na entrada
    # envenenaria ps/tvir e a sintese espectral (bug TQ0213L064).
    t, q, clipped = clip_physical(t, q)
    if clipped["t_clipped"] or clipped["q_clipped"]:
        import warnings

        warnings.warn(
            "analise de entrada: T clipado em %d ponto(s) (min=%.3g, max=%.3g K) "
            "e q em %d ponto(s) (min=%.3g, max=%.3g kg/kg)"
            % (clipped["t_clipped"], clipped["t_min"], clipped["t_max"],
               clipped["q_clipped"], clipped["q_min"], clipped["q_max"]),
            RuntimeWarning, stacklevel=2)
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
    if delta.size != cfg.kmax_out:
        raise ValueError(
            "DeltaSigma %s tem %d camadas mas a resolucao pede %d "
            "(use -r T%s%04dL%03d ou --kmax %d)"
            % (delta_sigma_path, delta.size, cfg.kmax_out,
               cfg.grid, cfg.mend_out, delta.size, delta.size)
        )
    sig_inter, sig_layer = new_sigma(delta)
    p_out = sig_layer[:, None, None] * ps_out[None, :, :]

    remap = vert_sigma_inter({"u": u, "v": v, "t": t, "q": q}, analysis.pressure_hpa, p_out)
    tvir_sig = remap["t"] * (1.0 + cfg.ctvr * remap["q"])

    kmax_in = u.shape[0]
    div = np.empty_like(u)
    vort = np.empty_like(u)
    for k in range(kmax_in):
        div[k], vort[k] = spectral_div_vort(u[k], v[k], nlat, nlon, cfg.mend_out)
    # div/vort nao sao T/q: sem o clip fisico (que existe so para os campos
    # fisicos) e sem a extrapolacao de q -- usamos so o slot "t" interpolado.
    div_sig = vert_sigma_inter({"u": div, "v": div, "t": div, "q": div},
                               analysis.pressure_hpa, p_out,
                               physical=False)["t"]
    vort_sig = vert_sigma_inter({"u": vort, "v": vort, "t": vort, "q": vort},
                                analysis.pressure_hpa, p_out,
                                physical=False)["t"]

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
    # q: positividade na grade sintetizada (ringing de Gibbs apos a truncagem).
    # tvir/div/vort nao recebem esse passo (nao tem o problema).
    coef_q, q_positive_info = q_positive_coefficients(
        remap["q"], cfg.mend_out, nlat, nlon,
        max_iter=int(cfg.q_positivity), smooth=smooth)

    header = ICHeader(forecast_day=0, time_of_day=0.0,
                      date_inp=_date_tuple(cfg.date), date_cur=_date_tuple(cfg.date),
                      sigma_inter=sig_inter, sigma_layers=sig_layer)
    path = write_ic(cfg.output_path(), header,
                    file_order(coef_topo), file_order(coef_lnps),
                    coef_tvir.T, coef_div.T, coef_vort.T, coef_q.T)
    return PipelineResult(
        path=path, cfg=cfg,
        fields={"ps": ps_out, "topo": topo, "sigma": sig_layer,
                "q_positivity": q_positive_info})


def _date_tuple(text: str):
    return (int(text[8:10]), int(text[4:6]), int(text[6:8]), int(text[0:4]))
