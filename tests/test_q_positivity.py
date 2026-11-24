"""Positividade de q depois da truncagem espectral (ringing de Gibbs).

O caso real T213L64 mostrou q negativo em 23 niveis (min -2e-4 kg/kg) mesmo com
a vertical correta: o campo continuo e' nao negativo, mas a sintese com
mend=213 faz ringing. Aqui testamos o passo de recorte iterativo + garantia
constante do pipeline.
"""

from __future__ import annotations

import numpy as np
import pytest

from bamci.config import Config
from bamci.gribio import Analysis
from bamci.ic_format import read_ic
from bamci.levels import grid_shape
from bamci.pipeline import (Q_POSITIVITY_TOL, q_positive_coefficients, run)
from bamci.spectral import from_file_order, synthesis

MEND = 42
NLON, NLAT = grid_shape(MEND)          # grade do modelo para T42


def _sharp_q():
    """q com degrau forte (0/0.02) no centro: trunca em T42 e faz ringing."""
    lat = np.linspace(90.0, -90.0, NLAT)
    lon = np.arange(NLON) * 360.0 / NLON
    latg, long = np.meshgrid(lat, lon, indexing="ij")
    patch = (np.abs(long - 180.0) < 30.0) & (np.abs(latg) < 30.0)
    return np.where(patch, 0.02, 0.0)[None, :, :]


def _q_min(coef_file_order):
    qg = synthesis(from_file_order(coef_file_order, MEND), MEND, NLAT, NLON)
    return float(np.min(qg)), qg


def test_q_positivity_t42():
    q = _sharp_q()
    c0, info0 = q_positive_coefficients(q, MEND, NLAT, NLON, max_iter=0)
    m0, _ = _q_min(c0[0])
    assert m0 < -Q_POSITIVITY_TOL          # sem o passo, q negativo
    assert info0["min"] is None            # N=0 desligado (nao sintetiza no helper)

    c3, info3 = q_positive_coefficients(q, MEND, NLAT, NLON, max_iter=3)
    m3, qg3 = _q_min(c3[0])
    assert m3 >= -Q_POSITIVITY_TOL         # criterio e' a GRADE sintetizada
    assert info3["min"] >= -Q_POSITIVITY_TOL
    assert 1 <= info3["iterations"] <= 3
    assert info3["levels_below_tol"] == 0
    # estabilidade: a correcao e' pequena e o campo nao explode
    # (o recorte pode aumentar um pouco o pico positivo por Gibbs; fica < 0.05)
    assert np.all(np.isfinite(qg3))
    assert 0.0 < info3["lift_max"] < 1e-2
    assert qg3.max() <= 0.05
    assert qg3.max() < 10.0 * 0.02


def _analysis_with_sharp_q():
    lat = np.linspace(90.0, -90.0, 64)
    lon = np.arange(128) * 360.0 / 128
    latg, long = np.meshgrid(lat, lon, indexing="ij")
    p = np.array([1000.0, 850.0, 500.0])
    t = np.broadcast_to(280.0 + 5.0 * np.cos(np.radians(latg)) * np.cos(np.radians(long)),
                        (p.size, 64, 128)).copy()
    q2d = np.where(np.abs(long - 180.0) < 25.0, 0.02, 0.0)
    q = np.broadcast_to(q2d, (p.size, 64, 128)).copy()
    zero = np.zeros((64, 128))
    return Analysis(
        pressure_hpa=p, temperature=t, humidity=q,
        u_wind=np.broadcast_to(5.0 + zero, (p.size, 64, 128)).copy(),
        v_wind=np.broadcast_to(zero, (p.size, 64, 128)).copy(),
        surface_pressure=np.full((64, 128), 1013.0),
        sea_level_pressure=np.full((64, 128), 1013.0),
        topography=np.zeros((64, 128)), latitudes=lat, longitudes=lon,
    )


@pytest.mark.parametrize("smooth_fields", [False, True])
def test_pipeline_q_positivity_does_not_touch_tvir(tmp_path, smooth_fields):
    delta = tmp_path / "DeltaSigma.L003"
    delta.write_text("0.333333 0.333333 0.333334\n", encoding="utf-8")
    analysis = _analysis_with_sharp_q()

    results = {}
    for n in (0, 3):
        out_dir = tmp_path / ("n%d" % n)
        out_dir.mkdir()
        cfg = Config(case="ERA5", date="2026102300", data_out=out_dir)
        cfg.set_resolution("TQ0042L003")
        cfg.smooth_fields = smooth_fields
        cfg.q_positivity = n
        results[n] = run(cfg, analysis=analysis, delta_sigma_path=str(delta))

    data0 = read_ic(results[0].path)
    data3 = read_ic(results[3].path)
    # tvir/div/vort/lnps nao sao tocados pelo passo de q
    for name in ("tvir", "divg", "vort", "lnps", "topo"):
        np.testing.assert_array_equal(np.asarray(data0[name]), np.asarray(data3[name]))
    # q mudou e a grade sintetizada ficou nao negativa
    assert not np.allclose(np.asarray(data0["sphu"]), np.asarray(data3["sphu"]))

    def q_mins(data):
        kmax = int(data["kmax"])
        out = []
        for k in range(kmax):
            qg = synthesis(from_file_order(np.asarray(data["sphu"])[:, k], MEND),
                           MEND, NLAT, NLON)
            out.append(float(qg.min()))
        return np.asarray(out)

    mins0 = q_mins(data0)
    mins3 = q_mins(data3)
    assert mins0.min() < -Q_POSITIVITY_TOL
    assert mins3.min() >= -Q_POSITIVITY_TOL

    info = results[3].fields["q_positivity"]
    assert info["min"] >= -Q_POSITIVITY_TOL
    assert info["levels_below_tol"] == 0
    assert 1 <= info["iterations"] <= 3
    # N=0: passo desligado, sem sintese no helper
    assert results[0].fields["q_positivity"]["min"] is None
