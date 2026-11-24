"""Smoke tests da API publica de interpolacao vertical.

A implementacao fiel ao Fortran esta em tests/test_vertical_sigma.py; aqui so se
garante que a API exportada continua importavel e coerente.
"""

import numpy as np
import pytest

from bamci import interp_logp, new_ps, new_sigma, vert_sigma_inter


def test_public_api_roundtrip():
    p_in = np.array([1000.0, 850.0, 700.0, 500.0])
    field = np.tile(np.log(p_in)[:, None, None], (1, 1, 1))
    p_out = np.array([950.0, 800.0, 600.0])[:, None, None]
    got = interp_logp(field, p_in, p_out)
    np.testing.assert_allclose(got[:, 0, 0], np.log([950.0, 800.0, 600.0]), atol=1e-10)


def test_new_sigma_and_new_ps_smoke():
    inter, layer = new_sigma(np.array([0.5, 0.5]))
    assert inter[0] == 1.0 and inter[-1] == 0.0
    assert layer.shape == (2,)
    ps = new_ps(np.full((1, 1), 1013.25), np.array([1000.0, 500.0]),
                np.full((2, 1, 1), 285.0), np.zeros((1, 1)))
    assert ps[0, 0] > 900.0


# ---------------------------------------------------------------------------
# Extrapolacao nos topos/fundos e clip fisico (bug TQ0213L064)
# ---------------------------------------------------------------------------
import warnings

from bamci.vertical import clip_physical

P_ERA5 = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,
                   650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,
                   150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1], float)


def _era5_profiles():
    knots_p = np.array([1000, 850, 700, 500, 300, 200, 100, 50, 20, 10, 5, 3, 1], float)
    knots_t = np.array([288, 281, 270, 255, 230, 215, 210, 220, 240, 250, 255, 250, 235], float)
    t = np.interp(np.log(P_ERA5), np.log(knots_p[::-1]), knots_t[::-1])
    q = 0.012 * (P_ERA5 / 1000.0) ** 3
    return t, q


def _fields_at_levels(t, q, ny=2, nx=2):
    shape = (P_ERA5.size, ny, nx)
    return {
        "u": np.zeros(shape),
        "v": np.zeros(shape),
        "t": np.broadcast_to(t[:, None, None], shape).copy(),
        "q": np.broadcast_to(q[:, None, None], shape).copy(),
    }


def test_l64_topo_acima_da_entrada_era5():
    t, q = _era5_profiles()
    p_out = np.geomspace(1000.0, 0.1, 64)[:, None, None] * np.ones((1, 2, 2))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = vert_sigma_inter(_fields_at_levels(t, q), P_ERA5, p_out)
    assert np.all(np.isfinite(out["t"]))
    assert np.all(np.isfinite(out["q"]))
    assert 150.0 <= out["t"].min() and out["t"].max() <= 350.0
    assert out["q"].min() >= 0.0 and out["q"].max() <= 0.05
    # T isosterica e q -> 0 acima do nivel mais alto da entrada (1 hPa)
    above = p_out[:, 0, 0] < P_ERA5[-1]
    assert above.sum() > 0
    np.testing.assert_allclose(out["t"][above], t[-1], rtol=0.0, atol=1e-12)
    expected_top = q[-1] * p_out[above, 0, 0]
    np.testing.assert_allclose(out["q"][above][:, 0, 0], expected_top,
                               rtol=1e-12, atol=0.0)
    # q nao aumenta com a altura
    assert np.all(np.diff(out["q"][:, 0, 0]) <= 1e-12)
    assert not any("clipado" in str(w.message) for w in caught)


def test_l42_nao_regride():
    t, q = _era5_profiles()
    # topo real do L42 com ps~1013 hPa e' ~2.07 hPa: todo o alvo fica dentro
    # do intervalo da ERA5 (1..1000 hPa), sem extrapolacao
    p_out = np.geomspace(1013.0, 2.07, 42)[:, None, None] * np.ones((1, 2, 2))
    out = vert_sigma_inter(_fields_at_levels(t, q), P_ERA5, p_out)
    assert np.all(np.isfinite(out["t"]))
    assert np.all(np.isfinite(out["q"]))
    assert 150.0 <= out["t"].min() and out["t"].max() <= 350.0
    assert out["q"].min() >= 0.0 and out["q"].max() <= 0.05
    # dentro do dominio o interpolador continua exato para campo log-linear
    p_in = np.array([1000.0, 850.0, 700.0, 500.0])
    z = -np.log(p_in)
    field = np.tile((2.0 + 3.0 * z)[:, None, None], (1, 2, 2))
    p_test = np.array([950.0, 800.0, 600.0])[:, None, None] * np.ones((1, 2, 2))
    got = interp_logp(field, p_in, p_test)
    np.testing.assert_allclose(got[:, 0, 0], 2.0 + 3.0 * (-np.log(p_test[:, 0, 0])),
                               atol=1e-10)


def test_clip_fisico_avisa_e_limita():
    p_in = np.array([1000.0, 850.0, 700.0, 500.0, 300.0, 100.0, 10.0, 1.0])
    t = np.linspace(290.0, 220.0, p_in.size)
    q = np.linspace(0.012, 0.0, p_in.size)
    t[0] = 100.0          # absurdo no fundo
    q[0] = -0.01          # q negativo no fundo
    p_out = np.array([1000.0, 900.0, 700.0, 500.0, 300.0, 100.0, 10.0, 1.0])[:, None, None] \
        * np.ones((1, 2, 2))
    fields = {
        "u": np.zeros((p_in.size, 2, 2)), "v": np.zeros((p_in.size, 2, 2)),
        "t": np.broadcast_to(t[:, None, None], (p_in.size, 2, 2)).copy(),
        "q": np.broadcast_to(q[:, None, None], (p_in.size, 2, 2)).copy(),
    }
    with pytest.warns(RuntimeWarning, match="clipado"):
        out = vert_sigma_inter(fields, p_in, p_out)
    assert out["t"].min() >= 150.0 and out["t"].max() <= 350.0
    assert out["q"].min() >= 0.0 and out["q"].max() <= 0.05
    # physical=False preserva os valores (usado para div/vort)
    out2 = vert_sigma_inter(fields, p_in, p_out, physical=False)
    assert out2["t"].min() < 150.0
    assert out2["q"].min() < 0.0


def test_clip_physical_nao_finitos():
    t, q, info = clip_physical(
        np.array([np.nan, 100.0, 400.0]), np.array([np.nan, -1.0, 0.1]))
    assert np.all(np.isfinite(t)) and np.all(np.isfinite(q))
    assert t.min() >= 150.0 and t.max() <= 350.0
    assert q.min() >= 0.0 and q.max() <= 0.05
    assert info["t_clipped"] == 3 and info["q_clipped"] == 3


def test_interp_logp_bordas_sem_nan():
    p_in = np.array([1000.0, 500.0, 100.0, 1.0])
    field = np.array([1.0, 2.0, 3.0, 4.0])[:, None, None]
    p_out = np.array([2000.0, 1500.0, 0.5, 0.01])[:, None, None]
    got = interp_logp(field, p_in, p_out)
    assert np.all(np.isfinite(got))
    np.testing.assert_allclose(got[:, 0, 0], [1.0, 1.0, 4.0, 4.0])

