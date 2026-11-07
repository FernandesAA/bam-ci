import numpy as np
import pytest

from bamci.vertical import (DLTDZ, new_ps, new_sigma, interp_logp,
                            surface_fix, vert_sigma_inter)


def test_new_sigma_interfaces_and_layers():
    d = np.array([0.2, 0.3, 0.5])
    inter, layer = new_sigma(d)
    np.testing.assert_allclose(inter, [1.0, 0.8, 0.5, 0.0])
    assert np.all(np.diff(inter) < 0)
    assert inter[-1] == 0.0
    assert np.all((layer > 0) & (layer < 1))
    assert np.all(np.diff(layer) < 0)          # camadas de baixo para cima


def test_new_sigma_absorbs_rounding():
    d = np.array([0.3, 0.3, 0.3])              # soma 0.9
    inter, _ = new_sigma(d)
    assert inter[-1] == 0.0


def test_new_ps_at_sea_level():
    ps_in = np.array([[1013.25, 1013.25]])
    p = np.array([1000.0, 850.0, 700.0, 500.0])
    tv = np.full((4, 1, 2), 290.0)
    topo = np.zeros((1, 2))
    ps = new_ps(ps_in, p, tv, topo)
    assert ps.shape == (1, 2)
    assert np.all(np.isfinite(ps))
    assert abs(ps[0, 0] - 1013.25) < 5.0       # ~igual a pressao ao nivel do mar


def test_new_ps_decreases_with_topography():
    p = np.array([1000.0, 850.0, 700.0, 500.0])
    tv = np.full((4, 1, 3), 285.0)
    ps_in = np.full((1, 3), 1013.25)
    topo = np.array([[0.0, 1000.0, 3000.0]])
    ps = new_ps(ps_in, p, tv, topo)[0]
    assert ps[0] > ps[1] > ps[2]
    assert ps[2] < 750.0                       # 3000 m ~ 700 hPa


def test_interp_is_exact_for_log_linear_field():
    p_in = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])
    z = -np.log(p_in)
    field = np.tile((2.0 + 3.0 * z)[:, None, None], (1, 2, 3))
    z_target = np.array([[z[0] + 0.3, z[1] + 0.2, z[2] + 0.1]] * 2)
    p_out = np.exp(-z_target)[None, :, :]
    got = interp_logp(field, p_in, p_out)
    np.testing.assert_allclose(got[0], 2.0 + 3.0 * (-np.log(p_out[0])), atol=1e-10)


def test_interp_is_constant_for_constant_field():
    p_in = np.array([1000.0, 500.0, 100.0])
    field = np.full((3, 2, 2), 7.0)
    p_out = np.full((4, 2, 2), 700.0)
    p_out[1] = 400.0
    np.testing.assert_allclose(interp_logp(field, p_in, p_out), 7.0)


def test_interp_cubic_does_not_overshoot():
    p_in = np.array([1000.0, 850.0, 700.0, 500.0, 300.0, 100.0])
    prof = np.array([10.0, 10.0, 10.0, 20.0, 20.0, 20.0])
    field = np.tile(prof[:, None, None], (1, 1, 1))
    p_out = np.full((3, 1, 1), 600.0)
    p_out[1] = 750.0
    p_out[2] = 400.0
    got = interp_logp(field, p_in, p_out)
    assert got.min() >= 10.0 - 1e-12 and got.max() <= 20.0 + 1e-12


def test_interp_above_top_uses_top_value():
    p_in = np.array([1000.0, 500.0, 100.0])
    field = np.array([1.0, 2.0, 3.0])[:, None, None]
    got = interp_logp(field, p_in, np.array([50.0, 1.0])[:, None, None])
    np.testing.assert_allclose(got[:, 0, 0], [3.0, 3.0])


def test_surface_fix_first_layer_uses_input():
    p_in = np.array([1000.0, 900.0, 800.0])
    p_out = np.array([[1050.0], [950.0], [700.0]])        # 1050 e 950 abaixo de 1000
    t1 = np.array([[290.0], [285.0], [280.0]])
    q1 = np.array([[0.012], [0.010], [0.008]])
    t2 = np.array([[289.0], [284.0], [270.0]])
    q2 = np.array([[0.011], [0.009], [0.004]])
    tf, qf = surface_fix(t2, q2, p_in, p_out, t1, q1)
    assert tf[0, 0] == pytest.approx(290.0)                # camada 1 = t1(1)
    assert qf[0, 0] == pytest.approx(0.012)
    # ultima camada (p=700 > 800? nao) permanece interpolada
    assert tf[2, 0] == pytest.approx(270.0)


def test_surface_fix_last_layer_lapse_rate():
    p_in = np.array([1000.0, 500.0, 100.0])
    p_out = np.array([[1100.0], [1050.0], [1020.0]])   # abaixo de p_in[0]=1000
    t1 = np.array([[300.0], [280.0], [250.0]])
    q1 = np.array([[0.02], [0.01], [0.004]])
    t2 = np.zeros_like(p_out)
    q2 = np.zeros_like(p_out)
    tf, qf = surface_fix(t2, q2, p_in, p_out, t1, q1)
    dz = -np.log(p_out[-1, 0]) - (-np.log(p_in[0]))
    assert tf[-1, 0] == pytest.approx(300.0 * np.exp(DLTDZ * dz))
    assert tf[-1, 0] > 300.0        # abaixo da superficie aquece a -6.5 K/km


def test_vert_sigma_inter_shapes_and_surface():
    p_in = np.array([1000.0, 850.0, 700.0, 500.0])
    ny, nx, k2 = 2, 3, 5
    rng = np.random.default_rng(0)
    fields = {name: rng.normal(size=(4, ny, nx)) for name in "uvtq"}
    ps = np.full((ny, nx), 1010.0)
    delta = np.full(k2, 1.0 / k2)
    _, layers = new_sigma(delta)
    p_out = layers[:, None, None] * ps[None]
    out = vert_sigma_inter(fields, p_in, p_out)
    assert set(out) == set("uvtq")
    for name in out:
        assert out[name].shape == (k2, ny, nx)
    assert np.all(np.isfinite(out["t"]))
