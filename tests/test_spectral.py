import numpy as np
import pytest

from bamci.levels import gaussian_latitudes
from bamci.spectral import (analysis, file_order, from_file_order, gauss_weights,
                            hoskins_weights, legendre_norm, synthesis)


def test_gauss_weights_sum_two():
    w = gauss_weights(192)
    assert abs(w.sum() - 2.0) < 1e-12
    assert w[0] == pytest.approx(w[-1])


def test_legendre_normalisation_invariants():
    lat = gaussian_latitudes(96)
    mu = np.sin(np.radians(lat))[:48]
    w = gauss_weights(96)[:48]
    pln = legendre_norm(20, mu)
    np.testing.assert_allclose(pln[0, 0], 1.0 / np.sqrt(2.0), rtol=0, atol=1e-14)
    for m in range(4):
        for n in range(m, 4):
            total = np.sum(2.0 * w * pln[m, n] ** 2)
            assert abs(total - 1.0) < 1e-10, (m, n, total)


def test_mean_invariant():
    nlat, nlon, mend = 32, 64, 5
    f = np.full((nlat, nlon), 7.5)
    coef = analysis(f, mend)
    assert coef[0, 0].real == pytest.approx(np.sqrt(2.0) * 7.5)
    assert abs(coef[0, 0].imag) < 1e-12


def test_roundtrip_band_limited():
    nlat, nlon, mend = 32, 64, 6
    rng = np.random.default_rng(3)
    coef = np.zeros((mend + 1, mend + 1), dtype=complex)
    for m in range(0, 4):
        for n in range(m, min(mend, 5) + 1):
            coef[m, n] = rng.normal() + 1j * rng.normal()
    coef[0, 0] = 2.0
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            if m == 0:
                coef[m, n] = coef[m, n].real
    field = synthesis(coef, mend, nlat, nlon)
    back = analysis(field, mend)
    for m in range(4):
        for n in range(m, 6):
            assert abs(back[m, n] - coef[m, n]) < 1e-10, (m, n, back[m, n], coef[m, n])


def test_file_order_is_m_major_and_invertible():
    mend = 6
    rng = np.random.default_rng(1)
    coef = np.zeros((mend + 1, mend + 1), dtype=complex)
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            coef[m, n] = rng.normal() + 1j * rng.normal()
    vec = file_order(coef)
    assert vec.size == (mend + 1) * (mend + 2)
    assert vec[0] == coef[0, 0].real and vec[1] == coef[0, 0].imag
    np.testing.assert_allclose(from_file_order(vec, mend), coef)


def test_analysis_rejects_odd_grid():
    with pytest.raises(ValueError):
        analysis(np.zeros((31, 64)), 5)

def test_hoskins_weights_match_production_ratios():
    """Razoes medidas entre a CI de producao e a nossa (sem filtro)."""
    w = hoskins_weights(126, 0.12, 10)
    # razoes medidas sao aproximadas (coeficientes zonais isolados); tolerancia 35%
    for n, measured in ((32, 0.999), (64, 0.913), (100, 0.450), (126, 0.085)):
        assert abs(w[n] - measured) < 0.35 * max(measured, 0.1), (n, w[n], measured)
    assert w[0] == pytest.approx(1.0)
    assert np.all(np.diff(w) <= 0)
