from pathlib import Path

import numpy as np

from bamci.levels import (fft_friendly_size, gaussian_latitudes, grid_shape,
                          read_delta_sigma, read_hybrid, read_pressure_levels,
                          sigma_interfaces, sigma_layers)


def test_fft_friendly():
    assert fft_friendly_size(379) == 384
    assert fft_friendly_size(187) == 192
    assert fft_friendly_size(384) == 384


def test_grid_shape_matches_bam():
    assert grid_shape(126) == (384, 192)   # TQ126L042 do subsazonal
    assert grid_shape(62) == (192, 96)     # TQ062L028 legado


def test_gaussian_latitudes():
    lat = gaussian_latitudes(192)
    assert lat.shape == (192,)
    assert lat[0] > 89.0 and lat[-1] < -89.0
    np.testing.assert_allclose(lat, -lat[::-1], atol=1e-10)
    nodes, weights = np.polynomial.legendre.leggauss(192)
    assert abs(weights.sum() - 2.0) < 1e-12


def test_read_hybrid_pa_to_hpa(tmp_path):
    p = tmp_path / "HybridLevels.L002"
    p.write_text("0 1.0\n2000 0.5\n10000 0.0\n")
    a, b = read_hybrid(p)
    np.testing.assert_allclose(a, [0.0, 20.0, 100.0])   # 2000 Pa -> 20 hPa
    np.testing.assert_allclose(b, [1.0, 0.5, 0.0])


def test_read_pressure_levels(tmp_path):
    p = tmp_path / "PressureLevels.ERA5"
    p.write_text("  1000.0  975.0  950.0\n  925.0  900.0\n")
    np.testing.assert_allclose(read_pressure_levels(p), [1000, 975, 950, 925, 900])


def test_delta_sigma_normalised(tmp_path):
    p = tmp_path / "DeltaSigma.L002"
    p.write_text("0.500000 0.400000\n")
    d = read_delta_sigma(p)
    np.testing.assert_allclose(d.sum(), 1.0)
    inter = sigma_interfaces(d)
    np.testing.assert_allclose(inter, [0.0, d[0], 1.0])
    np.testing.assert_allclose(sigma_layers(d), 0.5 * (inter[:-1] + inter[1:]))
