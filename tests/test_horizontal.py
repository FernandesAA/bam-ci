import numpy as np
import pytest

from bamci.horizontal import interpolate_bilinear, regular_grid


def test_regular_grid():
    lat, lon = regular_grid(3, 4)
    np.testing.assert_allclose(lat, [90.0, 0.0, -90.0])
    np.testing.assert_allclose(lon, [0.0, 90.0, 180.0, 270.0])


def test_bilinear_reproduces_linear_field():
    src_lat, src_lon = regular_grid(37, 72, 90.0, -90.0)
    a, b, c = 3.0, 0.5, 0.25
    field = a + b * src_lat[:, None] + c * src_lon[None, :]
    dst_lat = np.array([60.0, 10.0, -45.0])
    dst_lon = np.array([12.5, 100.0, 300.0])   # sem cruzar a costura 0/360
    got = interpolate_bilinear(field, src_lat, src_lon, dst_lat, dst_lon)
    expect = a + b * dst_lat[:, None] + c * dst_lon[None, :]
    np.testing.assert_allclose(got, expect, atol=1e-10)


def test_longitude_is_periodic():
    src_lat = np.array([10.0, 0.0, -10.0])
    src_lon = np.arange(0.0, 360.0, 2.5)
    field = np.cos(np.radians(src_lon))[None, :].repeat(3, axis=0)
    got = interpolate_bilinear(field, src_lat, src_lon, np.array([0.0]), np.array([359.5]))
    assert got[0, 0] == pytest.approx(np.cos(np.radians(359.5)), abs=2e-3)
    assert got[0, 0] > 0.999


def test_same_grid_is_identity():
    lat, lon = regular_grid(19, 36)
    rng = np.random.default_rng(1)
    field = rng.normal(size=(19, 36))
    got = interpolate_bilinear(field, lat, lon, lat, lon)
    np.testing.assert_allclose(got, field)


def test_polar_extrapolation_is_constant():
    src_lat = np.array([80.0, 0.0, -80.0])
    src_lon = np.array([0.0, 180.0])
    field = np.array([[5.0, 5.0], [1.0, 1.0], [7.0, 7.0]])
    got = interpolate_bilinear(field, src_lat, src_lon, np.array([89.0, -89.0]), src_lon)
    np.testing.assert_allclose(got[0], [5.0, 5.0])
    np.testing.assert_allclose(got[1], [7.0, 7.0])


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        interpolate_bilinear(np.zeros((3, 5)), np.linspace(90, -90, 3),
                             np.linspace(0, 360, 4, endpoint=False),
                             np.array([0.0]), np.array([0.0]))
