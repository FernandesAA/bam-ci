import numpy as np
import pytest

from bamci.derivatives import A_EARTH, spectral_div_vort
from bamci.levels import gaussian_latitudes


def test_solid_body_rotation():
    nlat, nlon, mend = 96, 192, 40
    lat = gaussian_latitudes(nlat)
    lon = np.arange(nlon) * 360.0 / nlon
    u0 = 20.0
    phi = np.radians(lat)[:, None]
    u = np.full((nlat, nlon), 1.0) * (u0 * np.cos(phi))
    v = np.zeros((nlat, nlon))
    div, vort = spectral_div_vort(u, v, nlat, nlon, mend)
    expect = 2.0 * u0 * np.sin(phi) / A_EARTH
    mid = slice(4, nlat - 4)
    np.testing.assert_allclose(vort[mid], np.repeat(expect[mid], nlon, axis=1), rtol=0.02, atol=1e-8)
    assert np.max(np.abs(div)) < 1e-8


def test_zero_wind_gives_zero():
    nlat, nlon, mend = 32, 64, 8
    z = np.zeros((nlat, nlon))
    div, vort = spectral_div_vort(z, z, nlat, nlon, mend)
    assert np.max(np.abs(div)) == 0.0 and np.max(np.abs(vort)) == 0.0
