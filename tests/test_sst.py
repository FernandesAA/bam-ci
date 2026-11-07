import numpy as np
import pytest

from bamci.sst import UNDEF, to_model_grid, write_sstweekly


def test_to_model_grid_shape_and_order():
    lat = np.linspace(90, -90, 181)
    lon = np.arange(0, 360, 1.0)
    field = np.full((181, 360), 300.0)
    out = to_model_grid(field, lat, lon, nlat=48, nlon=96)
    assert out.shape == (48, 96)
    assert np.allclose(out[np.isfinite(out)], 300.0, atol=1e-6)


def test_land_becomes_undef():
    lat = np.linspace(90, -90, 19)
    lon = np.arange(0, 360, 10.0)
    field = np.full((19, 36), 300.0)
    field[:, :18] = np.nan          # metade "terra"
    out = to_model_grid(field, lat, lon, nlat=24, nlon=48)
    assert np.any(out <= UNDEF + 1)
    assert np.nanmax(out) == pytest.approx(300.0, abs=1e-6)


def test_write_two_records(tmp_path):
    ice = np.zeros((192, 384))
    sst = np.full((192, 384), 300.0)
    p = write_sstweekly(tmp_path / "SSTWeekly20261105.G00192", ice, sst)
    assert p.stat().st_size == 2 * 192 * 384 * 4
    raw = np.fromfile(p, dtype=">f4").reshape(2, 192, 384)
    assert raw[1, 0, 0] == pytest.approx(300.0)

def test_ice_from_sst_rule():
    from bamci.sst import ice_from_sst
    sst = np.array([270.0, 271.35, 272.6, 275.0])
    ice = ice_from_sst(sst)
    assert ice[0] == 1.0
    assert abs(ice[1] - 1.0) < 1e-12
    assert 0.0 < ice[2] < 1.0
    assert ice[3] == 0.0
