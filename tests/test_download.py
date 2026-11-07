from datetime import datetime
from pathlib import Path

import pytest

from bamci.download import (ERA5_LEVELS, NCEP_LEVELS, era5_paths,
                            era5_pressure_request, era5_surface_request,
                            fnl_url, gfs_url, ncep_path, parse_date)


def test_parse_date():
    assert parse_date("2026110500") == datetime(2026, 11, 5, 0)
    for bad in ("2026", "2026130500", "20261105", "abcdefghij"):
        with pytest.raises(ValueError):
            parse_date(bad)


def test_era5_pressure_request_is_native_grid():
    req = era5_pressure_request(datetime(2026, 11, 5, 0))
    assert "grid" not in req                         # sem regrid
    assert req["pressure_level"][0] == "1000" and req["pressure_level"][-1] == "1"
    assert len(req["pressure_level"]) == len(ERA5_LEVELS) == 37
    assert req["time"] == ["00:00"]
    assert req["data_format"] == "grib"
    assert req["product_type"] == ["reanalysis"]
    assert "temperature" in req["variable"]


def test_era5_legacy_flavor():
    req = era5_pressure_request(datetime(2026, 11, 5, 12), flavor="legacy")
    assert req["format"] == "grib" and "data_format" not in req
    assert req["time"] == ["12:00"]
    with pytest.raises(ValueError):
        era5_pressure_request(datetime(2026, 11, 5), flavor="v3")


def test_era5_surface_request():
    req = era5_surface_request(datetime(2026, 11, 5, 0))
    assert "sea_surface_temperature" in req["variable"]
    assert "pressure_level" not in req


def test_era5_paths(tmp_path):
    p = era5_paths(tmp_path, datetime(2026, 11, 5, 0))
    assert p["pressure"].name == "era5_pl_2026110500.grib"
    assert p["surface"].name == "era5_sl_2026110500.grib"


def test_fnl_url_filter_and_direct():
    d = datetime(2026, 11, 5, 0)
    filt = fnl_url(d)
    assert "filter_fnl.pl" in filt
    assert "file=fnl_20261105_00_00.grib2" in filt
    assert "dir=%2Ffnl.20261105" in filt
    direct = fnl_url(d, mode="direct")
    assert direct.endswith("/fnl.20261105/fnl_20261105_00_00.grib2")


def test_gfs_url():
    d = datetime(2026, 11, 5, 12)
    assert gfs_url(d).endswith("/gfs.20261105/12/atmos/gfs.t12z.pgrb2.0p25.f000")


def test_ncep_path(tmp_path):
    d = datetime(2026, 11, 5, 6)
    assert ncep_path(tmp_path, "fnl", d).name == "fnl_2026110506.grib2"
    assert ncep_path(tmp_path, "gfs", d).name == "gfs_2026110506.grib2"


def test_ncep_levels_are_26():
    assert len(NCEP_LEVELS) == 26 and NCEP_LEVELS[0] == 1000
