import numpy as np
import pytest

from bamci.gribio import (Analysis, RawField, build_analysis,
                          geopotential_to_height, specific_humidity_from_rh)


def _field(name, level_type, level, values, units=""):
    ny, nx = values.shape
    return RawField(short_name=name, level_type=level_type, level=level,
                    values=values, latitudes=np.linspace(90, -90, ny),
                    longitudes=np.linspace(0, 360, nx, endpoint=False), units=units)


def _ncep_fields(levels=(1000.0, 850.0, 500.0)):
    ny, nx = 3, 4
    base = np.linspace(250.0, 300.0, ny * nx).reshape(ny, nx)
    fields = []
    for k, p in enumerate(levels):
        fields.append(_field("t", "isobaricInhPa", p, base + k, "K"))
        fields.append(_field("u", "isobaricInhPa", p, np.full((ny, nx), 5.0 + k), "m/s"))
        fields.append(_field("v", "isobaricInhPa", p, np.full((ny, nx), -2.0), "m/s"))
        if p >= 100.0:
            fields.append(_field("r", "isobaricInhPa", p, np.full((ny, nx), 50.0), "%"))
    fields.append(_field("pres", "surface", 0.0, np.full((ny, nx), 101325.0), "Pa"))
    fields.append(_field("prmsl", "meanSea", 0.0, np.full((ny, nx), 101500.0), "Pa"))
    fields.append(_field("hgt", "surface", 0.0, np.full((ny, nx), 500.0), "gpm"))
    fields.append(_field("icec", "surface", 0.0, np.zeros((ny, nx))))
    fields.append(_field("t", "surface", 0.0, np.full((ny, nx), 290.0), "K"))
    return fields


def test_specific_humidity_formula():
    t = np.array([[300.0]])
    q = specific_humidity_from_rh(np.array([[50.0]]), t, 1000.0)
    es = 6.1078 * np.exp(17.2693882 * (300.0 - 273.16) / (300.0 - 35.86))
    e = 0.5 * es
    expect = 0.622 * e / (1000.0 - (1 - 0.622) * e)
    assert abs(q[0, 0] - expect) < 1e-15
    assert 0.005 < q[0, 0] < 0.02


def test_specific_humidity_uses_pressure_in_hpa():
    t = np.full((1, 1), 280.0)
    q_hpa = specific_humidity_from_rh(np.full((1, 1), 60.0), t, 850.0)
    q_pa = specific_humidity_from_rh(np.full((1, 1), 60.0), t, 85000.0)
    assert q_hpa[0, 0] > 0.004
    # passar a pressao em Pa subestima q por um fator ~100 (o erro do legado no grib2)
    assert 80.0 < q_hpa[0, 0] / q_pa[0, 0] < 120.0


def test_build_analysis_ncep():
    fields = _ncep_fields()
    an = build_analysis(fields, "nl", [1000.0, 850.0, 500.0])
    assert an.temperature.shape == (3, 3, 4)
    assert an.u_wind.shape == (3, 3, 4)
    assert an.humidity.shape == (3, 3, 4)
    assert an.surface_pressure[0, 0] == pytest.approx(1013.25)
    assert an.sea_level_pressure[0, 0] == pytest.approx(1015.0)
    assert an.topography[0, 0] == pytest.approx(500.0)
    assert an.sst[0, 0] == pytest.approx(290.0)
    assert np.all(an.ice == 0.0)
    assert an.pressure_hpa.tolist() == [1000.0, 850.0, 500.0]


def test_build_analysis_missing_level_raises():
    fields = [f for f in _ncep_fields() if not (f.short_name == "t" and f.level == 850.0
                                                and f.level_type == "isobaricInhPa")]
    with pytest.raises(KeyError):
        build_analysis(fields, "nl", [1000.0, 850.0, 500.0])


def test_geopotential_to_height():
    assert geopotential_to_height(np.array([980.665]))[0] == pytest.approx(100.0)
