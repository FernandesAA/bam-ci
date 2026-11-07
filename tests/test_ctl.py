import numpy as np
import pytest

from bamci.ctl import CtlError, parse_ctl

CTL = """dset ^analise.dat
title Analise BAM
undef -9.99e8
options big_endian yrev
xdef 3 linear 0 1
ydef 2 linear -10 20
zdef 2 levels 1000 500
tdef 2 linear 00Z01JAN2020 1dy
vars 2
u 2 1000 500 m/s vento zonal
t 2 1000 500 k temperatura
endvars
"""


def test_axes_and_vars():
    ctl = parse_ctl(CTL)
    assert (ctl.nx, ctl.ny, ctl.nz, ctl.nt) == (3, 2, 2, 2)
    np.testing.assert_allclose(ctl.xdef.values(), [0, 1, 2])
    np.testing.assert_allclose(ctl.ydef.values(), [-10, 10])
    np.testing.assert_allclose(ctl.zdef.values(), [1000, 500])
    assert [v.name for v in ctl.variables] == ["u", "t"]
    assert ctl.variable("u").levels == (1000.0, 500.0)
    assert "yrev" in ctl.options


def test_time_axis_days():
    ctl = parse_ctl(CTL)
    times = list(ctl.tdef.times())
    assert [t.day for t in times] == [1, 2]
    assert times[0].year == 2020


def test_grid_flips_lat_with_yrev():
    ctl = parse_ctl(CTL)
    grid = ctl.grid()
    np.testing.assert_allclose(grid["lat"], [10, -10])
    assert grid["time"].shape == (2,)


def test_rejects_unknown_key():
    with pytest.raises(CtlError):
        parse_ctl(CTL.replace("title Analise BAM", "titulo errado"))


def test_rejects_var_count_mismatch():
    with pytest.raises(CtlError):
        parse_ctl(CTL.replace("vars 2", "vars 3"))
