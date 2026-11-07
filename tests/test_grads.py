from pathlib import Path

import numpy as np
import pytest

from bamci.grads import GradsFile

CTL = """dset ^analise.dat
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


def _write_case(tmp_path: Path, values: np.ndarray) -> Path:
    ctl_path = tmp_path / "analise.ctl"
    ctl_path.write_text(CTL, encoding="utf-8")
    (tmp_path / "analise.dat").write_bytes(
        np.asarray(values, dtype=">f4").tobytes()
    )
    return ctl_path


def test_read_field_and_yrev(tmp_path):
    # layout: t -> var -> nivel -> y -> x  (2*2*6 = 24 valores por tempo)
    t1 = np.arange(1, 13, dtype=float)   # u lev1(6), u lev2(6)
    t1 = np.concatenate([t1, np.arange(21, 33, dtype=float)])  # t lev1, t lev2
    t2 = t1 + 100.0
    ctl = _write_case(tmp_path, np.concatenate([t1, t2]))

    g = GradsFile(ctl)
    raw_u1 = np.arange(1, 7, dtype=float).reshape(2, 3)  # sul -> norte
    got = g.read("u", t=1, z=1)
    np.testing.assert_allclose(got, raw_u1[::-1, :])      # yrev aplicado
    np.testing.assert_allclose(g.read("u", t=1, z=2), np.arange(7, 13).reshape(2, 3)[::-1, :])
    np.testing.assert_allclose(g.read("t", t=2, z=1), (np.arange(21, 27) + 100).reshape(2, 3)[::-1, :])
    assert g.read_volume("u", t=1).shape == (2, 2, 3)


def test_undef_becomes_nan(tmp_path):
    vals = np.arange(1, 25, dtype=float)
    vals[2] = -9.99e8
    ctl = _write_case(tmp_path, vals)
    field = GradsFile(ctl).read("u", t=1, z=1)
    assert np.isnan(field).sum() == 1


def test_truncated_file_raises(tmp_path):
    ctl = _write_case(tmp_path, np.arange(1, 25, dtype=float)[:10])
    with pytest.raises(ValueError):
        GradsFile(ctl).read("t", t=2, z=2)


def test_template_files(tmp_path):
    ctl_text = CTL.replace("dset ^analise.dat", "dset ^analise.%y4%m2%d2.dat").replace(
        "options big_endian yrev", "options big_endian yrev template"
    )
    ctl_path = tmp_path / "tpl.ctl"
    ctl_path.write_text(ctl_text, encoding="utf-8")
    block = np.arange(1, 25, dtype=float)
    (tmp_path / "analise.20200101.dat").write_bytes(np.asarray(block, dtype=">f4").tobytes())
    (tmp_path / "analise.20200102.dat").write_bytes(np.asarray(block + 100, dtype=">f4").tobytes())
    g = GradsFile(ctl_path)
    np.testing.assert_allclose(g.read("u", t=2, z=1), (np.arange(1, 7) + 100).reshape(2, 3)[::-1, :])
