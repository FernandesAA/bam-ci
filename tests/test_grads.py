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


def _write_mixed_case(tmp_path: Path, values: np.ndarray):
    """Caso sintetico: 3-D (42 niveis) intercalado com duas variaveis 2-D."""
    levels = " ".join(str(1000 - 20 * i) for i in range(42))
    text = (
        "dset ^mixed.dat\n"
        "undef -9.99e8\n"
        "options big_endian yrev\n"
        "xdef 4 linear 0 90\n"
        "ydef 2 linear -45 90\n"
        "zdef 42 levels %s\n"
        "tdef 2 linear 00Z05NOV2026 3hr\n"
        "vars 4\n"
        "pslc 0 99 pressao\n"
        "tvir 42 %s temperatura\n"
        "prec 0 99 precipitacao\n"
        "uvel 42 %s vento\n"
        "endvars\n" % (levels, levels, levels)
    )
    ctl_path = tmp_path / "mixed.ctl"
    ctl_path.write_text(text, encoding="utf-8")
    (tmp_path / "mixed.dat").write_bytes(np.asarray(values, dtype=">f4").tobytes())
    return ctl_path


def test_variavel_2d_ocupa_um_nivel_e_nao_desalinha(tmp_path):
    nx, ny, nt, nlev3 = 4, 2, 2, 42
    cell = nx * ny
    # por tempo: pslc(1) + tvir(42) + prec(1) + uvel(42) = 86 registros
    records_per_time = 1 + nlev3 + 1 + nlev3
    values = np.arange(1, nt * records_per_time * cell + 1, dtype=float)
    ctl = _write_mixed_case(tmp_path, values)

    g = GradsFile(ctl)
    assert g._layout == {
        "pslc": (0, 1),
        "tvir": (1 * 4 * cell, 42),
        "prec": ((1 + nlev3) * 4 * cell, 1),
        "uvel": ((1 + nlev3 + 1) * 4 * cell, 42),
    }
    assert g.time_block_bytes == records_per_time * 4 * cell

    def expected(t: int, record: int) -> np.ndarray:
        start = ((t - 1) * records_per_time + record) * cell
        return values[start:start + cell].reshape(ny, nx)[::-1, :]

    two_d = g.read("pslc", t=1)
    assert two_d.shape == (ny, nx)          # (nlat, nlon)
    np.testing.assert_allclose(two_d, expected(1, 0))
    assert g.read_volume("pslc", t=1).shape == (1, ny, nx)

    np.testing.assert_allclose(g.read("tvir", t=1, z=1), expected(1, 1))
    np.testing.assert_allclose(g.read("tvir", t=1, z=42), expected(1, 42))
    np.testing.assert_allclose(g.read("prec", t=1), expected(1, 43))
    np.testing.assert_allclose(g.read("uvel", t=1, z=1), expected(1, 44))
    np.testing.assert_allclose(g.read("uvel", t=1, z=42), expected(1, 85))
    np.testing.assert_allclose(g.read("uvel", t=2, z=42), expected(2, 85))

    assert g.read_volume("uvel", t=1).shape == (nlev3, ny, nx)
    assert g.read_time(t=1)["prec"].shape == (1, ny, nx)


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
