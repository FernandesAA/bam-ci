import struct
from pathlib import Path

import numpy as np
import pytest

from bamci.ic_format import (ICHeader, check_ic, coeff_bytes, ic_file_size,
                             n_coeffs, pair_index, read_ic, write_ic)


def _mk(mend=6, kmax=3, seed=0):
    rng = np.random.default_rng(seed)
    n = n_coeffs(mend)
    return dict(
        topo=rng.normal(size=n),
        lnps=rng.normal(size=n),
        tvir=rng.normal(size=(n, kmax)),
        divg=rng.normal(size=(n, kmax)) * 1e-5,
        vort=rng.normal(size=(n, kmax)) * 1e-5,
        sphu=np.abs(rng.normal(size=(n, kmax))) * 1e-3,
    )


def test_pair_index_file_order():
    idx = pair_index(4)
    # m externo: (0,0),(0,1),(0,2),(0,3),(0,4),(1,1),(1,2)...
    assert [idx[0, n] for n in range(5)] == [0, 1, 2, 3, 4]
    assert idx[1, 1] == 5 and idx[1, 4] == 8
    assert idx[2, 2] == 9
    assert idx[2, 1] == -1           # fora do triangulo
    assert sorted(idx[idx >= 0].ravel().tolist()) == list(range(15))


def test_write_read_roundtrip(tmp_path):
    fields = _mk()
    header = ICHeader(date_inp=(0, 10, 1, 2026), date_cur=(0, 10, 1, 2026),
                      a_hybr=np.zeros(4), b_hybr=np.linspace(1, 0, 4))
    path = write_ic(tmp_path / "ic.unf", header, **fields)
    assert path.stat().st_size == ic_file_size(6, 3, "hybrid")
    back = read_ic(path)
    np.testing.assert_allclose(back["topo"], fields["topo"], rtol=0, atol=1e-5)
    np.testing.assert_allclose(back["tvir"], fields["tvir"], rtol=0, atol=1e-5)
    assert back["kmax"] == 3
    assert back["date_inp"] == (0, 10, 1, 2026)


def test_sigma_header_layout(tmp_path):
    header = ICHeader(sigma_inter=np.array([0.0, 0.3, 0.6, 1.0]),
                      sigma_layers=np.array([0.15, 0.45, 0.8]))
    fields = _mk()
    path = write_ic(tmp_path / "ic_sigma.unf", header, **fields)
    assert path.stat().st_size == ic_file_size(6, 3, "sigma")
    back = read_ic(path)
    np.testing.assert_allclose(back["a_hybr"], [0.0, 0.3, 0.6, 1.0])


def test_size_matches_production():
    assert ic_file_size(126, 42) == 11_055_828


def test_check_ic_detects_imaginary_zonal(tmp_path):
    fields = _mk(mend=4, kmax=2)
    fields["topo"] = fields["topo"].copy()
    idx = pair_index(4)
    fields["topo"][2 * idx[0, 1] + 1] = 5.0        # seno de m=0 proibido (n=1)
    header = ICHeader(a_hybr=np.zeros(3), b_hybr=np.linspace(1, 0, 3))
    path = write_ic(tmp_path / "bad.unf", header, **fields)
    problems = check_ic(path, mend=4)
    assert any("seno de m=0" in p for p in problems)
