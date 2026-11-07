import numpy as np

from bamci.compare import compare_ic, format_report
from bamci.ic_format import ICHeader, n_coeffs, write_ic


def _mk(mend=6, kmax=3, scale=1.0, seed=0):
    rng = np.random.default_rng(seed)
    n = n_coeffs(mend)
    return dict(topo=rng.normal(scale=scale, size=n),
                lnps=rng.normal(scale=scale, size=n),
                tvir=rng.normal(scale=scale, size=(n, kmax)),
                divg=rng.normal(scale=1e-6 * scale, size=(n, kmax)),
                vort=rng.normal(scale=1e-6 * scale, size=(n, kmax)),
                sphu=np.abs(rng.normal(scale=1e-4 * scale, size=(n, kmax))))


def test_compare_identical_is_zero(tmp_path):
    f = _mk()
    h = ICHeader(sigma_inter=np.linspace(1, 0, 4), sigma_layers=np.linspace(0.9, 0.1, 3))
    p1 = write_ic(tmp_path / "a.unf", h, **f)
    p2 = write_ic(tmp_path / "b.unf", h, **f)
    res = compare_ic(p1, p2, mend=6, nlat=8, nlon=12)
    assert res["compatible"]
    for m in res["fields"].values():
        assert m["coef_max"] == 0.0
        assert m["grid_rms"] == 0.0


def test_compare_detects_difference(tmp_path):
    f1 = _mk()
    f2 = _mk(scale=2.0, seed=1)
    h = ICHeader(sigma_inter=np.linspace(1, 0, 4), sigma_layers=np.linspace(0.9, 0.1, 3))
    p1 = write_ic(tmp_path / "a.unf", h, **f1)
    p2 = write_ic(tmp_path / "b.unf", h, **f2)
    res = compare_ic(p1, p2, mend=6, nlat=8, nlon=12)
    assert res["fields"]["topo"]["coef_max"] > 0
    assert "topo" in format_report(res)


def test_compare_incompatible(tmp_path):
    h = ICHeader(sigma_inter=np.linspace(1, 0, 4), sigma_layers=np.linspace(0.9, 0.1, 3))
    p1 = write_ic(tmp_path / "a.unf", h, **_mk(mend=6))
    p2 = write_ic(tmp_path / "b.unf", h, **_mk(mend=7, kmax=3))
    res = compare_ic(p1, p2, mend=6, nlat=8, nlon=12)
    assert res["compatible"] is False
