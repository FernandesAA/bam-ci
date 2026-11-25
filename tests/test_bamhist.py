"""Testes do leitor do historico do BAM (.dir + .fct) - bamci.bamhist.

Os casos sao SINTETICOS: os testes escrevem um .dir/.fct no formato inferido do
legado (markers de 4 bytes, cabecalho de 10 int32, uma variavel/nivel por
registro) e conferem a leitura contra a sintese espectral do proprio pacote.
"""

from __future__ import annotations

import struct
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest

from bamci import bamhist, mb09
from bamci.ctl import read_ctl
from bamci.derivatives import wind_from_div_vort
from bamci.grads import GradsFile
from bamci.ic_format import n_coeffs, pair_index
from bamci.levels import gaussian_latitudes, grid_shape
from bamci.spectral import file_order, from_file_order, synthesis

HEADER_SECONDS = 3600.0


# ---------------------------------------------------------------------------
# gerador do caso sintetico
# ---------------------------------------------------------------------------
def _table(mend, kmax, nlat, nlon):
    ncoeff = n_coeffs(mend)
    npts = nlat * nlon
    return [
        ("TOPOGRAPHY", "FIXD", ncoeff, 1, 10),
        ("LN SURFACE PRESSURE", "PROG", ncoeff, 1, 142),
        ("VIRTUAL TEMPERATURE", "PROG", ncoeff, kmax, 40),
        ("DIVERGENCE", "PROG", ncoeff, kmax, 50),
        ("VORTICITY", "PROG", ncoeff, kmax, 50),
        ("SPECIFIC HUMIDITY", "PROG", ncoeff, kmax, 0),
        ("LAND SEA MASK", "FIXD", npts, 1, 0),
    ]


def _coefficients(mend, kmax, seed=0, scale=1.0):
    """Coeficientes (ncoeff,) ou (kmax, ncoeff) com seno de m=0 nulo."""
    rng = np.random.default_rng(seed)
    ncoeff = n_coeffs(mend)
    idx = pair_index(mend)
    zonal = idx[0] >= 0
    def one():
        v = rng.normal(size=ncoeff) * scale
        v[2 * idx[0][zonal] + 1] = 0.0        # seno de m=0 = 0
        return v
    return one(), np.stack([one() for _ in range(kmax)])


def _write_case(root, *, stem="GFCTTEST", mend=3, kmax=4, lead_hours=3.0,
                seed=0, analysis=(0, 11, 19, 2026), data=None):
    """Escreve <stem>.dir + <stem>.fct sinteticos; devolve (dir, fct)."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    nlon, nlat = grid_shape(mend)
    table = _table(mend, kmax, nlat, nlon)
    delta = np.full(kmax, 1.0 / kmax)
    if data is None:
        topo, _ = _coefficients(mend, 1, seed=seed, scale=10.0)
        lnps = np.zeros(n_coeffs(mend))
        lnps[0] = np.sqrt(2.0) * np.log(100.0)
        _, tvir = _coefficients(mend, kmax, seed=seed + 1, scale=30.0)
        tvir[0] = 0.0
        tvir[0, 0] = np.sqrt(2.0) * 260.0
        _, divg = _coefficients(mend, kmax, seed=seed + 2, scale=1e-5)
        _, vort = _coefficients(mend, kmax, seed=seed + 3, scale=3e-5)
        _, sphu = _coefficients(mend, kmax, seed=seed + 4, scale=1e-4)
        sphu[0] = 0.0
        sphu[0, 0] = np.sqrt(2.0) * 2e-3
        lsm = np.ones(nlat * nlon)
        data = {"TOPOGRAPHY": topo, "LN SURFACE PRESSURE": lnps,
                "VIRTUAL TEMPERATURE": tvir, "DIVERGENCE": divg,
                "VORTICITY": vort, "SPECIFIC HUMIDITY": sphu,
                "LAND SEA MASK": lsm}

    when = datetime(*[int(v) for v in (analysis[3], analysis[1], analysis[2],
                                       analysis[0])]) + timedelta(hours=lead_hours)
    hh_i, mm_i, dd_i, yy_i = analysis
    head2 = ("0001 SEQU T%04d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d TAPE"
             % (mend, mend + 1, kmax, kmax, when.hour, when.day, when.month,
                when.year, hh_i, mm_i, dd_i, yy_i))
    lines = ["CPTEC SIGMA VERS 2.0", head2,
             "CPTEC AGCM REVIS 1.0 2000  T%04d%02d   COLD" % (mend, kmax)]
    for start in range(0, kmax, 5):
        lines.append("  " + "  ".join("%.8E" % v for v in delta[start:start + 5]))
    for name, kind, npts, nlev, code in table:
        lines.append("%-40s %-4s %8d %6d %6d" % (name, kind, npts, nlev, code))
    dir_path = root / (stem + ".dir.TQ")
    dir_path.write_text("\n".join(lines) + "\n", encoding="ascii")

    def rec(payload: bytes) -> bytes:
        n = struct.pack(">i", len(payload))
        return n + payload + n

    header = struct.pack(">ii4i4i", int(lead_hours // 24),
                         int(round((lead_hours % 24) * HEADER_SECONDS)),
                         hh_i, mm_i, dd_i, yy_i,
                         when.hour, when.month, when.day, when.year)
    chunks = [rec(header)]
    for name, kind, npts, nlev, code in table:
        arr = np.asarray(data[name], dtype=np.float64)
        if arr.ndim == 1:
            arr = arr[None]
        chunks.extend(rec(np.asarray(row, dtype=">f4").tobytes()) for row in arr)
    fct_path = root / (stem + ".fct.TQ")
    fct_path.write_bytes(b"".join(chunks))
    return dir_path, fct_path


# ---------------------------------------------------------------------------
# .dir
# ---------------------------------------------------------------------------
def test_parse_dir_header_sigma_and_table(tmp_path):
    dir_path, _ = _write_case(tmp_path, mend=3, kmax=4, lead_hours=9.0)
    info = bamhist.parse_dir(dir_path)
    assert info.version.startswith("CPTEC SIGMA")
    assert info.header.trunc == "T0003"
    assert info.header.mwv == 4 and info.mend == 3
    assert info.header.nz == 4 and info.kmax == 4
    assert len(info.variables) == 7
    assert info.names[0] == "TOPOGRAPHY" and info.names[-1] == "LAND SEA MASK"
    var = info.variables[2]
    assert (var.name, var.kind, var.npts, var.nlev, var.code) == \
        ("VIRTUAL TEMPERATURE", "PROG", n_coeffs(3), 4, 40)
    assert var.index == 2 and var.nbytes == 4 * n_coeffs(3)
    np.testing.assert_allclose(info.sigma.sum(), 1.0)
    assert info.sigma_interfaces[0] == 1.0 and info.sigma_interfaces[-1] == 0.0
    assert np.all(np.diff(info.sigma_interfaces) < 0)
    # datas: linha 2 leva a data corrente em (hh, dd, mm, yyyy) e a inicial em
    # (hh, mm, dd, yyyy); o .fct e' a referencia
    assert info.header.current_time == datetime(2026, 11, 19, 9)
    assert info.header.initial_time == datetime(2026, 11, 19, 0)


def test_dir_aliases_resolve():
    # nomes normalizados -> linhas da tabela
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        dir_path, _ = _write_case(tmp)
        info = bamhist.parse_dir(dir_path)
        assert info.variable("PS").name == "LN SURFACE PRESSURE"
        assert info.variable("T").name == "VIRTUAL TEMPERATURE"
        assert info.variable("Q").name == "SPECIFIC HUMIDITY"
        assert info.variable("DIV").name == "DIVERGENCE"
        assert info.variable("VORT").name == "VORTICITY"
        assert info.variable("TOPO").name == "TOPOGRAPHY"
        with pytest.raises(KeyError):
            info.variable("NAO EXISTE")
        assert info.to_dict()["version"].startswith("CPTEC")


def _dir_text(*, version="CPTEC SIGMA VERS 2.0", kmax=4, delta=None,
              table=True):
    delta = np.full(kmax, 1.0 / kmax) if delta is None else np.asarray(delta)
    lines = [version,
             "0001 SEQU T0003    4    4    4    3   19   11 2026    0   11   19 2026 TAPE",
             "CPTEC AGCM REVIS 1.0 2000  T000342   COLD"]
    for start in range(0, len(delta), 5):
        lines.append("  " + "  ".join("%.8E" % v for v in delta[start:start + 5]))
    if table:
        lines += [
            "TOPOGRAPHY                                FIXD     20   1  10",
            "LN SURFACE PRESSURE                       PROG     20   1 142",
            "VIRTUAL TEMPERATURE                       PROG     20   4  40",
            "LAND SEA MASK                             FIXD     72   1   0",
        ]
    return "\n".join(lines) + "\n"


def test_parse_dir_rejects_bad_files(tmp_path):
    bad = tmp_path / "bad1.dir"
    bad.write_text(_dir_text(version="OUTRA COISA"))
    with pytest.raises(bamhist.FctError, match="primeira linha"):
        bamhist.parse_dir(bad)

    bad2 = tmp_path / "bad2.dir"
    bad2.write_text(_dir_text(delta=np.full(4, 0.5)))       # soma 2
    with pytest.raises(bamhist.FctError, match="sigma nao soma 1"):
        bamhist.parse_dir(bad2)

    bad3 = tmp_path / "bad3.dir"
    bad3.write_text(_dir_text(table=False))
    with pytest.raises(bamhist.FctError, match="tabela de variaveis vazia"):
        bamhist.parse_dir(bad3)

    bad4 = tmp_path / "bad4.dir"
    bad4.write_text(_dir_text(delta=np.full(3, 1.0 / 3.0)))  # 3 valores, nz=4
    with pytest.raises(bamhist.FctError, match="valores de sigma"):
        bamhist.parse_dir(bad4)


# ---------------------------------------------------------------------------
# .fct
# ---------------------------------------------------------------------------
def test_fct_size_and_verify(tmp_path):
    dir_path, fct_path = _write_case(tmp_path, mend=3, kmax=4)
    info = bamhist.parse_dir(dir_path)
    expected = 8 + 40 + sum(v.nlev * (8 + v.nbytes) for v in info.variables)
    assert fct_path.stat().st_size == expected == bamhist.fct_size(info)
    fct = bamhist.read_fct(fct_path)
    assert fct.verify() == []
    assert fct.dir.path == dir_path
    assert [v.name for v in fct.variables] == list(info.names)
    assert fct.mend == 3 and (fct.nlat, fct.nlon) == grid_shape(3)[::-1]


def test_fct_header_dates(tmp_path):
    _, fct_path = _write_case(tmp_path, lead_hours=9.0)
    fct = bamhist.read_fct(fct_path)
    head = fct.header
    assert head.forecast_day == 0 and head.seconds == 9 * 3600
    assert head.date_inp == (0, 11, 19, 2026)
    assert head.date_cur == (9, 11, 19, 2026)
    assert head.lead_hours == 9.0
    assert head.analysis_time == datetime(2026, 11, 19, 0)
    assert head.valid_time == datetime(2026, 11, 19, 9)
    # o PostGrib le esse campo como REAL4 (bug do legado) -> valor absurdo
    assert head.time_of_day_real != 9 * 3600
    assert abs(head.time_of_day_real) < 1e-6


def test_fct_verify_detects_corruption(tmp_path):
    _, fct_path = _write_case(tmp_path)
    raw = bytearray(fct_path.read_bytes())
    raw[48:52] = struct.pack(">i", 999)              # 1o marcador de dado
    fct_path.write_bytes(bytes(raw))
    with pytest.raises(bamhist.FctError, match="marcador"):
        bamhist.read_fct(fct_path)
    fct = bamhist.read_fct(fct_path, check=False)
    problems = fct.verify()
    assert any("marcador" in p for p in problems)
    # e um arquivo truncado e' detectado pelo tamanho
    fct_path.write_bytes(bytes(raw)[:-8])
    assert any("tamanho" in p
               for p in bamhist.read_fct(fct_path, check=False).verify())


def test_fct_read_spectral_and_grid(tmp_path):
    data = None
    dir_path, fct_path = _write_case(tmp_path, kmax=3, data=data)
    fct = bamhist.read_fct(fct_path)
    topo = fct.spectral("TOPO")
    assert topo.shape == (1, n_coeffs(3))
    assert fct.is_spectral("T") and not fct.is_spectral("LSM")
    lsm = fct.grid("LSM")
    assert lsm.shape == (1,) + tuple(grid_shape(3)[::-1])
    np.testing.assert_allclose(lsm[0], 1.0)
    with pytest.raises(bamhist.FctError, match="pontos de grade"):
        fct.grid("T")
    with pytest.raises(bamhist.FctError, match="coeficientes espectrais"):
        fct.spectral("LSM")
    with pytest.raises(IndexError):
        fct.read("T", 99)


def test_fct_read_returns_float32_values(tmp_path):
    rng = np.random.default_rng(7)
    mend, kmax = 3, 2
    topo, others = _coefficients(mend, kmax, seed=1)
    data = {"TOPOGRAPHY": topo, "LN SURFACE PRESSURE": np.zeros(n_coeffs(mend)),
            "VIRTUAL TEMPERATURE": others, "DIVERGENCE": others * 0,
            "VORTICITY": others * 0, "SPECIFIC HUMIDITY": others * 0,
            "LAND SEA MASK": np.zeros(grid_shape(mend)[0] * grid_shape(mend)[1])}
    _, fct_path = _write_case(tmp_path, mend=mend, kmax=kmax, data=data)
    fct = bamhist.read_fct(fct_path)
    np.testing.assert_allclose(fct.spectral("TOPO")[0], topo, rtol=0, atol=1e-6)
    np.testing.assert_allclose(fct.spectral("T"), others, rtol=0, atol=1e-6)


def test_synth_stack_matches_synthesis(tmp_path):
    _, fct_path = _write_case(tmp_path, kmax=3)
    fct = bamhist.read_fct(fct_path)
    vec = fct.spectral("T")
    stack = bamhist._synth_stack(bamhist._split_pairs(vec, fct.mend),
                                 fct.mend, fct.nlat, fct.nlon)
    nlon, nlat = fct.nlon, fct.nlat
    for k in range(vec.shape[0]):
        ref = synthesis(from_file_order(vec[k], fct.mend), fct.mend, nlat, nlon)
        np.testing.assert_allclose(stack[k], ref, rtol=0, atol=1e-10)


def test_split_pairs_uses_ic_order():
    mend = 4
    idx = pair_index(mend)
    vec = np.arange(2 * n_coeffs(mend), dtype=float)
    coef = bamhist._split_pairs(vec[None], mend)[0]
    for m in range(mend + 1):
        for n in range(m, mend + 1):
            assert coef[m, n] == complex(vec[2 * idx[m, n]], vec[2 * idx[m, n] + 1])
    # e a ida-e-volta com a ordem do IC
    np.testing.assert_allclose(
        bamhist._split_pairs(file_order(coef)[None], mend)[0], coef, rtol=0, atol=0)


def test_fields_physical_meaning(tmp_path):
    mend, kmax = 3, 4
    nlon, nlat = grid_shape(mend)
    ncoeff = n_coeffs(mend)
    lnps = np.zeros(ncoeff)
    lnps[0] = np.sqrt(2.0) * np.log(100.0)          # ps = 1000 hPa
    tvir = np.zeros((kmax, ncoeff))
    tvir[:, 0] = np.sqrt(2.0) * 250.0               # 250 K em todos os niveis
    zeros = np.zeros((kmax, ncoeff))
    data = {"TOPOGRAPHY": np.zeros(ncoeff), "LN SURFACE PRESSURE": lnps,
            "VIRTUAL TEMPERATURE": tvir, "DIVERGENCE": zeros, "VORTICITY": zeros,
            "SPECIFIC HUMIDITY": zeros,
            "LAND SEA MASK": np.zeros(nlon * nlat)}
    _, fct_path = _write_case(tmp_path, mend=mend, kmax=kmax, data=data)
    fld = bamhist.read_fct(fct_path).fields()
    # o arquivo guarda float32: 1e-7 de erro relativo em lnps (e em ps)
    np.testing.assert_allclose(fld["ps"], 1000.0, rtol=1e-6)
    np.testing.assert_allclose(fld["t"], 250.0, rtol=1e-6)
    np.testing.assert_allclose(fld["q"], 0.0, atol=1e-20)
    np.testing.assert_allclose(fld["u"], 0.0, atol=1e-12)
    np.testing.assert_allclose(fld["topo"], 0.0, atol=1e-12)
    np.testing.assert_allclose(fld["lat"], gaussian_latitudes(nlat))
    assert fld["lat"][0] > fld["lat"][-1]           # norte -> sul
    assert fld["sigma_layers"].shape == (kmax,)


def test_fields_clip_q(tmp_path):
    mend, kmax = 3, 2
    ncoeff = n_coeffs(mend)
    q = np.zeros((kmax, ncoeff))
    q[:, 0] = -np.sqrt(2.0) * 1e-3                  # q negativo
    zeros = np.zeros((kmax, ncoeff))
    data = {"TOPOGRAPHY": np.zeros(ncoeff),
            "LN SURFACE PRESSURE": np.zeros(ncoeff),
            "VIRTUAL TEMPERATURE": zeros + 0.0, "DIVERGENCE": zeros,
            "VORTICITY": zeros, "SPECIFIC HUMIDITY": q,
            "LAND SEA MASK": np.zeros(grid_shape(mend)[0] * grid_shape(mend)[1])}
    _, fct_path = _write_case(tmp_path, mend=mend, kmax=kmax, data=data)
    fct = bamhist.read_fct(fct_path)
    assert float(fct.fields(clip_q=True)["q"].min()) >= 0.0
    assert float(fct.fields(clip_q=False)["q"].min()) < 0.0


def test_fct_files_sorted_by_lead(tmp_path):
    for lead in (9.0, 3.0, 6.0):
        _write_case(tmp_path, stem="GFCT%02d" % int(lead), lead_hours=lead)
    paths = bamhist.fct_files(tmp_path)
    assert [bamhist.read_fct(p).lead_hours for p in paths] == [3.0, 6.0, 9.0]
    assert len(bamhist.read_fct_series(tmp_path)) == 3
    with pytest.raises(FileNotFoundError):
        bamhist.fct_files(tmp_path / "nao_existe")


# ---------------------------------------------------------------------------
# .ctl equivalente
# ---------------------------------------------------------------------------
def test_write_grads_roundtrip(tmp_path):
    mend, kmax = 3, 3
    nlon, nlat = grid_shape(mend)
    lat = gaussian_latitudes(nlat)
    lon = np.arange(nlon) * 360.0 / nlon
    rng = np.random.default_rng(3)
    ps = 950.0 + rng.normal(size=(nlat, nlon))
    t = 280.0 + rng.normal(size=(kmax, nlat, nlon))
    q = 1e-2 + 1e-4 * rng.random((kmax, nlat, nlon))
    u = 5.0 * rng.normal(size=(kmax, nlat, nlon))
    v = 5.0 * rng.normal(size=(kmax, nlat, nlon))
    fields = {"PS": ps, "T": t, "Q": q, "U": u, "V": v}
    levels = np.linspace(0.9, 0.1, kmax)
    dat, ctl = bamhist.write_grads(fields, tmp_path / "hist", lat=lat, lon=lon,
                                  levels=levels,
                                  time=datetime(2026, 11, 19, 3))
    assert dat.stat().st_size == 4 * nlon * nlat * (1 + 4 * kmax)
    c = read_ctl(ctl)
    assert (c.nx, c.ny, c.nz, c.nt) == (nlon, nlat, kmax, 1)
    assert "yrev" in c.options and "big_endian" in c.options
    assert [v.name for v in c.variables] == ["pslc", "tvir", "umes", "uvel", "vvel"]
    assert c.tdef.times().__next__() == datetime(2026, 11, 19, 3)
    g = GradsFile(c)
    # o .dat vai de norte para sul e o .ctl tem yrev: GradsFile devolve sul->norte
    np.testing.assert_allclose(g.read("pslc")[::-1], ps, rtol=0, atol=1e-4)
    np.testing.assert_allclose(g.read("tvir", z=2)[::-1], t[1], rtol=0, atol=1e-4)
    np.testing.assert_allclose(g.read("uvel", z=kmax)[::-1], u[-1], rtol=0, atol=1e-6)


def test_fct_to_grads(tmp_path):
    _, fct_path = _write_case(tmp_path, kmax=3)
    fct = bamhist.read_fct(fct_path)
    dat, ctl = bamhist.fct_to_grads(fct, tmp_path / "h")
    assert dat.stat().st_size == 4 * fct.nlon * fct.nlat * (1 + 4 * 3)
    c = read_ctl(ctl)
    assert c.nz == 3
    np.testing.assert_allclose(c.zdef.values(), fct.dir.sigma_layers, rtol=0, atol=1e-6)
    g = GradsFile(c)
    ps = g.read("pslc")[::-1]
    np.testing.assert_allclose(ps, fct.fields(keys=("PS",))["ps"], rtol=0, atol=1e-3)


def test_write_grads_rejects_bad_shapes(tmp_path):
    lat = np.array([30.0, -30.0])
    lon = np.array([0.0, 180.0])
    with pytest.raises(ValueError, match="forma"):
        bamhist.write_grads({"T": np.zeros((2, 3, 3))}, tmp_path / "x", lat=lat, lon=lon)
    with pytest.raises(ValueError, match="niveis"):
        bamhist.write_grads({"T": np.zeros((2, 2, 2))}, tmp_path / "x", lat=lat,
                            lon=lon, levels=[1.0])


# ---------------------------------------------------------------------------
# integracao com o MB09 (--source bamfct)
# ---------------------------------------------------------------------------
def _series_case(root, *, hours, seed=0, drop=None, name="m"):
    """Uma 'rodada': um .fct por hora (lead em horas)."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    for hour in hours:
        if drop is not None and float(hour) in tuple(drop):
            continue
        _write_case(root, stem="GFCT%s%02d" % (name, int(hour)),
                    lead_hours=float(hour), seed=seed + int(hour))
    return root


def test_build_eofs_from_fct(tmp_path):
    ctrl = _series_case(tmp_path / "ctrl", hours=(3, 6, 9), seed=0, name="c")
    memb = _series_case(tmp_path / "m01", hours=(3, 6, 9), seed=100, name="m")
    eofs = mb09.build_eofs_from_fct(ctrl, [memb], n=1, window=(3, 9, 3),
                                    remove_mean=True)
    assert eofs["m_times"] == 3 and eofs["n"] == 1
    assert eofs["meta"]["source"] == "bamfct"
    assert eofs["meta"]["hours"] == [3.0, 6.0, 9.0]
    assert eofs["meta"]["n_ignored_times"] == 0
    assert eofs["meta"]["analysis_time"] == "2026111900"
    assert set(eofs["eof_fields"]) == {"P", "T", "Q", "wind"}
    assert (eofs["nlat"], eofs["nlon"]) == grid_shape(3)[::-1]


def test_build_eofs_from_fct_ignores_extra_control_times(tmp_path):
    ctrl = _series_case(tmp_path / "ctrl", hours=(3, 6, 9, 12), seed=0, name="c")
    memb = _series_case(tmp_path / "m01", hours=(3, 6, 9), seed=100, name="m")
    eofs = mb09.build_eofs_from_fct(ctrl, [memb], n=1, window=(3, 9, 3))
    assert eofs["m_times"] == 3
    assert eofs["meta"]["n_ignored_times"] == 1


def test_build_eofs_from_fct_missing_hour(tmp_path):
    ctrl = _series_case(tmp_path / "ctrl", hours=(3, 6, 9), seed=0, name="c")
    memb = _series_case(tmp_path / "m01", hours=(3, 6, 9), seed=100, name="m",
                        drop=(6,))
    with pytest.raises(ValueError, match="faltam os tempos 6"):
        mb09.build_eofs_from_fct(ctrl, [memb], n=1, window=(3, 9, 3))


def test_build_eofs_from_fct_requires_odd_n(tmp_path):
    ctrl = _series_case(tmp_path / "ctrl", hours=(3, 6, 9), seed=0, name="c")
    memb = _series_case(tmp_path / "m01", hours=(3, 6, 9), seed=100, name="m")
    with pytest.raises(ValueError, match="impar"):
        mb09.build_eofs_from_fct(ctrl, [memb], n=2, window=(3, 9, 3))


def test_cli_build_bamfct(tmp_path, capsys):
    ctrl = _series_case(tmp_path / "ctrl", hours=(3, 6, 9), seed=0, name="c")
    memb = _series_case(tmp_path / "m01", hours=(3, 6, 9), seed=100, name="m")
    out = tmp_path / "eof.npz"
    code = mb09.main(["build", "--source", "bamfct", "--control", str(ctrl),
                      "--members", str(memb), "--window", "3:9:3", "--n", "1",
                      "--remove-mean", "no", "--out", str(out)])
    assert code == 0
    assert out.exists()
    txt = capsys.readouterr().out
    assert "horas=3,6,9" in txt
    eofs = mb09.load_eofs(out)
    assert eofs["meta"]["source"] == "bamfct"
    assert eofs["m_times"] == 3


def test_cli_build_bamfct_discovers_work_layout(tmp_path):
    root = tmp_path / "work" / "ens126"
    _series_case(root / "ctrl" / "out", hours=(3, 6, 9), seed=0, name="c")
    _series_case(root / "m01" / "out", hours=(3, 6, 9), seed=100, name="m")
    out = tmp_path / "eof.npz"
    code = mb09.main(["build", "--source", "bamfct", "--work", str(tmp_path / "work"),
                      "--window", "3:9:3", "--n", "1", "--out", str(out)])
    assert code == 0 and out.exists()
    eofs = mb09.load_eofs(out)
    assert "ens126" in eofs["meta"]["control"]
    assert len(eofs["meta"]["members"]) == 1


# ---------------------------------------------------------------------------
# regressao do nucleo espectral compartilhado
# ---------------------------------------------------------------------------
def test_wind_from_div_vort_stacked_equals_per_level():
    rng = np.random.default_rng(0)
    mend, nlat, nlon, nlev = 8, 12, 24, 3
    div = (rng.normal(size=(nlev, mend + 1, mend + 1))
           + 1j * rng.normal(size=(nlev, mend + 1, mend + 1))) * 1e-5
    vort = div * 2.0
    u, v = wind_from_div_vort(div, vort, mend, nlat, nlon)
    assert u.shape == (nlev, nlat, nlon)
    for k in range(nlev):
        u1, v1 = wind_from_div_vort(div[k], vort[k], mend, nlat, nlon)
        np.testing.assert_allclose(u[k], u1, rtol=0, atol=1e-10)
        np.testing.assert_allclose(v[k], v1, rtol=0, atol=1e-10)
