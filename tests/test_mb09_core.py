"""Testes do nucleo do MB09 (bamci.mb09), com estados sinteticos.

Nada aqui depende de arquivo real: os Estados sao construidos em memoria e os
poucos ICs usados (roundtrip/escrita) sao gerados com bamci.ic_format.write_ic.
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest

from bamci import mb09
from bamci import mb09_tables as mb
from bamci.ic_format import ICHeader, n_coeffs, write_ic


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _grid(nlat=10, nlon=20, north_to_south=True):
    if north_to_south:
        lat = np.linspace(90.0, -90.0, nlat)
    else:
        lat = np.linspace(-90.0, 90.0, nlat)
    lon = np.arange(nlon, dtype=float) * 360.0 / nlon
    return lat, lon


def _mask(lat, lon, region):
    # region_mask recebe as coordenadas em forma broadcastable (nlat, 1)/(1, nlon)
    return mb.region_mask(np.asarray(lat)[:, None], np.asarray(lon)[None, :], region)


def _state(nlat=10, nlon=20, nlev=3, seed=0, q_value=1e-2):
    rng = np.random.default_rng(seed)
    lat, lon = _grid(nlat, nlon)
    return mb09.State(
        t=280.0 + rng.normal(0.0, 2.0, (nlev, nlat, nlon)),
        q=q_value + 0.1 * q_value * rng.random((nlev, nlat, nlon)),
        u=rng.normal(0.0, 5.0, (nlev, nlat, nlon)),
        v=rng.normal(0.0, 5.0, (nlev, nlat, nlon)),
        ps=1000.0 + rng.normal(0.0, 3.0, (nlat, nlon)),
        lat=lat, lon=lon,
    )


def _make_ic(tmp_path: Path, mend=6, kmax=3, style="hybrid", seed=0, name="ic.unf"):
    rng = np.random.default_rng(seed)
    n = n_coeffs(mend)
    dates = dict(date_inp=(0, 10, 1, 2026), date_cur=(0, 10, 1, 2026))
    if style == "hybrid":
        header = ICHeader(a_hybr=np.zeros(kmax + 1),
                          b_hybr=np.linspace(1.0, 0.0, kmax + 1), **dates)
    else:
        header = ICHeader(sigma_inter=np.linspace(1.0, 0.0, kmax + 1),
                          sigma_layers=0.5 * (np.linspace(1.0, 0.0, kmax + 1)[:-1]
                                              + np.linspace(1.0, 0.0, kmax + 1)[1:]),
                          **dates)
    # o coeficiente (0,0) carrega a media (media = c00/sqrt(2)); mante-lo em
    # ln(100) e 280 K da ps ~1000 hPa e Tv plausivel. Os demais coeficientes sao
    # apenas ruido pequeno.
    lnps_c = 0.02 * rng.normal(size=n)
    lnps_c[0] = np.sqrt(2.0) * np.log(100.0)
    tvir_c = 2.0 * rng.normal(size=(n, kmax))
    tvir_c[0, :] = np.sqrt(2.0) * 280.0
    topo_c = 20.0 * rng.normal(size=n)
    topo_c[0] = np.sqrt(2.0) * 200.0
    fields = dict(
        topo=topo_c,
        lnps=lnps_c,
        tvir=tvir_c,
        divg=1e-5 * rng.normal(size=(n, kmax)),
        vort=1e-5 * rng.normal(size=(n, kmax)),
        sphu=1e-3 * np.abs(rng.normal(size=(n, kmax))),
    )
    return write_ic(tmp_path / name, header, **fields)


def _time_series(pattern_t, pattern_q, pattern_ps, pattern_u, pattern_v,
                 lat, lon, coef):
    states = []
    for a in coef:
        states.append(mb09.State(
            t=a * pattern_t, q=a * pattern_q, ps=a * pattern_ps,
            u=a * pattern_u, v=a * pattern_v, lat=lat, lon=lon,
        ))
    return states


# ---------------------------------------------------------------------------
# 1. perturbacoes aleatorias
# ---------------------------------------------------------------------------
def test_random_deltas_amplitudes():
    nlat, nlon, nlev = 64, 128, 4
    state = _state(nlat, nlon, nlev)
    rng = np.random.default_rng(123)
    deltas = mb09.random_deltas(state, rng)

    assert set(deltas) == {"u", "v", "t", "ps", "q"}
    for key, target in (("u", 3.0), ("v", 3.0), ("t", 0.6), ("ps", 1.0)):
        values = deltas[key].ravel()
        rms = float(np.sqrt(np.mean(values ** 2)))
        assert abs(rms - target) <= 0.05 * target, (key, rms, target)
        assert abs(float(values.mean())) <= 0.05 * target

    qprof = np.asarray(mb.q_stddev_for_layers(nlev))
    for k in range(nlev):
        rms = float(np.sqrt(np.mean(deltas["q"][k] ** 2)))
        assert abs(rms - qprof[k]) <= max(0.05 * qprof[k], 1e-12)
        assert abs(float(deltas["q"][k].mean())) <= max(0.05 * qprof[k], 1e-12)

    # sigma por campo (inclusive vento como um unico sigma)
    d2 = mb09.random_deltas(state, np.random.default_rng(7),
                            sigma={"wind": 1.0, "T": 0.2, "P": 0.5, "Q": 0.0})
    assert abs(float(np.sqrt(np.mean(d2["u"] ** 2))) - 1.0) <= 0.05
    assert abs(float(np.sqrt(np.mean(d2["v"] ** 2))) - 1.0) <= 0.05
    assert abs(float(np.sqrt(np.mean(d2["t"] ** 2))) - 0.2) <= 0.02
    assert abs(float(np.sqrt(np.mean(d2["ps"] ** 2))) - 0.5) <= 0.05
    assert np.allclose(d2["q"], 0.0)


# ---------------------------------------------------------------------------
# 2. EOF dominante recuperado e sem matriz N x N
# ---------------------------------------------------------------------------
def test_build_eofs_recovers_dominant_mode_large_n():
    # N = 1 nivel * 100 x 200 = 20000 pontos por campo; um N x N em float64
    # teria 3.2 GB. O truque dual (11 x 11) roda com X de 11 x 20000.
    nlat, nlon, nlev = 100, 200, 1
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(11)
    pattern = rng.normal(size=(nlev, nlat, nlon))
    coef = np.linspace(0.1, 1.0, mb.M_TIMES)

    states = _time_series(pattern, np.zeros_like(pattern), np.zeros((nlat, nlon)),
                          np.zeros_like(pattern), np.zeros_like(pattern),
                          lat, lon, coef)
    diffs = {1: {"T": states}}
    eofs = mb09.build_eofs(diffs, preset="oper")

    mask = _mask(lat, lon, "ETR")
    mode = eofs["members"][0]["T"]["ETR"].mode
    assert mode.shape == (nlev, nlat, nlon)
    assert np.all(mode[~mask[None, :, :]] == 0.0)
    corr = np.corrcoef(mode[:, mask].ravel(), pattern[:, mask].ravel())[0, 1]
    assert corr > 0.99
    assert eofs["m_times"] == mb.M_TIMES
    assert eofs["n"] == 1


# ---------------------------------------------------------------------------
# 3. EOF ortonormal e autovalores decrescentes
# ---------------------------------------------------------------------------
def test_eofs_orthonormal_and_eigenvalues_sorted():
    nlat, nlon, nlev = 12, 24, 2
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(5)
    coef = np.linspace(0.2, 1.0, mb.M_TIMES)
    diffs = {}
    for member in (1, 2, 3):
        states = _time_series(
            rng.normal(size=(nlev, nlat, nlon)),
            np.abs(rng.normal(size=(nlev, nlat, nlon))) * 1e-4,
            rng.normal(size=(nlat, nlon)),
            rng.normal(size=(nlev, nlat, nlon)),
            rng.normal(size=(nlev, nlat, nlon)),
            lat, lon, coef * rng.uniform(0.8, 1.2, mb.M_TIMES),
        )
        diffs[member] = states

    eofs = mb09.build_eofs(diffs, preset="oper", keep_basis=True)
    for member in eofs["members"]:
        for field, regions in member.items():
            for region, mode in regions.items():
                evals = np.asarray(mode.eigenvalues)
                assert evals.size == mb.M_TIMES
                assert np.all(evals >= -1e-12)
                assert np.all(np.diff(evals) <= 1e-12)
                # modo normalizado sobre a regiao
                assert abs(float(np.linalg.norm(mode.mode.ravel())) - 1.0) < 1e-10
                # autovetores de G (11 x 11) ortonormais
                assert mode.eigenvectors.shape == (mb.M_TIMES, mb.M_TIMES)
                gram = mode.eigenvectors.conj().T @ mode.eigenvectors
                np.testing.assert_allclose(gram, np.eye(mb.M_TIMES), atol=1e-10)


# ---------------------------------------------------------------------------
# 4. aplicacao: 2n estados, simetria, q >= 0 e amplitudes alvo
# ---------------------------------------------------------------------------
def test_apply_modes_symmetry_and_target_amplitudes():
    nlat, nlon, nlev = 16, 32, 3
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(21)
    shape = (nlev, nlat, nlon)
    coef = np.linspace(0.1, 1.0, mb.M_TIMES)
    diffs = {}
    for member in range(1, mb.N_PERTURBATIONS + 1):
        diffs[member] = _time_series(
            rng.normal(size=shape), 1e-4 * np.abs(rng.normal(size=shape)),
            rng.normal(size=(nlat, nlon)),
            rng.normal(size=shape), rng.normal(size=shape),
            lat, lon, coef,
        )

    control = _state(nlat, nlon, nlev, q_value=1e-2)
    eofs = mb09.build_eofs(diffs, preset="trop")
    states = mb09.apply_modes(control, eofs, preset="trop")
    # sem preset explicito: usa o preset gravado no eofs
    fallback = mb09.apply_modes(control, eofs)
    np.testing.assert_allclose(fallback[0].t, states[0].t)

    assert len(states) == 2 * mb.N_PERTURBATIONS
    mask = _mask(lat, lon, "ETR")
    qprof = np.asarray(mb.q_stddev_for_layers(nlev))
    target = mb.AMPLITUDES_TARGET

    for i in range(0, len(states), 2):
        plus, minus = states[i], states[i + 1]
        dt = plus.t - control.t
        # ordem +p_k, -p_k e simetria em torno do controle
        np.testing.assert_allclose(dt, -(minus.t - control.t), atol=1e-10)
        assert abs(float(np.std(dt[:, mask])) - target["T"]) < 1e-6

        dps = plus.ps - control.ps
        np.testing.assert_allclose(dps, -(minus.ps - control.ps), atol=1e-10)
        assert abs(float(np.std(dps[mask])) - target["P"]) < 1e-6

        dwind = (plus.u - control.u) + 1j * (plus.v - control.v)
        assert abs(float(np.std(dwind[:, mask])) - target["wind"]) < 1e-6

        dq = plus.q - control.q
        for k in range(nlev):
            assert abs(float(np.std(dq[k][mask])) - qprof[k]) <= max(1e-12, 1e-6 * qprof[k])
        assert float(plus.q.min()) >= 0.0
        assert float(minus.q.min()) >= 0.0


# ---------------------------------------------------------------------------
# 5. escrita do ensemble (15 arquivos, tamanho do template)
# ---------------------------------------------------------------------------
def test_write_ensemble_15_files(tmp_path):
    control_path = _make_ic(tmp_path, mend=6, kmax=3, style="hybrid")
    control = mb09.State.from_ic_file(control_path)
    rng = np.random.default_rng(2)
    states = []
    for _ in range(2 * mb.N_PERTURBATIONS):
        delta = {
            "t": rng.normal(0.0, 0.1, control.t.shape),
            "q": rng.normal(0.0, 1e-6, control.q.shape),
            "u": rng.normal(0.0, 0.1, control.u.shape),
            "v": rng.normal(0.0, 0.1, control.v.shape),
            "ps": rng.normal(0.0, 0.1, control.ps.shape),
        }
        states.append(mb09.perturbed_state(control, delta))

    paths = mb09.write_ensemble(control_path, states, tmp_path / "ens", "tst",
                                q_profile="operacional_l42", remove_mean=False)
    assert len(paths) == mb.N_MEMBERS == 15
    assert [p.name for p in paths][:3] == ["tst_000.unf", "tst_001.unf", "tst_002.unf"]
    expected = control_path.stat().st_size
    for path in paths:
        assert path.exists()
        assert path.stat().st_size == expected
    assert len(list((tmp_path / "ens").glob("*.unf"))) == 15
    manifest = json.loads((tmp_path / "ens" / "members.json").read_text())
    assert len(manifest["members"]) == 14
    assert manifest["members"][0]["sign"] == "+"
    assert manifest["members"][1]["sign"] == "-"
    assert manifest["q_profile"] == "operacional_l42"
    assert manifest["remove_mean"] is False


# ---------------------------------------------------------------------------
# 6. roundtrip State -> IC -> State
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("style", ["hybrid", "sigma"])
def test_state_ic_roundtrip(tmp_path, style):
    src = _make_ic(tmp_path, mend=6, kmax=3, style=style)
    state = mb09.State.from_ic_file(src)
    out = tmp_path / ("rt_%s.unf" % style)
    state.to_ic_file(out, template=src)

    assert out.stat().st_size == src.stat().st_size
    back = mb09.State.from_ic_file(out)

    np.testing.assert_allclose(back.t, state.t, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(back.u, state.u, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(back.v, state.v, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(back.ps, state.ps, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(back.q, state.q, rtol=1e-6, atol=1e-12)
    # a topografia tambem volta (mesma transformada)
    np.testing.assert_allclose(back.topo, state.topo, rtol=1e-6, atol=1e-6)


# ---------------------------------------------------------------------------
# 7. fora da regiao a perturbacao e' nula
# ---------------------------------------------------------------------------
def test_region_perturbation_is_local():
    nlat, nlon, nlev = 20, 40, 2
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(3)
    coef = np.linspace(0.2, 1.0, mb.M_TIMES)
    diffs = {}
    for member in range(1, mb.N_PERTURBATIONS + 1):
        diffs[member] = _time_series(
            rng.normal(size=(nlev, nlat, nlon)), np.zeros((nlev, nlat, nlon)),
            np.zeros((nlat, nlon)),
            rng.normal(size=(nlev, nlat, nlon)),
            rng.normal(size=(nlev, nlat, nlon)),
            lat, lon, coef,
        )
    control = _state(nlat, nlon, nlev)
    eofs = mb09.build_eofs(diffs, preset="oper")
    states = mb09.apply_modes(control, eofs, preset="oper")

    mask = _mask(lat, lon, "ETR")
    outside = ~mask
    assert outside.any() and mask.any()
    for state in states:
        for field in ("t", "u", "v"):
            delta = getattr(state, field) - getattr(control, field)
            assert np.max(np.abs(delta[:, outside])) < 1e-12
            assert np.max(np.abs(delta[:, mask])) > 0.0


# ---------------------------------------------------------------------------
# extras: .npz de EOFs e build a partir de .ctl/.unf do PostGrib
# ---------------------------------------------------------------------------
def test_save_load_eofs_roundtrip(tmp_path):
    nlat, nlon, nlev = 10, 20, 2
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(9)
    coef = np.linspace(0.2, 1.0, mb.M_TIMES)
    diffs = {
        1: _time_series(rng.normal(size=(nlev, nlat, nlon)),
                        np.abs(rng.normal(size=(nlev, nlat, nlon))) * 1e-4,
                        rng.normal(size=(nlat, nlon)),
                        rng.normal(size=(nlev, nlat, nlon)),
                        rng.normal(size=(nlev, nlat, nlon)),
                        lat, lon, coef)
    }
    eofs = mb09.build_eofs(diffs, preset="trop", remove_mean=False,
                           q_profile="artigo_l28")
    path = mb09.save_eofs(tmp_path / "eof.npz", eofs)
    back = mb09.load_eofs(path)
    assert back["preset"] == eofs["preset"]
    assert back["n"] == eofs["n"] == 1
    # as escolhas de build sobrevivem ao .npz
    assert back["meta"]["remove_mean"] is False
    assert back["meta"]["q_profile"] == "artigo_l28"
    np.testing.assert_allclose(back["lat"], lat)
    np.testing.assert_allclose(back["lon"], lon)
    for field in eofs["eof_fields"]:
        for region in eofs["regions"]:
            a = eofs["members"][0][field][region]
            b = back["members"][0][field][region]
            np.testing.assert_allclose(b.mode, a.mode)
            np.testing.assert_allclose(b.eigenvalues, a.eigenvalues)
            np.testing.assert_allclose(b.sigma_i, a.sigma_i)

    # apply com o eofs recarregado reproduz o mesmo ensemble
    control = _state(nlat, nlon, nlev)
    s1 = mb09.apply_modes(control, eofs, preset="trop")
    s2 = mb09.apply_modes(control, back, preset="trop")
    for a, b in zip(s1, s2):
        np.testing.assert_allclose(a.t, b.t)
        np.testing.assert_allclose(a.u, b.u)


def _write_postgrib_case(tmp_path: Path, seed=0, yrev=False, nlat=8, nlon=16,
                         nlev=3, hours=None):
    """Escreve um .ctl/.unf no layout do PostGrib (M tempos no mesmo binario).

    'hours' e' a lista de horas de previsao (default FORECAST_TIMES_HOURS =
    6..36 passo 3). Os campos de cada hora sao gerados de forma DETERMINISTICA
    a partir de (seed, hora), de modo que subconjuntos do mesmo caso tenham
    valores identicos nas horas em comum (usado no teste da janela).
    ydef = latitudes sul -> norte e ps 2-D (1 registro). Com yrev, o disco fica
    norte -> sul; sem yrev, fica na propria ordem do ydef.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    hours = list(mb.FORECAST_TIMES_HOURS if hours is None else hours)
    ntimes = len(hours)
    lat = np.linspace(-60.0, 60.0, nlat)          # ydef sul -> norte
    lon = np.arange(nlon, dtype=float) * 360.0 / nlon
    fields = {name: [] for name in ("pslc", "tvir", "umes", "uvel", "vvel")}
    for hour in hours:
        rr = np.random.default_rng(int(seed) * 1000 + int(round(hour)))
        fields["pslc"].append(1000.0 + 5.0 * rr.normal(size=(1, nlat, nlon)))
        fields["tvir"].append(280.0 + 3.0 * rr.normal(size=(nlev, nlat, nlon)))
        fields["umes"].append(1e-2 + 1e-4 * rr.random((nlev, nlat, nlon)))
        fields["uvel"].append(8.0 * rr.normal(size=(nlev, nlat, nlon)))
        fields["vvel"].append(8.0 * rr.normal(size=(nlev, nlat, nlon)))
    order = ("pslc", "tvir", "umes", "uvel", "vvel")
    rows = slice(None, None, -1) if yrev else slice(None)
    # layout GrADS: t -> var -> nivel -> y -> x ; ps 2-D ocupa 1 registro
    blocks = []
    for it in range(ntimes):
        for name in order:
            blocks.append(fields[name][it][..., rows, :].ravel())
    raw = np.concatenate(blocks)
    (tmp_path / "gpost.unf").write_bytes(np.asarray(raw, dtype=">f4").tobytes())
    options = "yrev big_endian" if yrev else "big_endian"
    if len(hours) > 1:
        steps = np.diff(hours)
        if not np.allclose(steps, steps[0]):
            raise ValueError("hours precisa ser uniforme para o tdef linear")
        step_txt = "%ghr" % float(steps[0])
    else:
        step_txt = "3hr"
    first = datetime(2026, 1, 1) + timedelta(hours=float(hours[0]))
    tdef_first = first.strftime("%HZ%d%b%Y").upper()
    ctl = """dset ^gpost.unf
undef 1e20
options %s
xdef %d linear 0 %g
ydef %d levels %s
zdef %d levels %s
tdef %d linear %s %s
vars 5
pslc 0 99  hPa Surface Pressure
tvir %d 1 2 3 K Virtual Temperature
umes %d 1 2 3 kg/kg Specific Humidity
uvel %d 1 2 3 m/s Zonal Wind
vvel %d 1 2 3 m/s Meridional Wind
endvars
""" % (options, nlon, 360.0 / nlon, nlat, " ".join("%.3f" % v for v in lat),
       nlev, " ".join(str(k) for k in range(1, nlev + 1)), ntimes,
       tdef_first, step_txt, nlev, nlev, nlev, nlev)
    (tmp_path / "gpost.ctl").write_text(ctl, encoding="utf-8")
    return tmp_path / "gpost.ctl"


@pytest.mark.parametrize("yrev", [False, True])
def test_postgrib_reader_layout_and_yrev(tmp_path, yrev):
    ctl = _write_postgrib_case(tmp_path, seed=4, yrev=yrev)
    series = mb09._GradsSeries(ctl)
    state = series.state(1)
    # valores esperados: o layout t->var->nivel->y->x e' validado pelas medias
    assert abs(float(state.ps.mean()) - 1000.0) < 1.0
    assert abs(float(state.t.mean()) - 280.0) < 1.0
    assert abs(float(state.q.mean()) - 1e-2) < 1e-4
    assert abs(float(state.u.mean())) < 1.0
    assert state.shape == (3, 8, 16)
    # lat na ordem do ydef (sul -> norte), independente do yrev
    # ydef e' escrito com 3 casas (F10.3), como no .ctl do PostGrib
    np.testing.assert_allclose(series.lat, np.linspace(-60.0, 60.0, 8), atol=1e-3)
    # horas de previsao derivadas do tdef (analise implicita em 00Z do 1o dia)
    np.testing.assert_allclose(series.hours, mb.FORECAST_TIMES_HOURS, atol=1e-9)

    # yrev e sem yrev representam o MESMO campo fisico
    ctl2 = _write_postgrib_case(tmp_path / ("sem_yrev" if yrev else "com_yrev"),
                                seed=4, yrev=not yrev)
    state2 = mb09._GradsSeries(ctl2).state(1)
    np.testing.assert_allclose(state2.t, state.t, atol=1e-5)
    np.testing.assert_allclose(state2.u, state.u, atol=1e-5)


def test_build_from_postgrib_ctl(tmp_path):
    ctl0 = _write_postgrib_case(tmp_path / "c", seed=1)
    members = [_write_postgrib_case(tmp_path / ("m%d" % k), seed=10 + k)
               for k in range(1, 4)]
    eofs = mb09.build_eofs_from_files(ctl0, members, preset="trop")
    # --vars parcial: sobrescreve so o campo indicado e auto-detecta o resto
    assert mb09._parse_vars("ps=pslc,T=tvir") == {"P": "pslc", "T": "tvir"}
    eofs_override = mb09.build_eofs_from_files(
        ctl0, members, preset="trop", variables={"ps": "pslc"})
    assert eofs_override["n"] == 3
    assert eofs["n"] == 3
    assert eofs["m_times"] == mb.M_TIMES
    assert eofs["kmax"] == 3
    assert set(eofs["eof_fields"]) == {"P", "T", "Q", "wind"}
    mode = eofs["members"][0]["T"]["ETR"]
    assert mode.mode.shape == (3, 8, 16)
    assert np.all(np.isfinite(mode.mode))
    assert np.linalg.norm(mode.mode) == pytest.approx(1.0)
    assert np.all(mode.eigenvalues[:-1] >= mode.eigenvalues[1:])


# ---------------------------------------------------------------------------
# extras: janela 6..36 h (select_window e build)
# ---------------------------------------------------------------------------
def test_select_window_indices():
    times = np.arange(0.0, 40.0, 3.0)          # 0, 3, ..., 39 h
    idx = mb09.select_window(times)
    assert idx.tolist() == list(range(2, 13))  # 11 indices: 6..36 h
    assert times[idx].tolist() == list(range(6, 37, 3))
    assert mb09.select_window(times, 0, 39, 3).size == 14
    # lista explicita e' devolvida em ordem cronologica
    np.testing.assert_array_equal(mb09.select_times(times, [9, 6]), [2, 3])
    # janela com horarios ausentes: erro claro com os tempos disponiveis
    with pytest.raises(ValueError, match="disponiveis"):
        mb09.select_window([0.0, 3.0, 6.0])
    with pytest.raises(ValueError):
        mb09.select_window(times, 6, 36, 0)
    # parsing do CLI
    assert mb09._parse_window("6:36:3") == (6.0, 36.0, 3.0)
    assert mb09._parse_window("0:39:3") == (0.0, 39.0, 3.0)
    with pytest.raises(ValueError):
        mb09._parse_window("6:36")


def _window_case(root, hours, seeds):
    control = _write_postgrib_case(root / "c", seed=seeds[0], hours=hours)
    members = [_write_postgrib_case(root / ("m%d" % k), seed=seeds[k], hours=hours)
               for k in range(1, len(seeds))]
    return control, members


def test_build_window_ignores_extra_times(tmp_path):
    full = list(range(0, 40, 3))               # 0..39 h (3 tempos extras)
    sel = list(range(6, 37, 3))                # janela do MB09
    seeds = (31, 101, 102, 103)

    c_full, m_full = _window_case(tmp_path / "full", full, seeds)
    c_sel, m_sel = _window_case(tmp_path / "sel", sel, seeds)
    eof_full = mb09.build_eofs_from_files(c_full, m_full, preset="trop")
    eof_sel = mb09.build_eofs_from_files(c_sel, m_sel, preset="trop")

    assert eof_full["m_times"] == mb.M_TIMES
    assert eof_full["meta"]["n_ignored_times"] == 3
    assert eof_full["meta"]["hours"] == [float(h) for h in sel]
    assert eof_sel["meta"]["n_ignored_times"] == 0
    # tempos extras nao mudam os EOFs
    for field in eof_full["eof_fields"]:
        for region in eof_full["regions"]:
            a = eof_full["members"][0][field][region]
            b = eof_sel["members"][0][field][region]
            np.testing.assert_allclose(a.mode, b.mode, atol=1e-12)
            np.testing.assert_allclose(a.eigenvalues, b.eigenvalues, atol=1e-12)
            np.testing.assert_allclose(a.sigma_i, b.sigma_i, atol=1e-12)

    # --times explicito e analysis-time explicito reproduzem a mesma janela
    eof_times = mb09.build_eofs_from_files(c_full, m_full, preset="trop", hours=sel)
    eof_at = mb09.build_eofs_from_files(
        c_full, m_full, preset="trop", analysis_time="2026010100")
    np.testing.assert_allclose(eof_times["members"][0]["T"]["ETR"].mode,
                               eof_sel["members"][0]["T"]["ETR"].mode, atol=1e-12)
    np.testing.assert_allclose(eof_at["members"][0]["T"]["ETR"].mode,
                               eof_sel["members"][0]["T"]["ETR"].mode, atol=1e-12)

    # .ctl sem os horarios da janela: erro lista os tempos disponiveis
    c_short = _write_postgrib_case(tmp_path / "short", seed=31, hours=[0, 3, 6])
    with pytest.raises(ValueError, match="tempos disponiveis"):
        mb09.build_eofs_from_files(c_short, [c_short], preset="trop")


def test_build_cli_warns_about_extra_times(tmp_path, capsys):
    full = list(range(0, 40, 3))
    c_ctl, m_ctls = _window_case(tmp_path, full, (51, 201, 202, 203))
    out = tmp_path / "eof.npz"
    rc = mb09.main(["build", "--control", str(c_ctl),
                    "--members", *[str(p) for p in m_ctls],
                    "--preset", "trop", "--n", "3", "--out", str(out)])
    text = capsys.readouterr().out
    assert rc == 0
    assert out.exists()
    assert "aviso:" in text and "3 tempo" in text
    assert "horas=6,9,12" in text


# ---------------------------------------------------------------------------
# extras: n (numero de perturbacoes) configuravel, impar, default 7
# ---------------------------------------------------------------------------
def test_random_deltas_n_list():
    state = _state(nlat=8, nlon=16, nlev=2)
    rng = np.random.default_rng(0)
    single = mb09.random_deltas(state, rng)
    assert isinstance(single, dict) and set(single) == {"u", "v", "t", "ps", "q"}
    one = mb09.random_deltas(state, rng, n=1)
    assert isinstance(one, list) and len(one) == 1
    three = mb09.random_deltas(state, rng, n=3)
    assert len(three) == 3
    assert not np.allclose(three[0]["u"], three[1]["u"])
    assert all(set(d) == {"u", "v", "t", "ps", "q"} for d in three)


def test_n_must_be_odd_and_positive():
    state = _state()
    diffs = {1: [_state() for _ in range(mb.M_TIMES)]}
    for bad in (0, -1, 2, 4, 6):
        with pytest.raises(ValueError, match="n deve ser"):
            mb09.random_deltas(state, np.random.default_rng(0), n=bad)
    for bad in (0, 2, 4):
        with pytest.raises(ValueError, match="n deve ser"):
            mb09.build_eofs(diffs, preset="trop", n=bad)
    eofs = mb09.build_eofs(diffs, preset="trop", n=1)
    for bad in (0, 2, 4):
        with pytest.raises(ValueError, match="n deve ser"):
            mb09.apply_modes(_state(), eofs, preset="trop", n=bad)
    with pytest.raises(ValueError, match="n deve ser"):
        mb09.run_mb09("nao-existe.unf", [], "out", n=2)


@pytest.mark.parametrize("n", [1, 3, 5, 9])
def test_ensemble_size_for_n(tmp_path, n):
    control_path = _make_ic(tmp_path, mend=6, kmax=3, style="hybrid")
    control = mb09.State.from_ic_file(control_path)
    rng = np.random.default_rng(100 + n)
    shape = (control.nlev, control.nlat, control.nlon)
    coef = np.linspace(0.1, 1.0, mb.M_TIMES)
    diffs = {}
    for k in range(1, n + 1):
        diffs[k] = _time_series(
            rng.normal(size=shape), 1e-4 * np.abs(rng.normal(size=shape)),
            rng.normal(size=(control.nlat, control.nlon)),
            rng.normal(size=shape), rng.normal(size=shape),
            control.lat, control.lon, coef,
        )
    eofs = mb09.build_eofs(diffs, preset="trop", n=n)
    assert eofs["n"] == n
    states = mb09.apply_modes(control, eofs, preset="trop", n=n)
    assert len(states) == 2 * n
    for i in range(0, len(states), 2):
        plus, minus = states[i], states[i + 1]
        np.testing.assert_allclose(plus.t - control.t, -(minus.t - control.t), atol=1e-10)
        np.testing.assert_allclose(plus.u - control.u, -(minus.u - control.u), atol=1e-10)
        assert float(plus.q.min()) >= 0.0
        assert float(minus.q.min()) >= 0.0

    prefix = "ens%02d" % n
    out_dir = tmp_path / ("ens%d" % n)
    paths = mb09.write_ensemble(control_path, states, out_dir, prefix)
    assert len(paths) == 2 * n + 1
    assert paths[0].name == prefix + "_000.unf"
    assert paths[-1].name == prefix + "_%03d.unf" % (2 * n)
    expected = control_path.stat().st_size
    assert all(p.stat().st_size == expected for p in paths)
    manifest = json.loads((out_dir / "members.json").read_text())
    assert manifest["n"] == n
    assert manifest["n_members"] == 2 * n + 1
    assert len(manifest["members"]) == 2 * n


def test_n_truncates_members_and_apply():
    nlat, nlon, nlev = 12, 24, 2
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(33)
    coef = np.linspace(0.1, 1.0, mb.M_TIMES)
    shape = (nlev, nlat, nlon)
    diffs = {}
    for k in range(1, 6):                     # 5 membros disponiveis
        diffs[k] = _time_series(
            rng.normal(size=shape), np.abs(rng.normal(size=shape)) * 1e-4,
            rng.normal(size=(nlat, nlon)),
            rng.normal(size=shape), rng.normal(size=shape),
            lat, lon, coef,
        )
    full = mb09.build_eofs(diffs, preset="trop", n=3)
    ref = mb09.build_eofs({k: diffs[k] for k in (1, 2, 3)}, preset="trop", n=3)
    assert full["n"] == 3
    assert full["meta"]["n_available"] == 5 and full["meta"]["n_used"] == 3
    for field in full["eof_fields"]:
        for region in full["regions"]:
            np.testing.assert_allclose(
                full["members"][0][field][region].mode,
                ref["members"][0][field][region].mode, atol=1e-12)
    with pytest.raises(ValueError, match="membro"):
        mb09.build_eofs(diffs, preset="trop", n=7)

    control = _state(nlat, nlon, nlev)
    e1 = mb09.apply_modes(control, full, preset="trop", n=1)
    e3 = mb09.apply_modes(control, full, preset="trop", n=3)
    assert len(e1) == 2 and len(e3) == 6
    np.testing.assert_allclose(e1[0].t, e3[0].t)
    np.testing.assert_allclose(e1[1].t, e3[1].t)
    with pytest.raises(ValueError, match="pede mais"):
        mb09.apply_modes(control, full, preset="trop", n=5)


def test_build_from_files_n_truncates_and_warns(tmp_path, capsys):
    ctl0 = _write_postgrib_case(tmp_path / "c", seed=1)
    members = [_write_postgrib_case(tmp_path / ("m%d" % k), seed=10 + k)
               for k in range(1, 6)]
    out = tmp_path / "eof.npz"
    rc = mb09.main(["build", "--control", str(ctl0),
                    "--members", *[str(p) for p in members],
                    "--preset", "trop", "--n", "3", "--out", str(out)])
    text = capsys.readouterr().out
    assert rc == 0
    assert "aviso:" in text and "usando os 3 primeiros" in text
    assert mb09.load_eofs(out)["n"] == 3

    # n maior que o numero de membros: erro (main devolve 1)
    rc = mb09.main(["build", "--control", str(ctl0),
                    "--members", *[str(p) for p in members[:3]],
                    "--preset", "trop", "--n", "5",
                    "--out", str(tmp_path / "e2.npz")])
    assert rc == 1


def test_cli_random_n_modes(tmp_path):
    control_path = _make_ic(tmp_path, mend=6, kmax=3, style="hybrid")
    rc = mb09.main(["random", "--ic", str(control_path),
                    "--out-dir", str(tmp_path / "r"), "--n", "3", "--seed", "1"])
    assert rc == 0
    files = sorted((tmp_path / "r").glob("mb09rand_*.unf"))
    assert len(files) == 3
    manifest = json.loads((tmp_path / "r" / "random.json").read_text())
    assert manifest["n"] == 3 and len(manifest["members"]) == 3

    rc = mb09.main(["random", "--ic", str(control_path),
                    "--out-dir", str(tmp_path / "r2"), "--n", "4"])
    assert rc == 1
    assert not list((tmp_path / "r2").glob("*.unf"))


# ---------------------------------------------------------------------------
# extras: remove_mean e q_profile (flags do operacional)
# ---------------------------------------------------------------------------
def test_remove_mean_changes_eof_and_meta():
    nlat, nlon, nlev = 12, 24, 2
    lat, lon = _grid(nlat, nlon)
    rng = np.random.default_rng(77)
    pat_var = rng.normal(size=(nlev, nlat, nlon))
    pat_mean = rng.normal(size=(nlev, nlat, nlon))
    coef = np.linspace(0.2, 1.0, mb.M_TIMES)
    states = []
    for a in coef:
        t = a * pat_var + 3.0 * pat_mean       # media temporal nao nula
        states.append(mb09.State(
            t=t, q=np.full_like(t, 1e-2), u=np.zeros_like(t), v=np.zeros_like(t),
            ps=np.full((nlat, nlon), 1000.0), lat=lat, lon=lon))
    diffs = {1: {"T": states}}

    yes = mb09.build_eofs(diffs, preset="oper", n=1, remove_mean=True)
    no = mb09.build_eofs(diffs, preset="oper", n=1, remove_mean=False)
    assert yes["meta"]["remove_mean"] is True
    assert no["meta"]["remove_mean"] is False
    assert yes["meta"]["q_profile"] == mb09.DEFAULT_Q_PROFILE

    mask = _mask(lat, lon, "ETR")
    my = yes["members"][0]["T"]["ETR"].mode
    mn = no["members"][0]["T"]["ETR"].mode
    corr = np.corrcoef(my[:, mask].ravel(), pat_var[:, mask].ravel())[0, 1]
    assert abs(corr) > 0.99                     # EOF classico = parte variavel
    assert not np.allclose(my, mn, atol=1e-3)   # sem remover a media difere
    assert abs(np.corrcoef(my[:, mask].ravel(), mn[:, mask].ravel())[0, 1]) < 0.9


def test_q_profile_validation_and_cli_meta(tmp_path, capsys):
    state = _state()
    with pytest.raises(ValueError, match="q_profile"):
        mb09.random_deltas(state, np.random.default_rng(0), q_profile="invalido")
    diffs = {1: [_state() for _ in range(mb.M_TIMES)]}
    with pytest.raises(ValueError, match="q_profile"):
        mb09.build_eofs(diffs, preset="trop", q_profile="invalido")
    eofs = mb09.build_eofs(diffs, preset="trop", n=1)
    with pytest.raises(ValueError, match="q_profile"):
        mb09.apply_modes(state, eofs, preset="trop", q_profile="invalido")
    # o perfil operacional exige 42 camadas na hora de calcular sigma_f
    with pytest.raises(ValueError, match="42"):
        mb09.apply_modes(state, eofs, preset="trop", q_profile="operacional_l42")

    # CLI build registra e persiste as duas escolhas no eof.npz
    ctl0 = _write_postgrib_case(tmp_path / "c", seed=1)
    members = [_write_postgrib_case(tmp_path / ("m%d" % k), seed=10 + k)
               for k in range(1, 4)]
    out = tmp_path / "eof.npz"
    rc = mb09.main(["build", "--control", str(ctl0),
                    "--members", *[str(p) for p in members],
                    "--preset", "trop", "--n", "3",
                    "--remove-mean", "no",
                    "--q-profile", "operacional_l42",
                    "--out", str(out)])
    text = capsys.readouterr().out
    assert rc == 0
    assert "remove_mean=no" in text
    assert "q_profile=operacional_l42" in text
    saved = mb09.load_eofs(out)
    assert saved["meta"]["remove_mean"] is False
    assert saved["meta"]["q_profile"] == "operacional_l42"


def test_operacional_l42_zeroes_upper_q():
    nlat, nlon, nlev = 4, 8, 42
    control = _state(nlat, nlon, nlev)
    rng = np.random.default_rng(8)
    shape = (nlev, nlat, nlon)
    coef = np.linspace(0.1, 1.0, mb.M_TIMES)
    diffs = {1: _time_series(
        rng.normal(size=shape), np.abs(rng.normal(size=shape)) * 1e-4,
        rng.normal(size=(nlat, nlon)),
        rng.normal(size=shape), rng.normal(size=shape),
        control.lat, control.lon, coef)}

    eofs = mb09.build_eofs(diffs, preset="trop", n=1,
                           q_profile="operacional_l42")
    assert eofs["meta"]["q_profile"] == "operacional_l42"
    states = mb09.apply_modes(control, eofs, preset="trop", n=1)
    dq = states[0].q - control.q
    # camadas 18..42 (1-based, indices 17..41) ficam com sigma_f = 0
    assert np.all(dq[17:] == 0.0)
    assert np.any(dq[:17] != 0.0)
    # o perfil vem de eofs["meta"] quando nao e' passado
    states2 = mb09.apply_modes(control, eofs, preset="trop", n=1)
    np.testing.assert_allclose(states2[0].q, states[0].q)
    # override explicito coincide; o random tambem usa o perfil operacional
    states3 = mb09.apply_modes(control, eofs, preset="trop", n=1,
                               q_profile="operacional_l42")
    np.testing.assert_allclose(states3[0].q, states[0].q)
    delta = mb09.random_deltas(control, np.random.default_rng(1),
                               q_profile="operacional_l42")
    assert np.all(delta["q"][17:] == 0.0)
    assert np.any(delta["q"][:17] != 0.0)

