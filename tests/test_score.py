"""Testes das metricas de verificacao do MB09 (bamci.score).

Cobrem as funcoes puras (pac, rmse, spread, BSS com REL/RES/UNC, PECA), a
selecao por regiao (bamci.mb09_tables.region_mask) e a CLI com .ctl/.unf
sinteticos no layout do PostGrib (t -> var -> nivel -> y -> x, yrev).
"""

import json
from pathlib import Path

import numpy as np
import pytest

import bamci.score as sc
from bamci import mb09_tables as mb


# ---------------------------------------------------------------------------
# 1. PAC
# ---------------------------------------------------------------------------
def test_pac_identico_e_oposto():
    a = np.array([1.0, 2.0, 3.0])
    assert sc.pac(a, a, 0.0) == pytest.approx(1.0)
    assert sc.pac(a, -a, 0.0) == pytest.approx(-1.0)


def test_pac_valor_conferido_a_mao():
    fcst = np.array([1.0, 2.0, 3.0])
    obs = np.array([1.0, 2.0, 4.0])
    # af = [1,2,3], ao = [1,2,4]; r = 17/sqrt(14*21)
    esperado = 17.0 / np.sqrt(14.0 * 21.0)
    assert sc.pac(fcst, obs, 0.0) == pytest.approx(esperado, rel=1e-12)


def test_pac_e_invariante_a_deslocamento_da_climatologia():
    fcst = np.array([1.0, 2.0, 3.0, 5.0])
    obs = np.array([2.0, 2.0, 4.0, 4.0])
    clim = np.full(4, 2.5)
    # deslocar fcst, obs e clim pela mesma constante nao muda as anomalias
    assert sc.pac(fcst + 10.0, obs + 10.0, clim + 10.0) == pytest.approx(
        sc.pac(fcst, obs, clim), rel=1e-12)


def test_pac_nan_quando_anomalia_e_nula():
    a = np.array([1.0, 2.0, 3.0])
    assert np.isnan(sc.pac(a, a, a))       # anomalia identicamente nula


def test_anomaly_simples():
    np.testing.assert_allclose(sc.anomaly([1.0, 3.0], [0.5, 0.5]), [0.5, 2.5])
    np.testing.assert_allclose(sc.anomaly([1.0, 3.0], 1.0), [0.0, 2.0])


# ---------------------------------------------------------------------------
# 2. RMSE e spread
# ---------------------------------------------------------------------------
def test_rmse_e_spread_exatos():
    members = np.array([[0.0, 0.0], [2.0, 2.0]])
    assert sc.rmse(np.array([1.0, 1.0]), np.array([0.0, 0.0])) == pytest.approx(1.0)
    assert sc.rmse(np.array([0.0, 0.0]), np.array([3.0, 4.0])) == pytest.approx(
        np.sqrt(12.5))
    media = sc.ensemble_mean(members)
    np.testing.assert_allclose(media, [1.0, 1.0])
    np.testing.assert_allclose(sc.spread(members), [1.0, 1.0])
    # ddof=1 no caso de 2 membros: sqrt(2)/... -> desvio amostral
    np.testing.assert_allclose(sc.spread(members, ddof=1), [np.sqrt(2.0)] * 2)


def test_spread_zero_quando_todos_os_membros_sao_iguais():
    members = np.ones((5, 3, 2))
    np.testing.assert_allclose(sc.spread(members), np.zeros((3, 2)))


def test_ensemble_mean_de_um_membro_e_o_proprio_campo():
    field = np.arange(6.0).reshape(2, 3)
    np.testing.assert_allclose(sc.ensemble_mean(field[None, ...]), field)


def test_erros_de_forma():
    with pytest.raises(ValueError):
        sc.spread(np.zeros((1, 3)))
    with pytest.raises(ValueError):
        sc.rmse(np.zeros(0), np.zeros(0))


# ---------------------------------------------------------------------------
# 3. BSS / REL / RES
# ---------------------------------------------------------------------------
def test_bss_previsao_perfeita():
    prob = np.array([0.0, 1.0, 0.0, 1.0, 1.0, 0.0])
    obs = prob.copy()
    out = sc.brier_skill_score(prob, obs)
    assert out["bss"] == pytest.approx(1.0)
    assert out["rel"] == pytest.approx(0.0, abs=1e-15)
    assert out["res"] == pytest.approx(out["unc"])
    assert out["n"] == 6


def test_bss_probabilidade_climatologica_constante():
    obs = np.array([0.0, 1.0, 0.0, 1.0, 1.0])
    prob = np.full(5, obs.mean())          # 0.6 em todas as caixas
    out = sc.brier_skill_score(prob, obs)
    assert out["bss"] == pytest.approx(0.0)
    assert out["rel"] == pytest.approx(0.0)
    assert out["res"] == pytest.approx(0.0)
    assert out["base_rate"] == pytest.approx(0.6)


def test_bss_identidade_res_rel_unc():
    """BS = REL - RES + UNC e BSS = (RES - REL)/UNC com uma probabilidade por caixa.

    A identidade de Wilks pressupoe probabilidade constante dentro da caixa; por
    isso o teste usa probabilidades de um ensemble de 3 membros (0, 1/3, 2/3, 1),
    cada uma sozinha na sua caixa de 0.1.
    """
    rng = np.random.default_rng(7)
    n = 400
    prob = rng.choice([0.0, 1 / 3, 2 / 3, 1.0], size=n)
    obs = (rng.random(n) < prob).astype(float)
    out = sc.brier_skill_score(prob, obs)
    fechamento = (out["res"] - out["rel"]) / out["unc"]
    assert fechamento == pytest.approx(out["bss"], abs=1e-12)
    brier = np.mean((prob - obs) ** 2)
    assert brier == pytest.approx(out["rel"] - out["res"] + out["unc"], abs=1e-12)


def test_bss_identidade_aproximada_com_probabilidade_continua():
    """Com p continuo a identidade so vale ate o termo intra-caixa (Wilks)."""
    rng = np.random.default_rng(8)
    prob = rng.random(2000)
    obs = (rng.random(2000) < prob).astype(float)
    out = sc.brier_skill_score(prob, obs)
    brier = np.mean((prob - obs) ** 2)
    assert abs(brier - (out["rel"] - out["res"] + out["unc"])) < 0.01
    assert out["rel"] >= 0.0 and out["unc"] > 0.0


def test_bss_igual_para_probabilidade_discreta_de_3_membros():
    # 3 membros -> probabilidades 0, 1/3, 2/3, 1 (caixas de 0.1 agrupam isso)
    prob = np.array([0.0, 1 / 3, 2 / 3, 1.0])
    obs = np.array([0.0, 0.0, 1.0, 1.0])
    out = sc.brier_skill_score(prob, obs)
    assert -1.0 <= out["bss"] <= 1.0
    assert out["n"] == 4


def test_bss_caixas_customizadas():
    prob = np.array([0.05, 0.05, 0.95, 0.95])
    obs = np.array([0.0, 0.0, 1.0, 1.0])
    out = sc.brier_skill_score(prob, obs, thresholds=[0.0, 0.5, 1.0])
    # REL = 0.5*(0.05)^2 + 0.5*(0.05)^2 = 0.0025; RES = UNC = 0.25
    assert out["rel"] == pytest.approx(0.0025)
    assert out["bss"] == pytest.approx(0.99)
    assert out["bins"] == [0.0, 0.5, 1.0]


def test_bss_unc_zero_retorna_nan():
    out = sc.brier_skill_score(np.array([0.2, 0.3]), np.array([0.0, 0.0]))
    assert np.isnan(out["bss"])
    assert out["unc"] == 0.0


def test_bss_erros_de_argumento():
    with pytest.raises(ValueError):
        sc.brier_skill_score([0.5], [0.0, 1.0])            # tamanhos diferentes
    with pytest.raises(ValueError):
        sc.brier_skill_score([1.5], [0.0])                 # prob fora de [0,1]
    with pytest.raises(ValueError):
        sc.brier_skill_score([0.5], [0.5])                 # obs nao binaria
    with pytest.raises(ValueError):
        sc.brier_skill_score([0.5], [0.0], thresholds=[0.5])  # 1 borda
    with pytest.raises(ValueError):
        sc.brier_skill_score([0.5], [0.0], thresholds=[0.5, 0.2])


def test_ensemble_probability():
    members = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    np.testing.assert_allclose(sc.ensemble_probability(members, 0.5), [2 / 3, 0.0])
    np.testing.assert_allclose(
        sc.ensemble_probability(members, 0.5, above=False), [1 / 3, 1.0])


# ---------------------------------------------------------------------------
# 4. PECA
# ---------------------------------------------------------------------------
def test_peca_alinhado_melhora_rmse():
    rng = np.random.default_rng(11)
    truth = rng.normal(size=(3, 4))
    truth = truth - truth.mean()
    erro = rng.normal(size=(3, 4))
    erro = erro - erro.mean()
    mean = truth + erro                     # erro do membro medio = "erro"
    members = np.stack([mean + erro, mean - erro])   # perturbacoes +-e
    out = sc.peca(members, truth)
    assert out["sin"] == pytest.approx(1.0, abs=1e-12)
    # media com sinal de um par +p/-p e zero (por simetria)
    assert out["sin_signed"] == pytest.approx(0.0, abs=1e-12)
    opt = out["opt"]
    assert opt["rmse_after"] < opt["rmse_before"]
    assert opt["gain"] > 0.0
    assert opt["peca"] == pytest.approx(1.0, abs=1e-12)
    assert opt["rmse_after"] == pytest.approx(0.0, abs=1e-9)


def test_peca_ortogonal_nao_ganha():
    rng = np.random.default_rng(12)
    truth = rng.normal(size=(3, 4))
    erro = rng.normal(size=(3, 4))
    mean = truth + erro
    v = rng.normal(size=(3, 4))
    v = v - v.mean()
    e0 = erro - erro.mean()
    # projeta v no complemento ortogonal do erro
    p = v - (v * e0).sum() / (e0 * e0).sum() * e0
    p = p - p.mean()
    out = sc.peca(np.stack([mean + p, mean - p]), truth)
    assert abs(out["sin"]) < 1e-12
    assert out["opt"]["gain"] == pytest.approx(0.0, abs=1e-12)
    assert abs(out["opt"]["peca"]) < 1e-12


def test_peca_signed_preserva_o_sinal():
    rng = np.random.default_rng(13)
    truth = rng.normal(size=(3, 4))
    erro = rng.normal(size=(3, 4))
    mean = truth + erro
    v = erro - erro.mean()
    # coeficientes somando zero => perturbacoes +2v, -v, -v (media do ensemble
    # continua sendo 'mean'); correlacoes com o erro: +1, -1, -1
    members = np.stack([mean + 2.0 * v, mean - v, mean - v])
    out_abs = sc.peca(members, truth)
    out_sig = sc.peca(members, truth, signed=True)
    np.testing.assert_allclose(out_abs["sin_per_member"], [1.0, -1.0, -1.0], atol=1e-12)
    assert out_sig["sin"] == pytest.approx(out_sig["sin_signed"])
    assert out_sig["sin_signed"] == pytest.approx(-1.0 / 3.0)
    assert out_abs["sin"] == pytest.approx(1.0)


def test_peca_modos_e_erros():
    rng = np.random.default_rng(14)
    truth = rng.normal(size=(2, 3))
    members = np.stack([truth, truth])
    assert sc.peca(members, truth, mode="sin")["opt"] is None
    assert sc.peca(members, truth, mode="opt")["sin"] is None
    with pytest.raises(ValueError):
        sc.peca(members, truth, mode="ambos")
    with pytest.raises(ValueError):
        sc.peca(members, truth[:1])           # forma incompativel
    with pytest.raises(ValueError):
        sc.peca(members[:1], truth)           # so 1 membro


# ---------------------------------------------------------------------------
# 5. Regiao
# ---------------------------------------------------------------------------
def test_select_region_escolhe_o_subconjunto_certo():
    lat = np.array([-45.0, 0.0, 45.0])
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    field = np.arange(12.0).reshape(3, 4)

    np.testing.assert_allclose(sc.select_region(field, lat, lon, "global"), field.ravel())
    np.testing.assert_allclose(sc.select_region(field, lat, lon, "nh"), [8.0, 9.0, 10.0, 11.0])
    np.testing.assert_allclose(sc.select_region(field, lat, lon, "tr"), [4.0, 5.0, 6.0, 7.0])
    np.testing.assert_allclose(sc.select_region(field, lat, lon, "sh"), [0.0, 1.0, 2.0, 3.0])
    # NSA: lat 0 e lon 270 (20S-20N, 100W-10W); SSA: lat -45 e lon 270
    np.testing.assert_allclose(sc.select_region(field, lat, lon, "nsa"), [7.0])
    np.testing.assert_allclose(sc.select_region(field, lat, lon, "ssa"), [3.0])
    # aliases: maiusculo, None e tupla
    np.testing.assert_allclose(sc.select_region(field, lat, lon, "NH"), [8.0, 9.0, 10.0, 11.0])
    np.testing.assert_allclose(sc.select_region(field, lat, lon, None), field.ravel())
    np.testing.assert_allclose(
        sc.select_region(field, lat, lon, (-10.0, 10.0, 0.0, 360.0)), [4.0, 5.0, 6.0, 7.0])


def test_select_region_em_campo_3d_e_broadcast():
    lat = np.array([-45.0, 45.0])
    lon = np.array([0.0, 180.0])
    field = np.arange(8.0).reshape(2, 2, 2)      # (nivel, lat, lon)
    out = sc.select_region(field, lat, lon, "nh")
    assert out.shape == (2, 2)
    np.testing.assert_allclose(out, [[2.0, 3.0], [6.0, 7.0]])


def test_select_region_erros():
    lat = np.array([-45.0, 45.0])
    lon = np.array([0.0, 180.0])
    with pytest.raises(KeyError):
        sc.select_region(np.zeros((2, 2)), lat, lon, "xx")
    with pytest.raises(ValueError):
        sc.select_region(np.zeros((2, 2)), lat, lon, (10.0, 20.0, 0.0, 10.0))
    with pytest.raises(ValueError):
        sc.select_region(np.zeros(4), lat, lon, "nh")


def test_regiao_bate_com_a_mascara_do_mb09_tables():
    lat = np.array([-80.0, -40.0, -20.0, 0.0, 20.0, 40.0, 80.0])
    lon = np.linspace(0.0, 355.0, 72)
    field = np.ones((lat.size, lon.size))
    for nome in ("nh", "sh", "tr", "nsa", "ssa"):
        mask = mb.region_mask(lat[:, None], lon[None, :], nome.upper())
        assert sc.select_region(field, lat, lon, nome).size == int(mask.sum())


# ---------------------------------------------------------------------------
# 6. CLI com .ctl/.unf sintetico
# ---------------------------------------------------------------------------
def _write_case(tmp_path: Path, arrays, lat, lon, ref: Path | None = None) -> Path:
    """Escreve um .ctl/.unf no layout do PostGrib (t -> var -> nivel -> y -> x).

    arrays: dict nome -> (nt, nlev, nlat, nlon) para 3-D ou (nt, nlat, nlon)
    para 2-D. 'pslc' e' o unico 2-D; os demais sao 3-D com nlev niveis.
    ydef e' escrito sul -> norte e o disco norte -> sul (options yrev).
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    order = [n for n in ("pslc", "tvir", "umes", "uvel", "vvel") if n in arrays]
    first = np.asarray(arrays[order[0]], dtype=np.float64)
    nt = first.shape[0]
    nlev = 1
    for name in order:                      # o nlev vem de uma variavel 3-D
        arr = np.asarray(arrays[name])
        if arr.ndim == 4:
            nlev = arr.shape[1]
            break
    blocks = []
    for it in range(nt):
        for name in order:
            block = np.asarray(arrays[name][it], dtype=np.float64)
            if block.ndim == 2:
                block = block[None, ...]
            blocks.append(block[..., ::-1, :].ravel())      # yrev
    raw = np.concatenate(blocks)
    (tmp_path / "gpost.unf").write_bytes(np.asarray(raw, dtype=">f4").tobytes())

    zlev = " ".join(str(k) for k in range(1, nlev + 1))
    lines = [
        "dset ^gpost.unf",
        "undef 1e20",
        "options yrev big_endian",
        "xdef %d linear 0 %g" % (lon.size, 360.0 / lon.size),
        "ydef %d levels %s" % (lat.size, " ".join("%.3f" % v for v in lat)),
        "zdef %d levels %s" % (nlev, zlev),
        "tdef %d linear 06Z01JAN2026 3hr" % nt,
        "vars %d" % len(order),
    ]
    for name in order:
        if name == "pslc":
            lines.append("pslc 0 99 hPa Surface Pressure")
        else:
            lines.append("%s %d %s K/units %s field" % (name, nlev, zlev, name))
    lines.append("endvars")
    ctl = tmp_path / "gpost.ctl"
    ctl.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return ctl


def _case_arrays(nt=2, nlev=2, nlat=3, nlon=4, seed=0):
    rng = np.random.default_rng(seed)
    truth = {
        "pslc": 1000.0 + rng.normal(0.0, 1.0, (nt, nlat, nlon)),
        "tvir": 280.0 + rng.normal(0.0, 2.0, (nt, nlev, nlat, nlon)),
        "umes": 1e-2 * np.ones((nt, nlev, nlat, nlon)),
        "uvel": rng.normal(0.0, 5.0, (nt, nlev, nlat, nlon)),
        "vvel": rng.normal(0.0, 5.0, (nt, nlev, nlat, nlon)),
    }
    delta = rng.normal(0.0, 1.0, (nt, nlev, nlat, nlon))
    members = [truth, truth]
    members[0] = {k: v.copy() for k, v in truth.items()}
    members[1] = {k: v.copy() for k, v in truth.items()}
    members[0]["tvir"] = truth["tvir"] + delta
    members[1]["tvir"] = truth["tvir"] - delta
    return truth, members


def test_cli_gera_json_valido(tmp_path, capsys):
    lat = np.array([-45.0, 0.0, 45.0])
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    truth, members = _case_arrays()
    truth_ctl = _write_case(tmp_path / "truth", truth, lat, lon)
    m1 = _write_case(tmp_path / "m1", members[0], lat, lon)
    m2 = _write_case(tmp_path / "m2", members[1], lat, lon)
    out = tmp_path / "metrics.json"

    rc = sc.main(["--ens", str(m1), str(m2), "--truth", str(truth_ctl),
                  "--var", "tvir", "--level", "2", "--region", "nh",
                  "--leads", "1,2", "--out", str(out)])
    assert rc == 0
    stdout = capsys.readouterr().out
    assert "MB09 score" in stdout and "lead" in stdout

    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["n_members"] == 2
    assert payload["region"] == "nh"
    assert payload["variable"] == "tvir"
    assert payload["level"] == 2
    assert payload["leads"] == [1, 2]
    assert set(payload["metrics"]) == {"1", "2"}
    for lead in ("1", "2"):
        m = payload["metrics"][lead]
        assert set(m) >= {"rmse", "spread", "pac", "bss", "peca_sin", "peca_opt"}
        assert m["n_points"] == 4                 # so a linha de 45N
        assert m["n_members"] == 2
        assert np.isfinite(m["rmse"])
        assert set(m["bss"]) >= {"bss", "rel", "res", "unc", "n"}

    # confere o RMSE do primeiro lead contra o calculo direto (nivel 2 = il 1,
    # regiao NH = ultima linha da grade, 45N)
    mean = 0.5 * (members[0]["tvir"][0] + members[1]["tvir"][0])
    esperado = float(np.sqrt(np.mean((mean[1, 2, :] - truth["tvir"][0, 1, 2, :]) ** 2)))
    # o binario e float32: tolerancia absoluta de 1e-4
    assert payload["metrics"]["1"]["rmse"] == pytest.approx(esperado, abs=1e-4)


def test_cli_sem_out_so_stdout(tmp_path, capsys):
    lat = np.array([-45.0, 0.0, 45.0])
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    truth, members = _case_arrays()
    truth_ctl = _write_case(tmp_path / "t", truth, lat, lon)
    m1 = _write_case(tmp_path / "a", members[0], lat, lon)
    m2 = _write_case(tmp_path / "b", members[1], lat, lon)
    rc = sc.main(["--ens", str(m1), str(m2), "--truth", str(truth_ctl),
                  "--var", "pslc", "--region", "global"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "MB09 score" in out
    linhas = [ln for ln in out.splitlines() if ln.strip().split()[:1] in (["1"], ["2"])]
    assert len(linhas) == 2
    assert all("nan" not in ln.split()[1] for ln in linhas)   # rmse finito


def test_cli_erros_de_entrada(tmp_path, capsys):
    lat = np.array([-45.0, 0.0, 45.0])
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    truth, members = _case_arrays()
    truth_ctl = _write_case(tmp_path / "t", truth, lat, lon)
    m1 = _write_case(tmp_path / "a", members[0], lat, lon)

    assert sc.main(["--ens", str(m1), "--truth", str(truth_ctl),
                    "--var", "tvir", "--leads", "9"]) == 1
    assert "fora de" in capsys.readouterr().err
    assert sc.main(["--ens", str(m1), "--truth", str(truth_ctl),
                    "--var", "naoexiste"]) == 1
    assert "ausente" in capsys.readouterr().err
    assert sc.main(["--ens", str(m1), "--truth", str(truth_ctl),
                    "--var", "tvir", "--leads", "1,1"]) == 1
    assert "repetidos" in capsys.readouterr().err


def test_cli_grade_incompativel(tmp_path, capsys):
    lat = np.array([-45.0, 0.0, 45.0])
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    truth, members = _case_arrays()
    truth_ctl = _write_case(tmp_path / "t", truth, lat, lon)
    m1 = _write_case(tmp_path / "a", members[0], lat, lon)
    m2 = _write_case(tmp_path / "b", members[1], -lat, lon)
    rc = sc.main(["--ens", str(m1), str(m2), "--truth", str(truth_ctl),
                  "--var", "tvir"])
    assert rc == 1
    assert "grades diferentes" in capsys.readouterr().err
