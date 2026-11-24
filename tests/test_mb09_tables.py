"""Tabelas do MB09 (Mendonca & Bonatti 2009): caixas, presets, amplitudes, Q.

Verificam os valores transcritos do artigo (Tabela 1, p. 1442; Tabela 2 e
regioes, p. 1443) e as convencoes documentadas em docs/mb09.md. Nada aqui
depende de dados externos, exceto os arquivos DeltaSigma.* do proprio pacote.
"""

from pathlib import Path

import numpy as np
import pytest

from bamci import mb09_tables as mb

DELTASIGMA_DIR = Path(__file__).resolve().parents[1] / "DeltaSigma"

#: Tabela 1 do artigo, transcrita de forma independente do modulo (x1e-3).
TABELA_1_ARTIGO = [
    0.77, 0.78, 0.78, 0.78, 0.80, 0.82, 0.88, 0.98,
    1.14, 1.27, 1.37, 1.35, 1.18, 1.05, 0.90, 0.75,
    0.49, 0.26, 0.12, 0.05, 0.02, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00,
]


def _delta(name: str):
    return [float(t) for t in (DELTASIGMA_DIR / name).read_text().split()]


# ---------------------------------------------------------------------------
# Regioes
# ---------------------------------------------------------------------------
def test_seis_regioes_na_ordem_do_artigo():
    assert list(mb.REGIONS) == ["NH", "SH", "TR", "ETR", "NSA", "SSA"]
    assert mb.REGIONS_ORDER == ("NH", "SH", "TR", "ETR", "NSA", "SSA")
    assert set(mb.REGION_LABELS) == set(mb.REGIONS)


@pytest.mark.parametrize(
    "nome, esperado",
    [
        ("NH", (20.0, 90.0, 0.0, 360.0)),
        ("SH", (-90.0, -20.0, 0.0, 360.0)),
        ("TR", (-20.0, 20.0, 0.0, 360.0)),
        ("ETR", (-45.0, 30.0, 0.0, 360.0)),
        ("NSA", (-20.0, 20.0, 260.0, 350.0)),
        ("SSA", (-60.0, -20.0, 250.0, 340.0)),
    ],
)
def test_caixas_exatas(nome, esperado):
    assert mb.REGIONS[nome] == esperado


@pytest.mark.parametrize("nome", list(mb.REGIONS))
def test_limites_fisicos(nome):
    lat_min, lat_max, lon_min, lon_max = mb.REGIONS[nome]
    assert -90.0 <= lat_min < lat_max <= 90.0
    assert 0.0 <= lon_min <= 360.0 and 0.0 <= lon_max <= 360.0
    # caixas do MB09 nao cruzam 0/360 em 0-360 (100W-10W e 110W-20W)
    assert lon_min <= lon_max


# ---------------------------------------------------------------------------
# Presets (Tabela 2, p. 1443)
# ---------------------------------------------------------------------------
def test_presets_do_artigo():
    assert mb.PRESETS["oper"] == (("ETR",), ("T", "U", "V"))
    assert mb.PRESETS["ext1"] == (("NH", "SH", "TR"), ("T", "U", "V"))
    assert mb.PRESETS["trop"] == (("ETR",), ("P", "T", "Q", "U", "V"))
    assert mb.PRESETS["ext2"] == (("NH", "SH", "TR"), ("P", "T", "Q", "U", "V"))
    assert mb.PRESETS["etsa"] == (
        ("NH", "SH", "TR", "NSA", "SSA"),
        ("P", "T", "Q", "U", "V"),
    )


def test_preset_default_e_etsa():
    assert mb.DEFAULT_PRESET == "etsa"
    assert mb.preset() == mb.PRESETS["etsa"]
    assert mb.preset(None) == mb.PRESETS["etsa"]


@pytest.mark.parametrize("nome", list(mb.PRESETS))
def test_presets_consistentes(nome):
    regioes, campos = mb.PRESETS[nome]
    assert set(regioes) <= set(mb.REGIONS)
    assert len(set(regioes)) == len(regioes)
    assert set(campos) <= set(mb.FIELDS)
    assert mb.regions_of(nome) == regioes
    assert mb.fields_of(nome) == campos


def test_oper_e_ext1_nao_perturbam_p_e_q():
    for nome in ("oper", "ext1"):
        regioes, campos = mb.PRESETS[nome]
        assert "P" not in campos and "Q" not in campos


def test_preset_inexistente():
    with pytest.raises(KeyError):
        mb.preset("mb09-inexistente")


# ---------------------------------------------------------------------------
# Amplitudes (p. 1441-1442)
# ---------------------------------------------------------------------------
def test_amplitudes():
    assert mb.AMPLITUDES_RANDOM == {"wind": 3.0, "T": 0.6, "P": 1.0}
    assert mb.AMPLITUDES_TARGET == {"wind": 5.0, "T": 1.5, "P": 1.0}
    assert mb.AMPLITUDES["random"] == mb.AMPLITUDES_RANDOM
    assert mb.AMPLITUDES["target"] == mb.AMPLITUDES_TARGET
    assert mb.AMPLITUDE_UNITS == {"wind": "m/s", "T": "K", "P": "hPa", "Q": "kg/kg"}


def test_amplitude_alvo_maior_que_random():
    for chave, valor in mb.AMPLITUDES_RANDOM.items():
        assert mb.AMPLITUDES_TARGET[chave] >= valor


def test_constantes_temporais():
    assert mb.N_PERTURBATIONS == 7
    assert mb.N_MEMBERS == 15
    assert mb.DISCARD_HOURS == 6
    assert mb.FORECAST_HOURS == 36
    assert mb.OUTPUT_INTERVAL_HOURS == 3
    assert mb.FORECAST_TIMES_HOURS[0] == 6
    assert mb.FORECAST_TIMES_HOURS[-1] == 36
    assert mb.M_TIMES == 11


# ---------------------------------------------------------------------------
# Perfil de Q
# ---------------------------------------------------------------------------
def test_tabela1_transcrita():
    assert len(mb.Q_STDDEV_L28_X1E3) == 28
    assert list(mb.Q_STDDEV_L28_X1E3) == TABELA_1_ARTIGO


def test_tabela1_escala_1e3():
    assert len(mb.Q_STDDEV_L28) == 28
    for bruto, fisico in zip(mb.Q_STDDEV_L28_X1E3, mb.Q_STDDEV_L28):
        assert fisico == pytest.approx(bruto * 1.0e-3, rel=0.0, abs=0.0)
    assert mb.Q_STDDEV_L28[10] == pytest.approx(1.37e-3)
    assert max(mb.Q_STDDEV_L28) == pytest.approx(1.37e-3)


def test_tabela1_nao_negativa_e_unimodal():
    valores = np.asarray(mb.Q_STDDEV_L28)
    assert np.all(valores >= 0.0)
    pico = int(np.argmax(valores))
    assert pico == 10                        # sigma 11, ordem superficie->topo
    assert np.all(np.diff(valores[:pico + 1]) >= 0.0)     # sobe ate o pico
    assert np.all(np.diff(valores[pico:]) <= 0.0)         # decresce depois


def test_tabela1_zeros_no_topo():
    valores = np.asarray(mb.Q_STDDEV_L28)
    assert np.all(valores[21:] == 0.0)       # sigma 22..28
    assert valores[20] > 0.0                 # sigma 21 = 0.02


def test_q28_reproduz_a_tabela():
    q28 = np.asarray(mb.q_stddev_for_layers(28))
    assert q28.shape == (28,)
    assert np.allclose(q28, mb.Q_STDDEV_L28, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("n", [9, 18, 28, 42, 64])
def test_q_para_outras_resolucoes(n):
    q = np.asarray(mb.q_stddev_for_layers(n))
    assert q.shape == (n,)
    assert np.all(q >= 0.0)
    assert np.all(q <= max(mb.Q_STDDEV_L28) + 1e-15)
    pico = int(np.argmax(q))
    assert np.all(np.diff(q[:pico + 1]) >= 0.0)
    assert np.all(np.diff(q[pico:]) <= 0.0)
    assert q[-1] == 0.0                      # topo seco


def test_q42_zeros_no_topo_e_base_preservada():
    q42 = np.asarray(mb.q_stddev_for_layers(42))
    assert q42[0] == pytest.approx(0.77e-3)  # camada 1 = superficie
    assert np.all(q42[33:] == 0.0)           # topo (a partir de sigma ~0.78)


def test_q_com_deltasigma_real():
    d28 = _delta("DeltaSigma.L028")
    d42 = _delta("DeltaSigma.L042")
    q42 = np.asarray(mb.q_stddev_for_layers(42, delta_sigma=d42, source_delta_sigma=d28))
    assert q42.shape == (42,)
    assert q42[-1] == 0.0
    assert 1.3e-3 < q42.max() < 1.4e-3
    uniforme = np.asarray(mb.q_stddev_for_layers(42))
    # a aproximacao uniforme e proxima, mas nao identica (sigma real != uniforme)
    assert np.max(np.abs(q42 - uniforme)) < 5e-4
    assert not np.allclose(q42, uniforme)


def test_q_com_sigma_layers_explicito():
    sigma42 = mb.sigma_centers_from_delta(_delta("DeltaSigma.L042"))
    q = mb.q_stddev_for_layers(42, sigma_layers=sigma42)
    assert isinstance(q, list) and len(q) == 42
    via_delta = mb.q_stddev_for_layers(42, delta_sigma=_delta("DeltaSigma.L042"))
    assert np.allclose(q, via_delta, rtol=0.0, atol=0.0)


def test_q_precedencia_sigma_layers_sobre_delta_sigma():
    sigma_uniforme = mb.sigma_centers_uniform(42)
    uniforme = mb.q_stddev_for_layers(42, sigma_layers=sigma_uniforme)
    via_delta = mb.q_stddev_for_layers(42, delta_sigma=_delta("DeltaSigma.L042"))
    ambos = mb.q_stddev_for_layers(42, sigma_uniforme, delta_sigma=_delta("DeltaSigma.L042"))
    assert np.allclose(ambos, uniforme, rtol=0.0, atol=0.0)
    assert not np.allclose(ambos, via_delta)


def test_q_erros_de_argumento():
    with pytest.raises(ValueError):
        mb.q_stddev_for_layers(0)
    with pytest.raises(ValueError):
        mb.q_stddev_for_layers(42, sigma_layers=[0.1] * 41)
    with pytest.raises(ValueError):
        mb.q_stddev_for_layers(42, delta_sigma=[0.1] * 41)
    with pytest.raises(ValueError):
        mb.q_stddev_for_layers(28, source_delta_sigma=[0.1] * 20)
    with pytest.raises(ValueError):
        mb.q_stddev_for_layers(28, source_sigma_layers=[0.1] * 20)
    with pytest.raises(ValueError):
        # sigma_layers invertido (topo -> superficie)
        mb.q_stddev_for_layers(3, sigma_layers=[0.9, 0.5, 0.1])
    with pytest.raises(ValueError):
        mb.sigma_centers_from_delta([0.5, -0.5])


# ---------------------------------------------------------------------------
# Sigma
# ---------------------------------------------------------------------------
def test_sigma_centros_uniformes():
    assert np.allclose(mb.sigma_centers_uniform(2), [0.25, 0.75])
    s42 = mb.sigma_centers_uniform(42)
    assert s42[0] == pytest.approx(0.5 / 42)
    assert s42[-1] == pytest.approx(41.5 / 42)
    assert np.all(np.diff(s42) > 0.0)


def test_sigma_centros_de_deltasigma():
    d = _delta("DeltaSigma.L042")
    s = mb.sigma_centers_from_delta(d)
    assert s.shape == (42,)
    assert np.all(np.diff(s) > 0.0)
    assert s[0] == pytest.approx(d[0] / 2.0)
    assert s[-1] == pytest.approx(1.0 - d[-1] / 2.0)
    # soma normalizada: espessuras dobradas dao os mesmos centros
    dobrado = mb.sigma_centers_from_delta([2.0 * v for v in d])
    assert np.allclose(s, dobrado, rtol=0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Mascaras
# ---------------------------------------------------------------------------
def test_mascara_pontos_tipicos():
    # (lat, lon, regiao, esperado)
    casos = [
        (45.0, 200.0, "NH", True),
        (45.0, 200.0, "SH", False),
        (-45.0, 300.0, "SH", True),
        (0.0, 0.0, "TR", True),
        (10.0, 350.0, "TR", True),
        (-40.0, 180.0, "ETR", True),
        (-50.0, 180.0, "ETR", False),
        (35.0, 180.0, "ETR", False),
        (0.0, 300.0, "NSA", True),
        (0.0, 260.0, "NSA", False),      # 100W: borda inicial meio-aberta
        (0.0, 259.0, "NSA", False),
        (0.0, 250.0, "SSA", False),
        (-40.0, 300.0, "SSA", True),
        (-40.0, 200.0, "SSA", False),
        (-40.0, 240.0, "SSA", False),
    ]
    for lat, lon, regiao, esperado in casos:
        obtido = bool(mb.region_mask(lat, lon, regiao))
        assert obtido is esperado, (lat, lon, regiao, obtido, esperado)


def test_mascara_bordas_meio_abertas():
    """Convencao D11: linha de grade na emenda interna pertence ao maximo."""
    # NH/TR se encontram em 20N: 20 pertence a TR, nao a NH
    assert bool(mb.region_mask(20.0, 0.0, "TR"))
    assert not bool(mb.region_mask(20.0, 0.0, "NH"))
    assert bool(mb.region_mask(20.1, 0.0, "NH"))
    # TR/SH se encontram em 20S: -20 pertence a SH, nao a TR
    assert bool(mb.region_mask(-20.0, 0.0, "SH"))
    assert not bool(mb.region_mask(-20.0, 0.0, "TR"))
    assert not bool(mb.region_mask(20.1, 0.0, "TR"))
    # limite de dominio fechado: polos e Greenwich nao podem ficar descobertos
    assert bool(mb.region_mask(90.0, 0.0, "NH"))
    assert not bool(mb.region_mask(90.1, 0.0, "NH"))
    assert bool(mb.region_mask(-90.0, 0.0, "SH"))
    assert bool(mb.region_mask(45.0, 0.0, "NH"))
    assert bool(mb.region_mask(0.0, 0.0, "TR"))
    # lon limite de NSA: 10W = 350 em 0-360 (borda final, fechada)
    assert bool(mb.region_mask(0.0, 350.0, "NSA"))
    assert not bool(mb.region_mask(0.0, 350.1, "NSA"))
    assert bool(mb.region_mask(0.0, -10.0, "NSA"))            # -10 == 350
    assert bool(mb.region_mask(0.0, 360.0, "NSA")) is False   # 360 -> 0


def test_mascara_sem_duplicata_e_sem_buraco_na_lat():
    lat = np.arange(-90.0, 90.0001, 0.25)
    lon = np.zeros_like(lat)
    total = np.zeros_like(lat, dtype=int)
    for regiao in ("NH", "SH", "TR"):
        total += mb.region_mask(lat, lon, regiao).astype(int)
    assert np.all(total >= 1), "latitude sem cobertura"
    assert np.all(total <= 1), "latitude duplicada entre NH/SH/TR"
    # o equador e as emendas estao cobertos
    for valor in (-90.0, -20.0, 0.0, 20.0, 90.0):
        assert total[np.isclose(lat, valor)][0] == 1


def test_mascara_converte_longitude_negativa():
    assert bool(mb.region_mask(-40.0, -60.0, "SSA"))     # -60 == 300
    assert not bool(mb.region_mask(-40.0, -120.0, "SSA"))   # -120 == 240 -> fora
    assert mb.region_mask(-40.0, -60.0, "SSA") == mb.region_mask(-40.0, 300.0, "SSA")


def test_mascara_caixa_que_cruza_0_360():
    # caixa sintetica 300E-20E (cruza Greenwich); nenhuma do MB09 cruza
    regiao = (-10.0, 10.0, 300.0, 20.0)
    assert bool(mb.region_mask(0.0, 350.0, regiao))
    assert bool(mb.region_mask(0.0, 10.0, regiao))
    assert not bool(mb.region_mask(0.0, 300.0, regiao))   # borda inicial aberta
    assert not bool(mb.region_mask(0.0, 100.0, regiao))
    assert not bool(mb.region_mask(0.0, 250.0, regiao))


def test_mascara_vetorial_e_broadcast():
    lat = np.array([[-45.0, 0.0], [45.0, -80.0]])
    lon = np.array([[300.0, 200.0], [10.0, 180.0]])
    m = mb.region_mask(lat, lon, "SH")
    assert m.shape == (2, 2)
    assert m.tolist() == [[True, False], [False, True]]
    # lat escalar + lon vetorial
    m2 = mb.region_mask(0.0, np.array([0.0, 90.0, 180.0, 270.0, 360.0]), "TR")
    assert m2.tolist() == [True] * 5
    # latitude fora de qualquer caixa
    assert not mb.region_mask(0.0, 0.0, "NH")
    assert not mb.region_mask(-45.0, 0.0, "NH")


def test_mascara_regiao_invalida():
    with pytest.raises(KeyError):
        mb.region_mask(0.0, 0.0, "XX")
    with pytest.raises(ValueError):
        mb.region_mask(0.0, 0.0, (0.0, 10.0, 20.0))


def test_mascara_preset_etsa():
    lat = np.array([-40.0, 0.0, 45.0])
    lon = np.array([300.0, 300.0, 200.0])
    mascaras = mb.region_mask_preset(lat, lon)
    assert list(mascaras) == ["NH", "SH", "TR", "NSA", "SSA"]
    # sobreposicao proposital: o ponto (-40, 300) cai em SH, TR? nao (TR e
    # 20S-20N) -- cai em SH, SSA e ETR
    assert mascaras["SH"].tolist() == [True, False, False]
    assert mascaras["SSA"].tolist() == [True, False, False]
    assert mascaras["NH"].tolist() == [False, False, True]
    # (0, 300) cai em TR e NSA ao mesmo tempo (base do D9)
    assert mascaras["TR"].tolist() == [False, True, False]
    assert mascaras["NSA"].tolist() == [False, True, False]
    assert list(mb.region_mask_preset(lat, lon, "ext2")) == ["NH", "SH", "TR"]


# ---------------------------------------------------------------------------
# Perfis de Q: artigo (L28) e operacional (L42)
# ---------------------------------------------------------------------------
#: Perfil operacional transcrito de forma independente (umipert_eof.nml, g/kg).
Q_OPER_ARTIGO_NML = [
    0.770, 0.780, 0.780, 0.820, 0.800, 0.980, 1.030, 1.080, 1.130, 1.180,
    1.050, 0.786, 0.523, 0.260, 0.195, 0.130, 0.065,
] + [0.0] * 25


def test_perfil_operacional_42_valores_verbatim():
    assert len(mb.Q_STDDEV_OPER_L42_X1E3) == 42
    assert list(mb.Q_STDDEV_OPER_L42_X1E3) == Q_OPER_ARTIGO_NML


def test_perfil_operacional_zero_acima_da_camada_17():
    valores = np.asarray(mb.Q_STDDEV_OPER_L42_X1E3)
    assert np.all(valores[:17] > 0.0)          # camadas 1..17
    assert np.all(valores[17:] == 0.0)         # camadas 18..42
    assert valores[16] == pytest.approx(0.065)  # ultima camada perturbada
    assert valores[11] == pytest.approx(0.786)  # (g/kg, dividido por 1000 depois)
    assert valores[13] == pytest.approx(0.260)  # reaparece o 0.260 da Tabela 1


def test_perfil_operacional_escala_1e3():
    assert len(mb.Q_STDDEV_OPER_L42) == 42
    for bruto, fisico in zip(mb.Q_STDDEV_OPER_L42_X1E3, mb.Q_STDDEV_OPER_L42):
        assert fisico == pytest.approx(bruto * 1.0e-3, rel=0.0, abs=0.0)
    assert mb.Q_STDDEV_OPER_L42[0] == pytest.approx(7.7e-4)
    assert mb.Q_STDDEV_OPER_L42[17] == 0.0


def test_q_profiles_mapeia_os_dois_perfis():
    assert set(mb.Q_PROFILES) == {"artigo_l28", "operacional_l42"}
    assert mb.Q_PROFILES["artigo_l28"] is mb.Q_STDDEV_L28
    assert mb.Q_PROFILES["operacional_l42"] is mb.Q_STDDEV_OPER_L42
    assert mb.Q_PROFILE_DEFAULT == "artigo_l28"


def test_q_profile_artigo_e_o_default():
    base = mb.q_stddev_for_layers(42)
    assert base == mb.q_stddev_for_layers(42, profile="artigo_l28")
    # as duas opcoes divergem em L42: o operacional zera a partir da camada 18
    oper = mb.q_stddev_for_layers(42, profile="operacional_l42")
    assert oper != base
    assert np.all(np.asarray(oper[17:]) == 0.0)
    assert np.any(np.asarray(base[17:]) > 0.0)   # o artigo ainda perturba 18..33


def test_q_profile_operacional_devolve_verbatim():
    q = mb.q_stddev_for_layers(42, profile="operacional_l42")
    assert q == list(mb.Q_STDDEV_OPER_L42)
    assert len(q) == 42
    # aceita (e valida) as espessuras de L42 sem interpolar
    d42 = _delta("DeltaSigma.L042")
    assert mb.q_stddev_for_layers(42, profile="operacional_l42", delta_sigma=d42) == q


@pytest.mark.parametrize("n", [1, 9, 28, 41, 43, 64])
def test_q_profile_operacional_exige_42_camadas(n):
    with pytest.raises(ValueError) as exc:
        mb.q_stddev_for_layers(n, profile="operacional_l42")
    assert "42 camadas" in str(exc.value)
    assert "artigo_l28" in str(exc.value)


def test_q_profile_operacional_valida_sigma_de_42():
    with pytest.raises(ValueError):
        mb.q_stddev_for_layers(42, profile="operacional_l42",
                               sigma_layers=[0.5] * 41)


def test_q_profile_invalido():
    with pytest.raises(ValueError) as exc:
        mb.q_stddev_for_layers(42, profile="tabela1")
    assert "artigo_l28" in str(exc.value)
    # caixa/espacos sao normalizados: nao e erro
    assert (mb.q_stddev_for_layers(42, profile="ARTIGO_L28 ")
            == mb.q_stddev_for_layers(42))


# ---------------------------------------------------------------------------
# Regioes por indice de grade (operacional 384x192, j=1 = norte)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "regiao, esperado",
    [
        ("HN", (0, 75, 0, 384)),
        ("TR", (75, 117, 0, 384)),
        ("HS", (117, 192, 0, 384)),
        ("NAS", (74, 108, 275, 372)),
        ("SAS", (108, 161, 265, 362)),
    ],
)
def test_region_indices_exatos(regiao, esperado):
    assert mb.region_indices(192, 384, regiao) == esperado
    assert mb.region_indices(192, 384, regiao.lower()) == esperado


@pytest.mark.parametrize(
    "alias, operacional",
    [("NH", "HN"), ("SH", "HS"), ("TR", "TR"), ("NSA", "NAS"), ("SSA", "SAS")],
)
def test_region_indices_aceita_nome_do_artigo(alias, operacional):
    assert mb.region_indices(192, 384, alias) == mb.region_indices(192, 384, operacional)
    assert mb.region_indices(192, 384, alias.lower()) == mb.INDEX_REGIONS[operacional]


def test_region_indices_faixas_1based_conferidas():
    # J 1..75 / 76..117 / 118..192 e I 276..372 / 266..362 -> 0-based
    assert mb.INDEX_REGIONS["HN"] == (1 - 1, 75, 1 - 1, 384)
    assert mb.INDEX_REGIONS["TR"] == (76 - 1, 117, 1 - 1, 384)
    assert mb.INDEX_REGIONS["HS"] == (118 - 1, 192, 1 - 1, 384)
    assert mb.INDEX_REGIONS["NAS"] == (75 - 1, 108, 276 - 1, 372)
    assert mb.INDEX_REGIONS["SAS"] == (109 - 1, 161, 266 - 1, 362)


def test_region_indices_erros():
    with pytest.raises(KeyError):
        mb.region_indices(192, 384, "ETR")
    with pytest.raises(ValueError) as exc:
        mb.region_indices(384, 192, "HN")
    assert "192x384" in str(exc.value)
    with pytest.raises(ValueError):
        mb.region_indices(96, 192, "HN")
    with pytest.raises(ValueError):
        mb.region_mask_by_index(192, 383, "HN")


@pytest.mark.parametrize(
    "regiao, esperado",
    [
        ("HN", 384 * 75),
        ("TR", 384 * 42),
        ("HS", 384 * 75),
        ("NAS", 97 * 34),
        ("SAS", 97 * 53),
    ],
)
def test_region_mask_by_index_contagem(regiao, esperado):
    mask = mb.region_mask_by_index(192, 384, regiao)
    assert mask.shape == (192, 384)
    assert mask.dtype == bool
    assert int(mask.sum()) == esperado


def test_region_mask_by_index_e_retangular():
    mask = mb.region_mask_by_index(192, 384, "NAS")
    j0, j1, i0, i1 = mb.region_indices(192, 384, "NAS")
    assert np.all(mask[j0:j1, i0:i1])
    assert int(mask.sum()) == (j1 - j0) * (i1 - i0)
    # fora da caixa nao ha True
    fora = mask.copy()
    fora[j0:j1, i0:i1] = False
    assert not fora.any()


def _grid_latlon(nlat=192, nlon=384):
    """Latitudes/longitudes de celula (centros) da grade do operacional."""
    lat = 90.0 - (np.arange(nlat) + 0.5) * 180.0 / nlat      # j=1 -> norte
    lon = (np.arange(nlon) + 0.5) * 360.0 / nlon
    return lat, lon


@pytest.mark.parametrize("op, artigo", [("HN", "NH"), ("TR", "TR"), ("HS", "SH")])
def test_mascara_por_indice_identica_a_latlon_nas_bandas(op, artigo):
    """As 3 faixas (HN/TR/HS) reproduzem exatamente as caixas do artigo."""
    lat, lon = _grid_latlon()
    idx = mb.region_mask_by_index(192, 384, op)
    ll = mb.region_mask(lat[:, None], lon[None, :], artigo)
    assert np.array_equal(idx, ll)
    # as faixas sao disjuntas e cobrem todas as latitudes
    faixas = {op: mb.region_mask_by_index(192, 384, op)[:, 0] for op in ("HN", "TR", "HS")}
    assert sum(int(v.sum()) for v in faixas.values()) == len(lat)
    assert not (faixas["HN"] & faixas["TR"]).any()
    assert not (faixas["TR"] & faixas["HS"]).any()


def test_mascara_por_indice_longitudes_dentro_de_2_graus():
    """NAS/SAS: longitudes dentro de 2 graus das caixas do artigo (lat nao)."""
    lat, lon = _grid_latlon()
    casos = [("NAS", "NSA", 260.0, 350.0), ("SAS", "SSA", 250.0, 340.0)]
    for op, artigo, lon_min, lon_max in casos:
        j0, j1, i0, i1 = mb.region_indices(192, 384, op)
        assert abs(lon[i0] - (lon_min + 360.0) % 360.0) <= 2.0
        assert abs(lon[i1 - 1] - lon_max) <= 2.0
        # a mascara por indice fica dentro da caixa do artigo dilatada em 2 graus
        dilatada = mb.region_mask(lat[:, None], lon[None, :],
                                  (-90.0, 90.0, lon_min - 2.0, lon_max + 2.0))
        idx = mb.region_mask_by_index(192, 384, op)
        assert np.all(dilatada[idx])


def test_caixas_por_indice_do_operacional_nao_sao_as_do_artigo_em_latitude():
    """Em latitude NAS/SAS do operacional diferem >2 graus das caixas do artigo.

    Documenta a divergencia (nao e bug): J 75..108 cobre 20.2N-10.8S (o artigo
    pede 20S-20N) e J 109..161 cobre 11.7S-60.5S (o artigo pede 20S-60S).
    """
    lat, _ = _grid_latlon()
    j0, j1, _, _ = mb.region_indices(192, 384, "NAS")
    assert abs(lat[j1 - 1] - (-20.0)) > 2.0        # ~9.2 graus
    assert abs(lat[j0] - 20.0) <= 2.0
    j0, j1, _, _ = mb.region_indices(192, 384, "SAS")
    assert abs(lat[j0] - (-20.0)) > 2.0            # ~8.3 graus
    assert abs(lat[j1 - 1] - (-60.0)) <= 2.0

