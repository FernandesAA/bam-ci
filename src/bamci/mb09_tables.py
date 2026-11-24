"""Tabelas do metodo de perturbacao EOF de Mendonca & Bonatti (2009).

Constantes puras do MB09 (Mon. Wea. Rev., 137(4), 1438-1459,
DOI 10.1175/2008MWR2581.1) usadas para gerar condicoes iniciais perturbadas:
caixas de perturbacao, presets de experimento, amplitudes prescritas e o perfil
vertical de desvio-padrao da umidade especifica (Tabela 1, p. 1442).

Especificacao completa, com as paginas do artigo e a marcacao das decisoes de
implementacao: 'docs/mb09.md'.

Convencoes (ver docs/mb09.md, secao 7):

* sigma = coordenada de camada; camada 1 e a MAIS PROXIMA DA SUPERFICIE e sigma
  cresce para o topo (D1 -- o artigo nao declara o sentido);
* longitudes em graus, 0 <= lon < 360, 0 = Greenwich;
* bordas das caixas: MEIO-ABERTAS nas emendas internas (para nao duplicar a
  linha de grade na fronteira entre NH e TR, ja que ETSA soma sem renormalizar)
  e FECHADAS no limite do dominio (lat -90/+90, lon 0);
* unidades: vento m/s, temperatura K, pressao de superficie hPa, umidade
  especifica kg/kg.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

# ---------------------------------------------------------------------------
# 1. Regioes de perturbacao (MB09, p. 1443)
# ---------------------------------------------------------------------------
#: nome -> (lat_min, lat_max, lon_min_360, lon_max_360)
#:
#: Texto do artigo: NH 20N-90N, SH 20S-90S, TR 20S-20N, ETR 45S-30N,
#: NSA 20S-20N/100W-10W, SSA 60S-20S/110W-20W. Todas as caixas tem
#: lon_min <= lon_max em 0-360 (100W-10W = 260-350; 110W-20W = 250-340);
#: nenhuma cruza 0/360, mas region_mask suporta caixas que cruzem.
REGIONS: Dict[str, Tuple[float, float, float, float]] = {
    "NH": (20.0, 90.0, 0.0, 360.0),
    "SH": (-90.0, -20.0, 0.0, 360.0),
    "TR": (-20.0, 20.0, 0.0, 360.0),
    "ETR": (-45.0, 30.0, 0.0, 360.0),
    "NSA": (-20.0, 20.0, 260.0, 350.0),
    "SSA": (-60.0, -20.0, 250.0, 340.0),
}

#: Ordem canonica das caixas: a ordem em que aparecem na p. 1443.
REGIONS_ORDER: Tuple[str, ...] = ("NH", "SH", "TR", "ETR", "NSA", "SSA")

#: Descricao textual de cada caixa (mesma redacao do artigo, p. 1443).
REGION_LABELS: Dict[str, str] = {
    "NH": "Northern Hemisphere: 20N-90N, 0-360",
    "SH": "Southern Hemisphere: 20S-90S, 0-360",
    "TR": "tropics: 20S-20N, 0-360",
    "ETR": "extended tropical region: 45S-30N, 0-360",
    "NSA": "northern South America: 20S-20N, 100W-10W",
    "SSA": "southern South America: 60S-20S, 110W-20W",
}

# ---------------------------------------------------------------------------
# 2. Experimentos (Tabela 2, p. 1443)
# ---------------------------------------------------------------------------
#: Campos perturbaveis, na ordem usada na Tabela 2.
FIELDS: Tuple[str, ...] = ("P", "T", "Q", "U", "V")

#: preset -> (regioes na ordem de combinacao, campos perturbados)
PRESETS: Dict[str, Tuple[Tuple[str, ...], Tuple[str, ...]]] = {
    "oper": (("ETR",), ("T", "U", "V")),
    "ext1": (("NH", "SH", "TR"), ("T", "U", "V")),
    "trop": (("ETR",), ("P", "T", "Q", "U", "V")),
    "ext2": (("NH", "SH", "TR"), ("P", "T", "Q", "U", "V")),
    "etsa": (("NH", "SH", "TR", "NSA", "SSA"), ("P", "T", "Q", "U", "V")),
}

#: Configuracao padrao (D10 -- o artigo nao elege uma; ETSA e a mais completa).
DEFAULT_PRESET = "etsa"

# ---------------------------------------------------------------------------
# 3. Amplitudes prescritas (p. 1441-1442)
# ---------------------------------------------------------------------------
#: Desvio-padrao das perturbacoes aleatorias iniciais (passo 1, p. 1441).
#: Vento em m/s, T em K, P (pressao de superficie) em hPa. Q usa a Tabela 1.
AMPLITUDES_RANDOM: Dict[str, float] = {"wind": 3.0, "T": 0.6, "P": 1.0}

#: Desvio-padrao alvo da reescala (passo 5, p. 1442). ZK99 usa 3/0.6; o CPTEC
#: adota 5.0 m/s e 1.5 K (Coutinho 1999; Daley & Mayer 1986) e 1.0 hPa para P
#: (Anderson et al. 2005).
AMPLITUDES_TARGET: Dict[str, float] = {"wind": 5.0, "T": 1.5, "P": 1.0}

#: Unidades das amplitudes acima (Q em kg/kg, perfil da Tabela 1).
AMPLITUDE_UNITS: Dict[str, str] = {"wind": "m/s", "T": "K", "P": "hPa", "Q": "kg/kg"}

#: Atalho: AMPLITUDES["random"] / AMPLITUDES["target"].
AMPLITUDES: Dict[str, Dict[str, float]] = {
    "random": dict(AMPLITUDES_RANDOM),
    "target": dict(AMPLITUDES_TARGET),
}

#: Constantes temporais do metodo (p. 1441).
N_PERTURBATIONS = 7          # n; 2n = 14 estados + controle = 15 membros
N_MEMBERS = 2 * N_PERTURBATIONS + 1
DISCARD_HOURS = 6            # primeiras 6 h descartadas (auto-ajuste)
FORECAST_HOURS = 36          # integracao curta
OUTPUT_INTERVAL_HOURS = 3    # saida a cada 3 h
FORECAST_TIMES_HOURS: Tuple[int, ...] = tuple(
    range(DISCARD_HOURS, FORECAST_HOURS + 1, OUTPUT_INTERVAL_HOURS)
)
M_TIMES = len(FORECAST_TIMES_HOURS)   # M = 11

# ---------------------------------------------------------------------------
# 4. Perfil de umidade especifica (Tabela 1, p. 1442)
# ---------------------------------------------------------------------------
#: Tabela 1 do MB09 exatamente como impressa (fatores x1e-3), ordem
#: superficie -> topo (D1), camada sigma 1..28.
Q_STDDEV_L28_X1E3: Tuple[float, ...] = (
    0.77, 0.78, 0.78, 0.78, 0.80, 0.82, 0.88, 0.98,
    1.14, 1.27, 1.37, 1.35, 1.18, 1.05, 0.90, 0.75,
    0.49, 0.26, 0.12, 0.05, 0.02, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00,
)

#: O mesmo perfil em kg/kg (a tabela publicada traz x1e-3).
Q_STDDEV_L28: Tuple[float, ...] = tuple(1.0e-3 * v for v in Q_STDDEV_L28_X1E3)

#: Perfil de Q da IMPLEMENTACAO OPERACIONAL (oensMB09_bam_xc50), ja em 42 camadas
#: sigma do modelo (nao e a Tabela 1 do artigo, que e L28). Fonte:
#: 'include/TQ0126L042/umipert_eof.nml' (STDQ em g/kg); a conversao para kg/kg
#: esta em 'eofhum.f:235-236' (divide por 1000). Camadas 1..17 tem valor; de 18
#: a 42 o operacional zera (topo nao perturbado).
Q_STDDEV_OPER_L42_X1E3: Tuple[float, ...] = (
    0.770, 0.780, 0.780, 0.820, 0.800, 0.980, 1.030, 1.080, 1.130, 1.180,
    1.050, 0.786, 0.523, 0.260, 0.195, 0.130, 0.065,
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.000, 0.000,
)

#: O mesmo perfil operacional em kg/kg.
Q_STDDEV_OPER_L42: Tuple[float, ...] = tuple(1.0e-3 * v for v in Q_STDDEV_OPER_L42_X1E3)

#: Perfis de desvio-padrao de Q disponiveis (kg/kg).
#:
#: * 'artigo_l28'       -- Tabela 1 do MB09 (p. 1442), 28 camadas sigma;
#:                         interpolado para as camadas do modelo (default);
#: * 'operacional_l42'  -- perfil da implementacao operacional em 42 camadas
#:                         sigma (umipert_eof.nml); so vale para n_layers=42.
Q_PROFILES: Dict[str, Tuple[float, ...]] = {
    "artigo_l28": Q_STDDEV_L28,
    "operacional_l42": Q_STDDEV_OPER_L42,
}

#: Perfil default de q (o do artigo; ver docs/mb09.md secao 4).
Q_PROFILE_DEFAULT = "artigo_l28"


# ---------------------------------------------------------------------------
# 5. Utilitarios de sigma
# ---------------------------------------------------------------------------
def sigma_centers_uniform(n_layers: int) -> np.ndarray:
    """Centros de camada em sigma para n_layers camadas uniformes.

    sigma_k = (k - 0.5) / n, com k = 1..n: camada 1 junto a superficie, sigma
    crescente para o topo (D1, ver docs/mb09.md). Aproximacao: o artigo nao da
    as coordenadas sigma (a Tabela 1 so traz os indices 1..28).
    """
    n = int(n_layers)
    if n <= 0:
        raise ValueError("n_layers deve ser > 0 (recebido %r)" % (n_layers,))
    k = np.arange(1, n + 1, dtype=np.float64)
    return (k - 0.5) / n


def sigma_centers_from_delta(delta_sigma: Sequence[float]) -> np.ndarray:
    """Centros de camada sigma a partir das espessuras DeltaSigma.Lnnn.

    Interfaces [0, d1, d1+d2, ..., 1]; centro = media das interfaces (mesma
    convencao de bamci.levels.sigma_layers). A soma e normalizada para 1 quando
    o arquivo traz residuo de arredondamento.
    """
    d = np.asarray(delta_sigma, dtype=np.float64).ravel()
    if d.size == 0:
        raise ValueError("delta_sigma vazio")
    if np.any(d <= 0.0):
        raise ValueError("delta_sigma tem espessura nao positiva")
    total = float(d.sum())
    if total <= 0.0:
        raise ValueError("delta_sigma soma %r" % (total,))
    d = d / total
    interfaces = np.concatenate([[0.0], np.cumsum(d)])
    return 0.5 * (interfaces[:-1] + interfaces[1:])


def q_stddev_for_layers(
    n_layers: int,
    sigma_layers: Optional[Sequence[float]] = None,
    *,
    source_sigma_layers: Optional[Sequence[float]] = None,
    delta_sigma: Optional[Sequence[float]] = None,
    source_delta_sigma: Optional[Sequence[float]] = None,
    profile: str = Q_PROFILE_DEFAULT,
) -> List[float]:
    """Perfil de desvio-padrao de Q (kg/kg) para as camadas do modelo.

    Dois perfis (ver Q_PROFILES e docs/mb09.md secao 4):

    * profile='artigo_l28' (DEFAULT): interpola linearmente a Tabela 1 do MB09
      (28 camadas sigma, p. 1442) nas coordenadas sigma do modelo -- a mesma
      ideia que os autores usaram ao interpolar os valores da ECMWF para as
      camadas do CPTEC (p. 1442);
    * profile='operacional_l42': devolve o perfil da implementacao operacional
      (umipert_eof.nml, 42 camadas sigma) SEM interpolar; exige n_layers == 42.

    Regra de ouro (perfil do artigo): a amplitude de q tem de sair do DeltaSigma
    do modelo. Sempre que o chamador puder, passe as coordenadas reais
    (sigma_layers=... ou delta_sigma=...) tanto para o destino quanto para a
    origem (Tabela 1 em L28). A grade uniforme e apenas fallback, porque a
    diferenca chega a ~30% num ponto da troposfera alta (docs/mb09.md, secao 4).

    Parametros
    ----------
    n_layers:
        numero de camadas do modelo (ex.: 42).
    sigma_layers:
        centros de camada sigma do modelo, tamanho n_layers, ordem
        superficie -> topo (ex.: bamci.mb09_tables.sigma_centers_from_delta(
        bamci.levels.read_delta_sigma("DeltaSigma/DeltaSigma.L042"))).
        Vence delta_sigma quando os dois sao informados.
    source_sigma_layers:
        centros de camada sigma da origem (28 valores). Vence
        source_delta_sigma.
    delta_sigma:
        alternativa a sigma_layers: espessuras sigma do modelo (ex.:
        DeltaSigma.L042). Converte para centros com sigma_centers_from_delta.
    source_delta_sigma:
        alternativa a source_sigma_layers: espessuras da origem (DeltaSigma.L028).
    profile:
        'artigo_l28' (default) ou 'operacional_l42'. Os argumentos de sigma so
        sao usados para validar o tamanho quando o perfil e o operacional.

    Retorna
    -------
    lista de n_layers floats em kg/kg, na ordem superficie -> topo.
    """
    n = int(n_layers)
    if n <= 0:
        raise ValueError("n_layers deve ser > 0 (recebido %r)" % (n_layers,))

    key = str(profile).strip().lower()
    if key not in Q_PROFILES:
        raise ValueError(
            "profile %r desconhecido; use um de %s" % (profile, sorted(Q_PROFILES))
        )
    if key == "operacional_l42":
        if n != 42:
            raise ValueError(
                "profile='operacional_l42' e o perfil da implementacao "
                "operacional em 42 camadas sigma (include/TQ0126L042/"
                "umipert_eof.nml, kg/kg em eofhum.f:235-236); ele nao pode ser "
                "interpolado para %d camadas. Use n_layers=42 ou "
                "profile='artigo_l28' (Tabela 1 em L28, interpolada)." % n
            )
        # valida (sem usar) as coordenadas passadas: precisam ser de L42
        _resolve_sigma(
            42,
            sigma_layers=sigma_layers,
            delta_sigma=delta_sigma,
            origem="destino (perfil operacional L42)",
        )
        return [float(v) for v in Q_STDDEV_OPER_L42]

    src_vals = np.asarray(Q_PROFILES[key], dtype=np.float64)
    src_sigma = _resolve_sigma(
        src_vals.size,
        sigma_layers=source_sigma_layers,
        delta_sigma=source_delta_sigma,
        origem="origem (Tabela 1)",
    )

    tgt_sigma = _resolve_sigma(
        n, sigma_layers=sigma_layers, delta_sigma=delta_sigma, origem="destino"
    )

    values = np.interp(tgt_sigma, src_sigma, src_vals)
    # np.interp nao extrapola, mas o clip garante o intervalo fisico [0, max].
    values = np.clip(values, 0.0, float(src_vals.max()))
    return [float(v) for v in values]


def _resolve_sigma(
    n_expected: int,
    *,
    sigma_layers: Optional[Sequence[float]],
    delta_sigma: Optional[Sequence[float]],
    origem: str,
) -> np.ndarray:
    """Centros de camada: sigma_layers > delta_sigma > grade uniforme."""
    if sigma_layers is not None:
        sigma = np.asarray(sigma_layers, dtype=np.float64).ravel()
        if sigma.size != n_expected:
            raise ValueError(
                "sigma_layers de %s tem %d camadas; esperado %d"
                % (origem, sigma.size, n_expected)
            )
        if np.any(np.diff(sigma) <= 0.0):
            raise ValueError(
                "sigma_layers de %s deve ser estritamente crescente" % (origem,)
            )
        return sigma
    if delta_sigma is not None:
        sigma = sigma_centers_from_delta(delta_sigma)
        if sigma.size != n_expected:
            raise ValueError(
                "delta_sigma de %s tem %d camadas; esperado %d"
                % (origem, sigma.size, n_expected)
            )
        return sigma
    return sigma_centers_uniform(n_expected)


# ---------------------------------------------------------------------------
# 6. Presets
# ---------------------------------------------------------------------------
def preset(name: Optional[str] = None) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """Devolve (regioes, campos) do preset name (default: etsa)."""
    key = DEFAULT_PRESET if name is None else str(name)
    if key not in PRESETS:
        raise KeyError(
            "preset MB09 desconhecido: %r (disponiveis: %s)" % (key, sorted(PRESETS))
        )
    regioes, campos = PRESETS[key]
    return tuple(regioes), tuple(campos)


def regions_of(name: Optional[str] = None) -> Tuple[str, ...]:
    """Somente as regioes (na ordem de combinacao) do preset name."""
    return preset(name)[0]


def fields_of(name: Optional[str] = None) -> Tuple[str, ...]:
    """Somente os campos perturbados do preset name."""
    return preset(name)[1]


# ---------------------------------------------------------------------------
# 7. Mascaras geograficas
# ---------------------------------------------------------------------------
def _region_tuple(
    region: Union[str, Sequence[float]]
) -> Tuple[float, float, float, float]:
    """Normaliza str -> tupla (lat_min, lat_max, lon_min, lon_max)."""
    if isinstance(region, str):
        key = region.upper()
        if key not in REGIONS:
            raise KeyError(
                "regiao MB09 desconhecida: %r (disponiveis: %s)"
                % (region, sorted(REGIONS))
            )
        return REGIONS[key]
    values = tuple(float(v) for v in region)
    if len(values) != 4:
        raise ValueError(
            "regiao deve ser nome ou (lat_min, lat_max, lon_min, lon_max)"
        )
    return (values[0], values[1], values[2], values[3])


def region_mask(lat, lon, region: Union[str, Sequence[float]]) -> np.ndarray:
    """Mascara booleana da caixa region para lat/lon (graus).

    Convencao de borda (D11, ver docs/mb09.md secao 5):

    * emenda interna meio-aberta: a linha de grade exatamente na fronteira
      pertence a caixa cujo MAXIMO e aquele valor. Assim lat=20 pertence a TR
      (e nao a NH), lat=-20 pertence a SH (e nao a TR) e lon=350 pertence a NSA
      (e nao sobra em nenhuma vizinha). Sem duplicar a linha entre NH e TR --
      importante porque ETSA soma as caixas sem renormalizar;
    * limite do dominio fechado: lat=-90 entra em SH, lat=+90 entra em NH (o
      polo nao pode ficar sem perturbacao, mesmo com a borda meio-aberta) e
      lon=0 pertence as caixas globais (o meridiano de Greenwich nao pode
      ficar descoberto).

    region e um nome de REGIONS ou a tupla
    (lat_min, lat_max, lon_min_360, lon_max_360).
    lon pode vir em qualquer convencao: e reduzido a [0, 360).
    Se lon_min > lon_max a caixa cruza 0/360 e a mascara e a uniao dos dois
    intervalos (nenhuma caixa do MB09 precisa disso, mas o suporte evita erro
    silencioso em caixas futuras).
    """
    lat_min, lat_max, lon_min, lon_max = _region_tuple(region)

    lat_a = np.asarray(lat, dtype=np.float64)
    lon_a = np.asarray(lon, dtype=np.float64)
    lat_b, lon_b = np.broadcast_arrays(lat_a, lon_a)
    lon_b = np.mod(lon_b, 360.0)

    # lat: (lat_min, lat_max], fechado no polo sul (-90)
    if lat_min <= -90.0:
        acima = lat_b >= lat_min
    else:
        acima = lat_b > lat_min
    inside_lat = acima & (lat_b <= lat_max)

    # lon: (lon_min, lon_max], fechado no meridiano 0 (inicio do dominio 0-360)
    if lon_min > lon_max:                       # caixa que cruza 0/360
        inside_lon = (lon_b > lon_min) | (lon_b <= lon_max)
    else:
        if lon_min <= 0.0:
            apos = lon_b >= lon_min
        else:
            apos = lon_b > lon_min
        inside_lon = apos & (lon_b <= lon_max)
    return inside_lat & inside_lon


def region_mask_preset(lat, lon, name: Optional[str] = None) -> Dict[str, np.ndarray]:
    """Mascaras de todas as caixas do preset name (dict nome -> mascara)."""
    return {r: region_mask(lat, lon, r) for r in regions_of(name)}


# ---------------------------------------------------------------------------
# 8. Regioes por indice de grade (grade da implementacao operacional 384x192)
# ---------------------------------------------------------------------------
#: Grade de referencia dos indices operacionais: (nlat, nlon) = (192, 384).
INDEX_GRID: Tuple[int, int] = (192, 384)

#: Caixas operacionais em indices 0-based, meio-aberto (j0, j1, i0, i1).
#:
#: Faixas originais (1-based, inclusivas; j = 1 e a linha mais ao NORTE)::
#:
#:     HN   J   1..75    I   1..384     -> (0,   75,   0, 384)
#:     TR   J  76..117   I   1..384     -> (75,  117,  0, 384)
#:     HS   J 118..192   I   1..384     -> (117, 192,  0, 384)
#:     NAS  J  75..108   I 276..372     -> (74,  108, 275, 372)
#:     SAS  J 109..161   I 266..362     -> (108, 161, 265, 362)
#:
#: Sao as caixas da implementacao operacional (oensMB09_bam_xc50) e valem
#: apenas para a grade 384x192: os indices NAO sao genericos.
INDEX_REGIONS: Dict[str, Tuple[int, int, int, int]] = {
    "HN": (0, 75, 0, 384),
    "TR": (75, 117, 0, 384),
    "HS": (117, 192, 0, 384),
    "NAS": (74, 108, 275, 372),
    "SAS": (108, 161, 265, 362),
}

#: Nome do artigo -> letra usada pelo operacional.
INDEX_REGION_ALIASES: Dict[str, str] = {
    "NH": "HN",
    "SH": "HS",
    "TR": "TR",
    "NSA": "NAS",
    "SSA": "SAS",
}

#: Ordem canonica das caixas por indice (norte -> sul, como no operacional).
INDEX_REGIONS_ORDER: Tuple[str, ...] = ("HN", "TR", "HS", "NAS", "SAS")


def _check_index_grid(nlat: int, nlon: int) -> Tuple[int, int]:
    """Valida a grade e devolve (nlat, nlon) como inteiros."""
    nlat_i, nlon_i = int(nlat), int(nlon)
    if (nlat_i, nlon_i) != INDEX_GRID:
        raise ValueError(
            "indices do operacional definidos apenas para a grade %dx%d "
            "(nlat, nlon); recebido %dx%d -- use region_mask(lat, lon, ...) "
            "para outras grades" % (INDEX_GRID[0], INDEX_GRID[1], nlat_i, nlon_i)
        )
    return nlat_i, nlon_i


def region_indices(nlat: int, nlon: int, region: str) -> Tuple[int, int, int, int]:
    """Caixa da regiao em indices de grade: (j0, j1, i0, i1), 0-based, meio-aberto.

    j cresce de norte para sul (j=1 e a linha mais ao norte) e i de oeste para
    leste, como na implementacao operacional (grade 384x192). Aceita a letra do
    operacional (HN, TR, HS, NAS, SAS) ou o nome do artigo (NH, SH, TR, NSA,
    SSA, maiusculo ou minusculo). Levanta ValueError se a grade nao for 384x192,
    porque os indices sao especificos dela.
    """
    _check_index_grid(nlat, nlon)
    key = str(region).strip().upper()
    key = INDEX_REGION_ALIASES.get(key, key)
    if key not in INDEX_REGIONS:
        raise KeyError(
            "regiao %r desconhecida; use uma de %s (ou os nomes do artigo %s)"
            % (region, list(INDEX_REGIONS), list(INDEX_REGION_ALIASES))
        )
    return INDEX_REGIONS[key]


def region_mask_by_index(nlat: int, nlon: int, region: str) -> np.ndarray:
    """Mascara booleana (nlat, nlon) da caixa region em indices de grade.

    Equivalente a region_mask para a grade 384x192, mas usando exatamente as
    faixas de indices do operacional (j=1 = norte). Ver region_indices.
    """
    nlat_i, nlon_i = _check_index_grid(nlat, nlon)
    j0, j1, i0, i1 = region_indices(nlat_i, nlon_i, region)
    mask = np.zeros((nlat_i, nlon_i), dtype=bool)
    mask[j0:j1, i0:i1] = True
    return mask


__all__ = [
    "REGIONS",
    "REGIONS_ORDER",
    "REGION_LABELS",
    "FIELDS",
    "PRESETS",
    "DEFAULT_PRESET",
    "AMPLITUDES_RANDOM",
    "AMPLITUDES_TARGET",
    "AMPLITUDES",
    "AMPLITUDE_UNITS",
    "N_PERTURBATIONS",
    "N_MEMBERS",
    "DISCARD_HOURS",
    "FORECAST_HOURS",
    "OUTPUT_INTERVAL_HOURS",
    "FORECAST_TIMES_HOURS",
    "M_TIMES",
    "Q_STDDEV_L28_X1E3",
    "Q_STDDEV_L28",
    "Q_STDDEV_OPER_L42_X1E3",
    "Q_STDDEV_OPER_L42",
    "Q_PROFILES",
    "Q_PROFILE_DEFAULT",
    "INDEX_GRID",
    "INDEX_REGIONS",
    "INDEX_REGIONS_ORDER",
    "INDEX_REGION_ALIASES",
    "sigma_centers_uniform",
    "sigma_centers_from_delta",
    "q_stddev_for_layers",
    "preset",
    "regions_of",
    "fields_of",
    "region_mask",
    "region_mask_preset",
    "region_indices",
    "region_mask_by_index",
]
