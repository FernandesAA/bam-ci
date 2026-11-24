"""Metricas de verificacao de ensemble do MB09 (Mendonca & Bonatti 2009).

Funcoes puras (numpy) para verificar previsoes por conjuntos, mais uma CLI que
le .ctl/.unf e escreve as metricas por lead/regiao.

Metricas (paginas do artigo em docs/mb09.md):

* pac        -- pattern anomaly correlation (p. 1444);
* rmse       -- raiz do erro quadratico medio do membro medio (p. 1444);
* spread     -- desvio-padrao dos membros em relacao a media do ensemble
                (p. 1444);
* brier_skill_score -- BSS e a decomposicao REL/RES/UNC (p. 1444);
* peca       -- perturbation versus error correlation analysis, modos "sin"
                (correlacao perturbacao-erro membro a membro, por default em
                modulo) e "opt" (combinacao otima + ganho de RMSE), Wei & Toth
                (2003) apud MB09 (p. 1444).

Convencoes (escolhas de implementacao, o artigo nao detalha):

* "members" tem forma (n_membros, ...); o membro medio e a media sobre o eixo 0
  e o espalhamento e o desvio-padrao sobre o eixo 0 com ddof=0 (desvio-padrao
  da populacao, como "std dev of ensemble members with respect to the ensemble
  mean");
* pac e a correlacao NAO centrada das anomalias (a anomalia ja e centrada pela
  climatologia), como no uso operacional;
* peca usa o erro do membro medio (ensemble_mean - truth) e correlacao de
  Pearson ponto a ponto; "sin" = media do MODULO das correlacoes dos n membros
  (signed=True usa o sinal; pares +p/-p cancelam a media com sinal); "opt" =
  regressao linear (sem intercepto) do erro sobre as perturbacoes
  x_m - ensemble_mean, com ganho de RMSE;
* brier_skill_score agrupa por caixas de probabilidade (default: 11 caixas de
  0.1, [0, 0.1), [0.1, 0.2), ..., [0.9, 1.0]) e devolve bss = (res - rel)/unc.

Uso como biblioteca::

    import bamci.score as sc
    m = sc.ensemble_mean(members)
    print(sc.rmse(m, truth), sc.spread(members).mean())
    print(sc.pac(m, truth, clim))
    print(sc.brier_skill_score(prob, obs_bin))
    print(sc.peca(members, truth))

Uso como CLI::

    python -m bamci.score --ens m1.ctl m2.ctl --truth anl.ctl \
        --var tvir --level 3 --region nh --leads 1,2,3 --out metrics.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from . import mb09_tables as mb

#: Regioes aceitas pela CLI ('global' = dominio inteiro).
CLI_REGIONS: Tuple[str, ...] = ("global", "nh", "sh", "tr", "nsa", "ssa")

#: Caixas de probabilidade default do BSS (11 caixas de 0.1, extremos fechados).
DEFAULT_PROB_BINS: Tuple[float, ...] = tuple(float(v) for v in np.linspace(0.0, 1.0, 11))

__all__ = [
    "CLI_REGIONS",
    "DEFAULT_PROB_BINS",
    "anomaly",
    "ensemble_mean",
    "spread",
    "rmse",
    "pac",
    "brier_skill_score",
    "ensemble_probability",
    "peca",
    "select_region",
    "verify_ensemble",
    "GradsField",
    "main",
]


# ---------------------------------------------------------------------------
# Funcoes puras
# ---------------------------------------------------------------------------
def anomaly(field, climatology) -> np.ndarray:
    """Anomalia = campo - climatologia (a climatologia pode ser escalar)."""
    return np.asarray(field, dtype=np.float64) - np.asarray(climatology, dtype=np.float64)


def ensemble_mean(members) -> np.ndarray:
    """Media do ensemble sobre o eixo 0 (membros)."""
    arr = np.asarray(members, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[0] < 1:
        raise ValueError("members deve ter ao menos 1 membro no eixo 0")
    return arr.mean(axis=0)


def spread(members, ddof: int = 0) -> np.ndarray:
    """Desvio-padrao dos membros em relacao ao membro medio (eixo 0).

    ddof=0 (default) e o desvio-padrao da populacao, definicao usual de
    espalhamento de ensemble; use ddof=1 para o estimador amostral.
    """
    arr = np.asarray(members, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[0] < 2:
        raise ValueError("spread precisa de ao menos 2 membros")
    return arr.std(axis=0, ddof=int(ddof))


def rmse(fcst, obs, axis=None) -> float:
    """Raiz do erro quadratico medio entre fcst e obs."""
    diff = np.asarray(fcst, dtype=np.float64) - np.asarray(obs, dtype=np.float64)
    if diff.size == 0:
        raise ValueError("rmse de campo vazio")
    return float(np.sqrt(np.mean(np.square(diff), axis=axis)))


def _pattern_corr(a: np.ndarray, b: np.ndarray, centered: bool = False) -> float:
    """Correlacao de dois campos achatados (nan se algum for constante)."""
    x = np.asarray(a, dtype=np.float64).ravel()
    y = np.asarray(b, dtype=np.float64).ravel()
    if x.size != y.size:
        raise ValueError("campos de tamanhos diferentes: %d e %d" % (x.size, y.size))
    if x.size < 2:
        return float("nan")
    if centered:
        x = x - x.mean()
        y = y - y.mean()
    den = float(np.sqrt(np.dot(x, x) * np.dot(y, y)))
    if den == 0.0:
        return float("nan")
    return float(np.dot(x, y) / den)


def pac(fcst, obs, clim) -> float:
    """Pattern anomaly correlation (nao centrada) das anomalias.

    r = sum(af*ao) / sqrt(sum(af^2) * sum(ao^2)), com af = fcst - clim e
    ao = obs - clim, calculada sobre todos os pontos do array (aplique a
    mascara da regiao antes, com select_region). Retorna nan quando um dos
    campos de anomalia e identicamente nulo.
    """
    af = anomaly(fcst, clim)
    ao = anomaly(obs, clim)
    return _pattern_corr(af, ao, centered=False)


def brier_skill_score(
    prob,
    obs_binary,
    thresholds: Optional[Sequence[float]] = None,
) -> Dict[str, object]:
    """BSS com decomposicao REL/RES/UNC agrupando por caixa de probabilidade.

    Parametros
    ----------
    prob:
        probabilidade prevista do evento, em [0, 1] (qualquer forma).
    obs_binary:
        0/1 observado do mesmo evento, mesma forma de prob.
    thresholds:
        bordas das caixas de probabilidade (crescentes, de 0 a 1). Default
        DEFAULT_PROB_BINS = [0, 0.1, ..., 1.0], ou seja 11 caixas
        [0,0.1), [0.1,0.2), ..., [0.9,1.0] (primeira e ultima fechadas).

    Decomposicao (Wilks 1995, apud MB09 p. 1444)::

        UNC = obar*(1 - obar)
        REL = sum_k n_k/N * (p_k - o_k)^2
        RES = sum_k n_k/N * (o_k - obar)^2
        BSS = (RES - REL) / UNC

    com p_k/o_k as medias de probabilidade e de observacao na caixa k.
    Previsao perfeita -> REL = 0 e BSS = 1; probabilidade constante igual a
    frequencia-base (previsao climatologica) -> BSS = 0.

    Retorna dict com bss, rel, res, unc, n (numero de pontos), base_rate e bins
    (as bordas usadas). bss e nan quando unc == 0 (evento sempre/nunca
    observado).
    """
    p = np.asarray(prob, dtype=np.float64).ravel()
    o = np.asarray(obs_binary, dtype=np.float64).ravel()
    if p.size != o.size:
        raise ValueError("prob e obs_binary com tamanhos diferentes: %d e %d"
                         % (p.size, o.size))
    if p.size == 0:
        raise ValueError("brier_skill_score de campo vazio")
    if np.any(~np.isfinite(p)) or np.any(~np.isfinite(o)):
        raise ValueError("prob/obs_binary tem valores nao finitos")
    if np.any(p < 0.0) or np.any(p > 1.0):
        raise ValueError("prob fora de [0, 1]")
    if np.any((o != 0.0) & (o != 1.0)):
        raise ValueError("obs_binary deve conter apenas 0 e 1")

    edges = np.asarray(DEFAULT_PROB_BINS if thresholds is None else thresholds,
                       dtype=np.float64)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("thresholds deve ser uma sequencia com ao menos 2 bordas")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("thresholds deve ser estritamente crescente")
    if edges[0] < 0.0 or edges[-1] > 1.0:
        raise ValueError("thresholds deve ficar em [0, 1]")

    n = int(p.size)
    obar = float(o.mean())
    unc = obar * (1.0 - obar)

    # caixa de cada probabilidade; a ultima caixa e fechada em 1.0
    idx = np.searchsorted(edges, p, side="right") - 1
    idx = np.clip(idx, 0, edges.size - 2)

    rel = 0.0
    res = 0.0
    for k in range(edges.size - 1):
        sel = idx == k
        nk = int(sel.sum())
        if nk == 0:
            continue
        pk = float(p[sel].mean())
        ok = float(o[sel].mean())
        rel += (nk / n) * (pk - ok) ** 2
        res += (nk / n) * (ok - obar) ** 2

    bss = (res - rel) / unc if unc > 0.0 else float("nan")
    return {
        "bss": float(bss),
        "rel": float(rel),
        "res": float(res),
        "unc": float(unc),
        "n": n,
        "base_rate": obar,
        "bins": [float(v) for v in edges],
    }


def ensemble_probability(members, threshold: float, above: bool = True) -> np.ndarray:
    """Fracao de membros com valor acima (ou abaixo) de threshold."""
    arr = np.asarray(members, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[0] < 1:
        raise ValueError("members deve ter ao menos 1 membro no eixo 0")
    return (arr > threshold).mean(axis=0) if above else (arr < threshold).mean(axis=0)


def peca(members, truth, mode: str = "both", signed: bool = False) -> Dict[str, object]:
    """Perturbation versus error correlation analysis (Wei & Toth 2003).

    members tem forma (n, ...) e truth a forma ... (mesmos pontos). Erro de
    referencia: e = membro_medio - truth; perturbacoes: p_m = x_m - membro_medio.

    * sin: media das n correlacoes de Pearson corr(p_m, e) -- a "average over
      14 individual PECA values" (p. 1444). Por default usa o MODULO de cada
      correlacao (signed=False), porque um ensemble MB09 e feito de pares
      +p/-p e as correlacoes com sinal se cancelariam (media ~0 por simetria);
      use signed=True para a media com sinal. sin_signed/sin_per_member ficam
      sempre disponiveis;
    * opt: regressao linear sem intercepto e ~ sum_m a_m p_m; devolve peca
      (correlacao entre a combinacao otima e o erro), coef, rmse_before,
      rmse_after (residuo e - q) e os ganhos absoluto/relativo de RMSE.

    mode: "sin", "opt" ou "both" (default); no modo escolhido apenas as
    metricas correspondentes sao calculadas (a outra volta None).
    """
    key = str(mode).lower()
    if key not in ("sin", "opt", "both"):
        raise ValueError("mode deve ser 'sin', 'opt' ou 'both' (recebido %r)" % (mode,))

    arr = np.asarray(members, dtype=np.float64)
    obs = np.asarray(truth, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[0] < 2:
        raise ValueError("peca precisa de ao menos 2 membros")
    if arr.shape[1:] != obs.shape:
        raise ValueError("truth com forma %s; esperado %s"
                         % (obs.shape, arr.shape[1:]))

    mean = arr.mean(axis=0)
    perms = arr - mean
    err = mean - obs

    n_members = int(arr.shape[0])
    n_points = int(err.size)
    p_flat = perms.reshape(n_members, -1)
    e_flat = err.ravel()

    out: Dict[str, object] = {
        "mode": key,
        "n_members": n_members,
        "n_points": n_points,
        "sin": None,
        "sin_per_member": None,
        "opt": None,
    }

    if key in ("sin", "both"):
        # membro praticamente igual ao membro medio: correlacao indefinida
        # (norma do campo perturbado desprezivel frente ao erro)
        e_norm = float(np.linalg.norm(e_flat))
        eps = 1.0e-12 * (e_norm if e_norm > 0.0 else 1.0)
        corrs = []
        for m in range(n_members):
            if np.linalg.norm(p_flat[m]) <= eps:
                corrs.append(float("nan"))
            else:
                corrs.append(_pattern_corr(p_flat[m], e_flat, centered=True))
        out["signed"] = bool(signed)
        out["sin_per_member"] = [float(c) for c in corrs]
        finite = [c for c in corrs if np.isfinite(c)]
        if finite:
            signed_mean = float(np.mean(finite))
            abs_mean = float(np.mean([abs(c) for c in finite]))
            out["sin_signed"] = signed_mean
            out["sin"] = signed_mean if signed else abs_mean
        else:
            out["sin_signed"] = float("nan")
            out["sin"] = float("nan")

    if key in ("opt", "both"):
        a_mat = p_flat.T                       # (n_pontos, n_membros)
        coef, *_ = np.linalg.lstsq(a_mat, e_flat, rcond=None)
        q = a_mat @ coef
        rmse_before = float(np.sqrt(np.mean(e_flat ** 2)))
        resid = e_flat - q
        rmse_after = float(np.sqrt(np.mean(resid ** 2)))
        out["opt"] = {
            "peca": _pattern_corr(q, e_flat, centered=True),
            "coef": [float(c) for c in coef],
            "rmse_before": rmse_before,
            "rmse_after": rmse_after,
            "gain": rmse_before - rmse_after,
            "gain_rel": (1.0 - rmse_after / rmse_before) if rmse_before > 0 else float("nan"),
        }

    return out


def select_region(values, lat, lon, region: Optional[str] = "global") -> np.ndarray:
    """Recorta values (..., nlat, nlon) na regiao usando mb.region_mask.

    region aceita os nomes da CLI (global, nh, sh, tr, nsa, ssa), qualquer nome
    de mb.REGIONS ou a tupla (lat_min, lat_max, lon_min, lon_max). O resultado
    tem a ultima dimensao achatada: (..., n_pontos).
    """
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("values deve ter ao menos 2 dimensoes (lat, lon)")
    name = "global" if region is None else region
    name = name.strip() if isinstance(name, str) else name
    if isinstance(name, str) and name.lower() == "global":
        mask = np.ones((arr.shape[-2], arr.shape[-1]), dtype=bool)
    else:
        key = name.upper() if isinstance(name, str) else name
        if isinstance(key, str) and key not in mb.REGIONS:
            raise KeyError("regiao %r desconhecida (use %s ou %s)"
                           % (region, list(CLI_REGIONS), sorted(mb.REGIONS)))
        lat = np.asarray(lat, dtype=np.float64)
        lon = np.asarray(lon, dtype=np.float64)
        if lat.size != arr.shape[-2] or lon.size != arr.shape[-1]:
            raise ValueError("lat/lon com %d/%d pontos; values tem %d/%d"
                             % (lat.size, lon.size, arr.shape[-2], arr.shape[-1]))
        mask = mb.region_mask(lat[:, None], lon[None, :], key)
    if not mask.any():
        raise ValueError("regiao %r nao tem nenhum ponto na grade" % (region,))
    return arr[..., mask]


# ---------------------------------------------------------------------------
# Verificacao por lead/regiao (usada pela CLI)
# ---------------------------------------------------------------------------
def verify_ensemble(
    members,
    truth,
    *,
    lat=None,
    lon=None,
    region: str = "global",
    clim=None,
    climb=None,
    event_k: float = 1.0,
    above: bool = True,
    bins: Optional[Sequence[float]] = None,
) -> Dict[str, object]:
    """Pacote de metricas MB09 para um lead.

    members (n, nlat, nlon); truth (nlat, nlon); clim/climb climatologia do
    campo e do desvio-padrao climatologico (mesma forma). Sem clim, a
    "climatologia" passa a ser a media espacial da verdade (e o desvio-padrao
    espacial das anomalias serve de escala do evento), o que e apenas uma
    conveniencia para rodar sem arquivo de climatologia -- registre isso nos
    metadados do JSON.

    O evento do BSS e anomalia > event_k * sigma_clim (ou < quando above=False),
    no espirito da p. 1444.
    """
    arr = np.asarray(members, dtype=np.float64)
    obs = np.asarray(truth, dtype=np.float64)
    if arr.shape[1:] != obs.shape:
        raise ValueError("truth com forma %s; esperado %s" % (obs.shape, arr.shape[1:]))

    if lat is None or lon is None:
        if region and str(region).lower() != "global":
            raise ValueError("region=%r exige lat/lon" % (region,))
        m = arr
        o = obs
        clim_m = None if clim is None else np.asarray(clim, dtype=np.float64)
        climb_m = None if climb is None else np.asarray(climb, dtype=np.float64)
    else:
        m = select_region(arr, lat, lon, region)
        o = select_region(obs, lat, lon, region)
        clim_m = None if clim is None else select_region(clim, lat, lon, region)
        climb_m = None if climb is None else select_region(climb, lat, lon, region)

    mean = m.mean(axis=0)
    spread_field = m.std(axis=0, ddof=0)

    if clim_m is None:
        clim_field = np.full(o.shape, float(o.mean()))
        clim_source = "media_espacial_da_verdade"
    else:
        clim_field = clim_m
        clim_source = "arquivo"

    if climb_m is None:
        sigma = float(np.std(anomaly(o, clim_field)))
    else:
        sigma = float(np.mean(climb_m))
    event_threshold = float(event_k) * sigma

    if above:
        prob = ensemble_probability(m, event_threshold, above=True)
        obs_event = (anomaly(o, clim_field) > event_threshold).astype(float)
    else:
        prob = ensemble_probability(m, -event_threshold, above=False)
        obs_event = (anomaly(o, clim_field) < -event_threshold).astype(float)
    bss = brier_skill_score(prob, obs_event, bins)

    p = peca(m, o, mode="both")

    return {
        "n_members": int(m.shape[0]),
        "n_points": int(o.size),
        "rmse": rmse(mean, o),
        "spread": float(np.mean(spread_field)),
        "pac": pac(mean, o, clim_field),
        "bss": bss,
        "peca_sin": p["sin"],
        "peca_sin_signed": p.get("sin_signed"),
        "peca_opt": None if p["opt"] is None else p["opt"]["peca"],
        "rmse_opt": None if p["opt"] is None else p["opt"]["rmse_after"],
        "event": {
            "k": float(event_k),
            "above": bool(above),
            "threshold": event_threshold,
            "sigma": sigma,
        },
        "climatology_source": clim_source,
    }


# ---------------------------------------------------------------------------
# Leitura de .ctl/.unf (reusa bamci.mb09._GradsSeries)
# ---------------------------------------------------------------------------
class GradsField:
    """Le UM campo de um .ctl/.unf reusando o leitor do MB09.

    Nao reimplementa leitura de GrADS: embrulha bamci.mb09._GradsSeries (layout
    do PostGrib, yrev, big/little endian). Como o construtor do MB09 resolve os
    5 campos canonicos (P,T,Q,U,V), o .ctl precisa te-los na vlist -- e o caso
    do PostGrib de producao; --var escolhe qual campo (canonico ou nao) e
    --level o nivel z (1-based) a ler.
    """

    def __init__(self, ctl_path: Union[str, Path], var: str, level: int = 1,
                 ps_units: str = "hpa"):
        from . import mb09  # import local: mantem as funcoes puras leves

        self.path = Path(ctl_path)
        self.series = mb09._GradsSeries(ctl_path, ps_units=ps_units)
        self.level = int(level)
        names = {v.name.lower(): v.name for v in self.series.ctl.variables}
        key = str(var).strip()
        if key.lower() in names:
            self.name = names[key.lower()]
        else:
            try:
                canonical = mb09._normalize_var_key(key)
            except KeyError:
                raise KeyError("variavel %r ausente no .ctl; disponiveis: %s"
                               % (var, sorted(names.values()))) from None
            self.name = self.series.vars[canonical]
        self.is_ps = self.series.vars.get("P") == self.name
        self.ps_units = str(ps_units).lower()
        self.lat = np.array(self.series.lat, dtype=np.float64)
        self.lon = np.array(self.series.lon, dtype=np.float64)

    @property
    def nt(self) -> int:
        return int(self.series.nt)

    def field(self, t: int) -> np.ndarray:
        """Campo 2-D (nlat, nlon) no tempo t (1-based)."""
        out = np.array(self.series._read(self.name, int(t), self.level),
                       dtype=np.float64)
        if self.is_ps and self.ps_units == "pa":
            out = out / 100.0
        return out


def _parse_leads(text: Optional[str]) -> Optional[List[int]]:
    if text is None:
        return None
    leads = []
    for chunk in str(text).replace(";", ",").split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        leads.append(int(chunk))
    if not leads:
        raise ValueError("--leads vazio")
    if len(set(leads)) != len(leads):
        raise ValueError("--leads com indices repetidos")
    return leads


def build_parser() -> argparse.ArgumentParser:
    """Parser da CLI de metricas."""
    p = argparse.ArgumentParser(
        prog="python -m bamci.score",
        description="Metricas de verificacao de ensemble do MB09 (PAC, RMSE, "
                    "spread, BSS/REL/RES, PECA) lendo .ctl/.unf.",
    )
    p.add_argument("--ens", nargs="+", required=True, metavar="CTL",
                   help="um .ctl por membro do ensemble")
    p.add_argument("--truth", required=True, metavar="CTL",
                   help=".ctl da verdade (analise/observacao)")
    p.add_argument("--var", required=True,
                   help="variavel a verificar (nome na vlist ou P/T/Q/U/V)")
    p.add_argument("--region", default="global", choices=list(CLI_REGIONS),
                   help="regiao de verificacao (default: global)")
    p.add_argument("--leads", default=None,
                   help="indices de tempo 1-based separados por virgula "
                        "(default: todos os tempos comuns)")
    p.add_argument("--level", type=int, default=1,
                   help="nivel z (1-based) lido do .ctl (default: 1)")
    p.add_argument("--clim", default=None, metavar="CTL",
                   help=".ctl da climatologia do campo (opcional)")
    p.add_argument("--clim-std", default=None, metavar="CTL",
                   help=".ctl do desvio-padrao climatologico (opcional)")
    p.add_argument("--event-k", type=float, default=1.0,
                   help="limiar do evento do BSS em multiplos do sigma "
                        "climatologico (default: 1.0, como no artigo)")
    p.add_argument("--below", action="store_true",
                   help="evento e anomalia < -k*sigma (default: > +k*sigma)")
    p.add_argument("--ps-units", default="hpa", choices=("hpa", "pa"),
                   help="unidade da pressao de superficie no .ctl (default: hpa)")
    p.add_argument("--out", default=None, metavar="JSON",
                   help="arquivo JSON de saida (default: so stdout)")
    return p


def _fmt(value, digits: int = 4) -> str:
    if value is None:
        return "-"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(v):
        return "nan"
    return ("%." + str(int(digits)) + "f") % v


def _table(rows: List[Dict[str, object]]) -> str:
    head = ("lead", "rmse", "spread", "pac", "bss", "rel", "res", "peca_sin", "peca_opt")
    widths = (5, 9, 9, 8, 9, 9, 9, 10, 10)
    line = " ".join(h.rjust(w) for h, w in zip(head, widths))
    out = [line, "-" * len(line)]
    for row in rows:
        out.append(" ".join(str(row[h]).rjust(w) for h, w in zip(head, widths)))
    return "\n".join(out)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entrada da CLI: devolve 0 em sucesso, 1 em erro de entrada."""
    args = build_parser().parse_args(argv)
    try:
        leads = _parse_leads(args.leads)
    except (TypeError, ValueError) as exc:
        print("erro: --leads invalido: %s" % exc, file=sys.stderr)
        return 1

    try:
        ens = [GradsField(path, args.var, args.level, args.ps_units)
               for path in args.ens]
        truth = GradsField(args.truth, args.var, args.level, args.ps_units)
        clim = (GradsField(args.clim, args.var, args.level, args.ps_units)
                if args.clim else None)
        clim_std = (GradsField(args.clim_std, args.var, args.level, args.ps_units)
                    if args.clim_std else None)
    except (KeyError, ValueError, OSError) as exc:
        print("erro: %s" % exc, file=sys.stderr)
        return 1

    ntimes = min([truth.nt] + [f.nt for f in ens] +
                 ([clim.nt] if clim else []) +
                 ([clim_std.nt] if clim_std else []))
    if leads is None:
        leads = list(range(1, ntimes + 1))
    bad = [t for t in leads if t < 1 or t > ntimes]
    if bad:
        print("erro: leads %s fora de 1..%d" % (bad, ntimes), file=sys.stderr)
        return 1

    for f in [truth] + ens:
        if not (np.allclose(f.lat, truth.lat) and np.allclose(f.lon, truth.lon)):
            print("erro: grades diferentes entre %s e %s" % (f.path, truth.path),
                  file=sys.stderr)
            return 1

    rows: List[Dict[str, object]] = []
    per_lead: Dict[str, object] = {}
    for t in leads:
        members = np.stack([f.field(t) for f in ens], axis=0)
        obs = truth.field(t)
        clim_field = clim.field(t) if clim else None
        climb_field = clim_std.field(t) if clim_std else None
        metrics = verify_ensemble(
            members, obs, lat=truth.lat, lon=truth.lon, region=args.region,
            clim=clim_field, climb=climb_field, event_k=args.event_k,
            above=not args.below,
        )
        per_lead[str(t)] = metrics
        rows.append({
            "lead": t,
            "rmse": _fmt(metrics["rmse"]),
            "spread": _fmt(metrics["spread"]),
            "pac": _fmt(metrics["pac"]),
            "bss": _fmt(metrics["bss"]["bss"]),
            "rel": _fmt(metrics["bss"]["rel"]),
            "res": _fmt(metrics["bss"]["res"]),
            "peca_sin": _fmt(metrics["peca_sin"]),
            "peca_opt": _fmt(metrics["peca_opt"]),
        })

    payload = {
        "method": "MB09 (Mendonca & Bonatti 2009, MWR 137:1438-1459)",
        "variable": args.var,
        "level": int(args.level),
        "region": args.region,
        "n_members": len(ens),
        "leads": list(leads),
        "ps_units": args.ps_units,
        "event": {"k": float(args.event_k), "above": not args.below},
        "climatology": ("arquivo" if (clim or clim_std) else "media_espacial_da_verdade"),
        "metrics": per_lead,
    }

    print("MB09 score | var=%s level=%d region=%s membros=%d leads=%s"
          % (args.var, args.level, args.region, len(ens),
             ",".join(str(t) for t in leads)))
    print(_table(rows))

    if args.out:
        out = Path(args.out)
        if out.parent and not out.parent.exists():
            out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                       encoding="utf-8")
        print("\njson: %s" % out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
