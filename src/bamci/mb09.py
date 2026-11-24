"""Nucleo do metodo de perturbacao por EOF do MB09 (Mendonca & Bonatti 2009).

Implementa o metodo descrito em 'docs/mb09.md' a partir das tabelas de
'bamci.mb09_tables': n=7 perturbacoes aleatorias -> janela util de 6 a 36 h
(M=11, passo 3 h) -> EOF por campo e por regiao (1o autovetor) -> reescala para
as amplitudes prescritas -> 2n+1 = 15 membros sobre a MESMA analise de controle.

API publica
-----------
State
    Estado em grade (t, q, u, v, ps, topo) com:
      * 'State.from_ic_file(path)'      -> State  (le um IC espectral do BAM)
      * 'state.to_ic_file(out, template=None)' -> Path (grava um IC espectral)
random_deltas(state, rng, sigma=None) -> dict
    Perturbacoes aleatorias gaussianas de media zero (vento 3.0 m/s, T 0.6 K,
    ps 1.0 hPa, q = Tabela 1 por camada). Unidades do Estado (ps em hPa).
perturbed_state(state, delta) -> State
    Soma o delta ao estado e clipa q em 0.
build_eofs(diffs, preset="etsa", weights=None, *, keep_basis=False) -> dict
    'diffs' = {membro: {campo: [State x M]}} (ou {membro: [State x M]}), com
    'campo' em {"P","T","Q","U","V"}; U/V sao combinados em vento complexo
    u+iv. Devolve os modos (1o autovetor) por membro/campo/regiao + autovalores
    + sigma_i, prontos para 'apply_modes'.
apply_modes(control, eofs, preset="etsa", sigma=None) -> list[State]
    2n estados na ordem 2k-1 = x_c + p_k (indice 2k-2) e 2k = x_c - p_k
    (indice 2k-1), k = 1..n.
write_ensemble(control_path, states, out_dir, prefix) -> list[Path]
    Grava o controle + os 2n estados (15 arquivos) e 'members.json'.
run_mb09(control_ic, member_ics, out_dir, preset, seed) -> list[Path]
    Build + apply em um comando, quando as series de previsao ja existem.

Convencoes e decisoes de implementacao (o artigo NAO especifica; ver D1..D10 em
docs/mb09.md):

* X tem forma (M, N) por campo/regiao, com N = n_pontos_da_regiao * n_niveis
  (niveis empilhados, D5). A media temporal e removida de X (D2). Sem pesos de
  area por padrao (D3; 'weights' fica como opcao explicita).
* truque dual: G = (1/M) X X^H e' (M, M) = 11x11; autovalores/autovetores de G
  e EOF = X^H u normalizado. Nunca se forma N x N.
* vento e' complexo u + i v (G hermitiana). Para o vento, 'sigma_i' e' o desvio
  padrao do campo complexo na regiao (equivale ao RMS do modulo depois de
  remover a media espacial) e o alvo 5.0 m/s e' interpretado nessa mesma norma
  (RMS do vento vetorial). O artigo nao define sigma_i para o vento.
* Q: reescala POR CAMADA (Tabela 1), como em docs/mb09.md; os demais campos usam
  um fator ESCALAR por campo/regiao.
* regioes sobrepostas (ETSA) somam de forma aditiva, sem renormalizar (D9).
* P: o Estado guarda ps em hPa; o IC guarda ln(ps/10) e a conversao e'
  d ln(ps/10) = d ps / ps (D6). T: perturbacao aplicada em Tv (D7), que e' o
  campo 't' do State. q e' clipado em 0 apos a soma (D8).
* regioes vem de 'bamci.mb09_tables.region_mask'; campos/regioes do preset
  default sao ETSA.

Vento no IC espectral
---------------------
O IC guarda vorticidade e divergencia (nao u/v). 'State.from_ic_file' reconstroi
u/v pelo potencial de velocidade e funcao de corrente (mesmo operador de
derivadas de 'bamci.derivatives'). Como o vento em grade nao e' uma funcao
escalar de banda limitada (o fator 1/cos(fi)), a ida u/v -> div/vort -> u/v NAO
e' exata: 'to_ic_file' grava, portanto, a DIFERENCA em relacao ao vento que ele
proprio reconstruiu do IC (divergencia e vorticidade sao lineares no vento).
Assim 'State -> IC -> State' e' exato para o estado nao perturbado, e uma
perturbacao de vento e' projetada em div/vort com o truncamento do IC.

Aviso de escopo: a integracao nao linear de 36 h e' externa; 'run_mb09' consome
as saidas ja integradas (IC-format ou GrADS/PostGrib) e faz build + apply.
'build_eofs_from_files' (CLI build) seleciona a janela 6..36 h passo 3 h
(M=11, default) via 'select_window'/'select_times'; tempos extras (0, 3 h ou
> 36 h) sao ignorados com um aviso, e horarios da janela ausentes no .ctl
geram erro listando os tempos disponiveis.
"""

from __future__ import annotations

import json
import math
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from . import mb09_tables as mb
from .derivatives import A_EARTH, dpln_columns, spectral_div_vort
from .ic_format import ICHeader, n_coeffs, read_ic, write_ic
from .levels import gaussian_latitudes, grid_shape
from .spectral import analysis, file_order, from_file_order, legendre_columns, synthesis

__all__ = [
    "State",
    "Mode",
    "random_deltas",
    "perturbed_state",
    "build_eofs",
    "build_eofs_from_files",
    "select_window",
    "select_times",
    "apply_modes",
    "write_ensemble",
    "run_mb09",
    "save_eofs",
    "load_eofs",
]

_FIELD_ATTR = {"P": "ps", "T": "t", "Q": "q", "U": "u", "V": "v"}
_EOF_FIELDS = ("P", "T", "Q", "wind")


# ---------------------------------------------------------------------------
# utilidades espectrais
# ---------------------------------------------------------------------------
def _mend_from_ncoeff(n):
    """Truncacao mend a partir do numero de coeficientes (M+1)(M+2)."""
    n = int(n)
    m = int(round((-3.0 + math.sqrt(1.0 + 4.0 * n)) / 2.0))
    if n_coeffs(m) != n:
        raise ValueError("ncoeff=%d nao e' (mend+1)(mend+2)" % (n,))
    return m


def _trig(m, lam, nlon):
    """(cos, sen, fator) da ordem zonal m, na convencao de spectral.synthesis."""
    if m == 0:
        return np.ones(nlon), np.zeros(nlon), 1.0
    return np.cos(m * lam), np.sin(m * lam), 2.0


def _synth_weighted(coef, mend, nlat, nlon, weight=None):
    """Sintese de um campo escalar a partir de coeficientes (m+1, n+1)."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m, pln_m in legendre_columns(mend, mu):
        c = coef[m, m:].copy()
        if weight is not None:
            c = c * weight[m, m:]
        y = pln_m.T @ c
        trig, sint, fac = _trig(m, lam, nlon)
        out += fac * (np.real(y)[:, None] * trig[None, :] + np.imag(y)[:, None] * sint[None, :])
    return out


def _synth_dmu(coef, mend, nlat, nlon):
    """Sintese de d(.)/d mu (derivada em latitude)."""
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    out = np.zeros((nlat, nlon))
    for m, dpln_m in dpln_columns(mend, mu):
        c = coef[m, m:]
        y = dpln_m.T @ c
        trig, sint, fac = _trig(m, lam, nlon)
        out += fac * (np.real(y)[:, None] * trig[None, :] + np.imag(y)[:, None] * sint[None, :])
    return out


def _wind_from_divvort(div_sh, vort_sh, mend, nlat, nlon):
    """u, v em grade a partir dos coeficientes espectrais de div/vort.

    Resolve o Laplaciano (psi = -a^2 zeta / n(n+1), chi = -a^2 delta / n(n+1)) e
    usa u = -(1/a) dpsi/dphi + (1/(a cos fi)) dchi/dlam,
        v =  (1/(a cos fi)) dpsi/dlam + (1/a) dchi/dphi.
    O modo (0,0) nao tem vento associado e fica zerado.
    """
    div_sh = np.asarray(div_sh, dtype=np.complex128)
    vort_sh = np.asarray(vort_sh, dtype=np.complex128)
    n = np.arange(mend + 1, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(n > 0, -A_EARTH ** 2 / (n * (n + 1.0)), 0.0)
    psi = vort_sh * fac[None, :]
    chi = div_sh * fac[None, :]

    lat = np.radians(gaussian_latitudes(nlat))
    coslat = np.cos(lat)[:, None]
    lam = 2.0 * np.pi * np.arange(nlon) / nlon
    wlam = np.zeros((mend + 1, mend + 1), dtype=np.complex128)
    for m in range(mend + 1):
        wlam[m, :] = 1j * m

    dpsi_dlam = _synth_weighted(psi, mend, nlat, nlon, wlam)
    dchi_dlam = _synth_weighted(chi, mend, nlat, nlon, wlam)
    dpsi_dmu = _synth_dmu(psi, mend, nlat, nlon)
    dchi_dmu = _synth_dmu(chi, mend, nlat, nlon)

    u = -(coslat / A_EARTH) * dpsi_dmu + dchi_dlam / (A_EARTH * coslat)
    v = dpsi_dlam / (A_EARTH * coslat) + (coslat / A_EARTH) * dchi_dmu
    return u, v


def _divvort_sh_from_wind(u, v, nlat, nlon, mend):
    """Coeficientes espectrais (div, vort) de um campo de vento em grade."""
    div, vort = spectral_div_vort(u, v, nlat, nlon, mend)
    return analysis(div, mend), analysis(vort, mend)


def _header_from_data(data: Mapping[str, object]) -> ICHeader:
    """Reconstroi ICHeader a partir do dict de bamci.ic_format.read_ic."""
    kmax = int(data["kmax"])
    a = np.asarray(data["a_hybr"], dtype=np.float64)
    b = np.asarray(data["b_hybr"], dtype=np.float64)
    base = dict(
        forecast_day=int(data["forecast_day"]),
        time_of_day=float(data["time_of_day"]),
        date_inp=tuple(int(v) for v in data["date_inp"]),
        date_cur=tuple(int(v) for v in data["date_cur"]),
    )
    if b.size == kmax + 1:
        return ICHeader(a_hybr=a, b_hybr=b, **base)
    if b.size == kmax:
        return ICHeader(sigma_inter=a, sigma_layers=b, **base)
    raise ValueError(
        "cabecalho com %d valores de b para kmax=%d (esperado %d ou %d)"
        % (b.size, kmax, kmax + 1, kmax)
    )


def _default_header(nlev: int) -> ICHeader:
    """Cabecalho sigma uniforme, usado quando um State nao tem template."""
    inter = np.linspace(1.0, 0.0, int(nlev) + 1)
    layers = 0.5 * (inter[:-1] + inter[1:])
    return ICHeader(sigma_inter=inter, sigma_layers=layers)


def _header_from_template(template, state: "State") -> ICHeader:
    if template is None:
        return state.header if state.header is not None else _default_header(state.nlev)
    if isinstance(template, ICHeader):
        return template
    return _header_from_data(read_ic(template))


# ---------------------------------------------------------------------------
# Estado em grade
# ---------------------------------------------------------------------------
@dataclass
class State:
    """Estado atmosferico em grade gaussiana (nivel, latitude, longitude).

    Campos de nivel, todos com forma (nlev, nlat, nlon):
      * 't'  temperatura virtual (K) -- e' o Tvir do IC (D7);
      * 'q'  umidade especifica (kg/kg);
      * 'u'  vento zonal (m/s);
      * 'v'  vento meridional (m/s).
    'ps' tem forma (nlat, nlon) e fica em hPa. 'topo' (nlat, nlon, m) e'
    opcional. 'lat' (nlat) e 'lon' (nlon) em graus; 'lon' na convencao 0..360.

    Campos privados guardam o que o IC espectral nao tem em grade:
      * '_mend', '_div_sh', '_vort_sh' (kmax, m+1, n+1) complexos;
      * '_u_ref', '_v_ref': vento reconstruido desses coeficientes. E' a base
        do delta de vento gravado por 'to_ic_file' (ver docstring do modulo).
    """

    t: np.ndarray
    q: np.ndarray
    u: np.ndarray
    v: np.ndarray
    ps: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    topo: Optional[np.ndarray] = None
    header: Optional[ICHeader] = None
    _mend: Optional[int] = None
    _div_sh: Optional[np.ndarray] = None
    _vort_sh: Optional[np.ndarray] = None
    _u_ref: Optional[np.ndarray] = None
    _v_ref: Optional[np.ndarray] = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.t = np.asarray(self.t, dtype=np.float64)
        self.q = np.asarray(self.q, dtype=np.float64)
        self.u = np.asarray(self.u, dtype=np.float64)
        self.v = np.asarray(self.v, dtype=np.float64)
        self.ps = np.asarray(self.ps, dtype=np.float64)
        self.lat = np.asarray(self.lat, dtype=np.float64).ravel()
        self.lon = np.asarray(self.lon, dtype=np.float64).ravel()
        if self.t.ndim != 3:
            raise ValueError("t deve ser (nlev, nlat, nlon); recebido %s" % (self.t.shape,))
        for name in ("q", "u", "v"):
            arr = getattr(self, name)
            if arr.shape != self.t.shape:
                raise ValueError("%s tem forma %s, esperada %s" % (name, arr.shape, self.t.shape))
        if self.ps.shape != self.t.shape[1:]:
            raise ValueError("ps tem forma %s, esperada %s" % (self.ps.shape, self.t.shape[1:]))
        if self.lat.size != self.t.shape[1] or self.lon.size != self.t.shape[2]:
            raise ValueError("lat/lon incompativeis com o grid %s" % (self.t.shape,))
        if self.topo is not None:
            self.topo = np.asarray(self.topo, dtype=np.float64)
            if self.topo.shape != self.ps.shape:
                raise ValueError("topo tem forma %s, esperada %s" % (self.topo.shape, self.ps.shape))

    # -- propriedades ------------------------------------------------------
    @property
    def nlev(self) -> int:
        return int(self.t.shape[0])

    @property
    def nlat(self) -> int:
        return int(self.t.shape[1])

    @property
    def nlon(self) -> int:
        return int(self.t.shape[2])

    @property
    def shape(self) -> Tuple[int, int, int]:
        return self.t.shape

    @property
    def tvir(self) -> np.ndarray:
        """Alias explicito de 't' (o IC guarda temperatura virtual, D7)."""
        return self.t

    @property
    def lnps(self) -> np.ndarray:
        """ln(ps/10 hPa), a variavel de pressao do IC."""
        return np.log(self.ps / 10.0)

    def copy(self) -> "State":
        return State(
            t=self.t.copy(), q=self.q.copy(), u=self.u.copy(), v=self.v.copy(),
            ps=self.ps.copy(), lat=self.lat.copy(), lon=self.lon.copy(),
            topo=None if self.topo is None else self.topo.copy(),
            header=self.header, _mend=self._mend,
            _div_sh=None if self._div_sh is None else self._div_sh.copy(),
            _vort_sh=None if self._vort_sh is None else self._vort_sh.copy(),
            _u_ref=None if self._u_ref is None else self._u_ref.copy(),
            _v_ref=None if self._v_ref is None else self._v_ref.copy(),
            meta=dict(self.meta),
        )

    # -- IC ---------------------------------------------------------------
    @classmethod
    def from_ic_file(cls, path: Union[str, Path]) -> "State":
        """Le um IC espectral do BAM e devolve o estado em grade."""
        data = read_ic(path)
        kmax = int(data["kmax"])
        ncoeff = int(data["ncoeff"])
        mend = _mend_from_ncoeff(ncoeff)
        nlon, nlat = grid_shape(mend)
        lat = gaussian_latitudes(nlat)
        lon = np.arange(nlon, dtype=np.float64) * 360.0 / nlon

        def grid(vec):
            return synthesis(from_file_order(vec, mend), mend, nlat, nlon)

        t = np.empty((kmax, nlat, nlon))
        q = np.empty((kmax, nlat, nlon))
        for k in range(kmax):
            t[k] = grid(np.asarray(data["tvir"])[:, k])
            q[k] = grid(np.asarray(data["sphu"])[:, k])
        div_sh = np.empty((kmax, mend + 1, mend + 1), dtype=np.complex128)
        vort_sh = np.empty_like(div_sh)
        u = np.empty((kmax, nlat, nlon))
        v = np.empty_like(u)
        for k in range(kmax):
            div_sh[k] = from_file_order(np.asarray(data["divg"])[:, k], mend)
            vort_sh[k] = from_file_order(np.asarray(data["vort"])[:, k], mend)
            u[k], v[k] = _wind_from_divvort(div_sh[k], vort_sh[k], mend, nlat, nlon)

        lnps = grid(np.asarray(data["lnps"]))
        ps = 10.0 * np.exp(lnps)
        topo = grid(np.asarray(data["topo"]))
        return cls(
            t=t, q=q, u=u, v=v, ps=ps, lat=lat, lon=lon, topo=topo,
            header=_header_from_data(data), _mend=mend,
            _div_sh=div_sh, _vort_sh=vort_sh, _u_ref=u.copy(), _v_ref=v.copy(),
        )

    def to_ic_file(self, out_path: Union[str, Path], template=None) -> Path:
        """Grava o estado como IC espectral do BAM.

        'template' pode ser um caminho de IC (copia o cabecalho/datas) ou um
        ICHeader. Sem template usa 'self.header' ou um cabecalho sigma uniforme.
        """
        if self._mend is None:
            raise ValueError("State sem truncacao espectral (_mend): use State.from_ic_file")
        mend = int(self._mend)
        kmax = self.nlev
        nlat, nlon = self.nlat, self.nlon
        ncoeff = n_coeffs(mend)
        header = _header_from_template(template, self)

        tvir = np.empty((ncoeff, kmax))
        sphu = np.empty((ncoeff, kmax))
        for k in range(kmax):
            tvir[:, k] = file_order(analysis(self.t[k], mend))
            sphu[:, k] = file_order(analysis(self.q[k], mend))

        lnps = file_order(analysis(self.lnps, mend))
        topo = (file_order(analysis(self.topo, mend)) if self.topo is not None
                else np.zeros(ncoeff))

        if self._div_sh is not None and self._u_ref is not None and self._v_ref is not None:
            div_sh = self._div_sh.copy()
            vort_sh = self._vort_sh.copy()
            du = self.u - self._u_ref
            dv = self.v - self._v_ref
            if np.any(du) or np.any(dv):
                for k in range(kmax):
                    dd, dz = _divvort_sh_from_wind(du[k], dv[k], nlat, nlon, mend)
                    div_sh[k] = div_sh[k] + dd
                    vort_sh[k] = vort_sh[k] + dz
        else:
            div_sh = np.empty((kmax, mend + 1, mend + 1), dtype=np.complex128)
            vort_sh = np.empty_like(div_sh)
            for k in range(kmax):
                div_sh[k], vort_sh[k] = _divvort_sh_from_wind(
                    self.u[k], self.v[k], nlat, nlon, mend)

        divg = np.stack([file_order(div_sh[k]) for k in range(kmax)], axis=1)
        vort = np.stack([file_order(vort_sh[k]) for k in range(kmax)], axis=1)
        return write_ic(out_path, header, topo, lnps, tvir, divg, vort, sphu)


# ---------------------------------------------------------------------------
# perturbacoes aleatorias e aplicacao
# ---------------------------------------------------------------------------
_SIGMA_ALIASES = {
    "u": "u", "v": "v", "wind": "wind", "vento": "wind",
    "t": "t", "tvir": "t", "temp": "t", "temperatura": "t",
    "p": "ps", "ps": "ps", "pressao": "ps",
    "q": "q", "sphu": "q", "umid": "q", "umidade": "q",
}


def _sigma_entry(key: str) -> str:
    k = str(key).lower()
    if k not in _SIGMA_ALIASES:
        raise KeyError("sigma desconhecido: %r (use u, v, wind, t, ps, q)" % (key,))
    return _SIGMA_ALIASES[k]


def random_deltas(state: State, rng, sigma: Optional[Mapping[str, object]] = None) -> dict:
    """Perturbacoes aleatorias gaussianas de media zero (passo 1 do MB09).

    Defaults (p. 1441-1442): vento 3.0 m/s (u e v independentes), T 0.6 K,
    ps 1.0 hPa e q com o perfil da Tabela 1 por camada sigma. 'sigma' permite
    sobrescrever por campo (chaves u, v, wind, t, ps, q); o valor de q pode ser
    escalar ou um vetor de nlev valores.

    Retorna dict com 'u','v','t','ps','q' nas unidades do State (ps em hPa).
    """
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    sig = {"u": 3.0, "v": 3.0, "t": 0.6, "ps": 1.0}
    qprof = np.asarray(mb.q_stddev_for_layers(state.nlev), dtype=np.float64)

    if sigma:
        for key, value in sigma.items():
            k = _sigma_entry(key)
            if k == "wind":
                sig["u"] = sig["v"] = float(value)
            elif k in ("u", "v", "t", "ps"):
                sig[k] = float(value)
            else:  # q
                arr = np.asarray(value, dtype=np.float64)
                if arr.ndim == 0:
                    qprof = np.full(state.nlev, float(arr))
                else:
                    qprof = np.broadcast_to(arr, (state.nlev,)).astype(np.float64)

    return {
        "u": rng.normal(0.0, sig["u"], state.u.shape),
        "v": rng.normal(0.0, sig["v"], state.v.shape),
        "t": rng.normal(0.0, sig["t"], state.t.shape),
        "ps": rng.normal(0.0, sig["ps"], state.ps.shape),
        "q": rng.normal(0.0, qprof[:, None, None], state.q.shape),
    }


_DELTA_ALIASES = {"tvir": "t", "p": "ps", "sphu": "q"}


def perturbed_state(state: State, delta: Mapping[str, np.ndarray]) -> State:
    """Aplica 'delta' ao estado (soma) e clipa q em 0 (D8).

    Chaves aceitas: u, v, t (ou tvir), ps (hPa), q. Delta e' a perturbacao,
    nao o campo absoluto.
    """
    norm: Dict[str, np.ndarray] = {}
    for key, value in delta.items():
        k = _DELTA_ALIASES.get(str(key), str(key))
        if k not in ("u", "v", "t", "ps", "q"):
            raise KeyError("delta desconhecido: %r (use u, v, t, ps, q)" % (key,))
        norm[k] = np.asarray(value, dtype=np.float64)

    out = state.copy()
    if "u" in norm:
        out.u = state.u + norm["u"]
    if "v" in norm:
        out.v = state.v + norm["v"]
    if "t" in norm:
        out.t = state.t + norm["t"]
    if "ps" in norm:
        out.ps = state.ps + norm["ps"]
    if "q" in norm:
        out.q = np.maximum(state.q + norm["q"], 0.0)
    return out


# ---------------------------------------------------------------------------
# EOF
# ---------------------------------------------------------------------------
@dataclass
class Mode:
    """1o modo EOF de uma serie (membro/campo/regiao) e metadados.

    'mode' tem forma (nlev, nlat, nlon); e' complexo para o vento (u+iv) e real
    para P/T/Q, com zeros fora da regiao. 'sigma_i' e' o desvio padrao do modo
    cru na regiao (escalar para P/T/vento; vetor por nivel para Q).
    'eigenvectors' so e' preenchido com keep_basis=True.
    """

    field: str
    region: str
    mode: np.ndarray
    eigenvalues: np.ndarray
    sigma_i: np.ndarray
    eigenvectors: Optional[np.ndarray] = None

    def scaled(self, target) -> np.ndarray:
        """Reescala o modo para 'target' (sigma_f) via sigma_f/sigma_i."""
        tgt = np.asarray(target, dtype=np.float64)
        if tgt.ndim == 1:
            tgt = tgt.reshape((-1, 1, 1))
        sig = np.asarray(self.sigma_i, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(sig > 0.0, tgt / sig, 0.0)
        return self.mode * ratio


def _fix_phase(mode: np.ndarray) -> np.ndarray:
    """Fixa a fase/sinal arbitrario do autovetor (ponto de maior modulo > 0)."""
    if mode.size == 0:
        return mode
    i = int(np.argmax(np.abs(mode)))
    v = mode[i]
    if np.iscomplexobj(mode):
        if v != 0:
            mode = mode * (np.conj(v) / abs(v))
    elif v < 0:
        mode = -mode
    return mode


def _eof_first(X: np.ndarray, weights=None, keep_basis: bool = False):
    """EOF dominante pelo truque dual: G = (1/M) X X^H (M x M).

    X: (M, N) real ou complexo; media temporal removida aqui. 'weights' e' um
    vetor (N,) de pesos metricos opcionais (aplica sqrt(w) em X e desfaz no
    modo). Devolve (modo normalizado, autovalores decrescentes, base ou None).
    """
    M = X.shape[0]
    Xc = X - X.mean(axis=0, keepdims=True)
    if weights is not None:
        w = np.asarray(weights, dtype=np.float64).ravel()
        if w.size != X.shape[1]:
            raise ValueError("weights tem %d valores; esperado %d" % (w.size, X.shape[1]))
        if np.any(w <= 0.0):
            raise ValueError("weights deve ser positivo em todos os pontos")
        Xw = Xc * np.sqrt(w)[None, :]
    else:
        w = None
        Xw = Xc

    G = (Xw @ Xw.conj().T) / float(M)
    G = 0.5 * (G + G.conj().T)
    evals, evecs = np.linalg.eigh(G)
    order = np.argsort(evals)[::-1]
    evals = np.clip(np.real(evals[order]), 0.0, None)
    evecs = evecs[:, order]

    mode = Xw.conj().T @ evecs[:, 0]
    if w is not None:
        mode = mode / np.sqrt(w)
    norm = float(np.linalg.norm(mode))
    if norm > 0.0:
        mode = mode / norm
    mode = _fix_phase(mode)
    return mode, evals, (evecs if keep_basis else None)


def _stack_states(states: Sequence[State], attr: str) -> np.ndarray:
    if not states:
        raise ValueError("serie vazia para o campo %r" % (attr,))
    arrays = [np.asarray(getattr(s, attr)) for s in states]
    shapes = {a.shape for a in arrays}
    if len(shapes) != 1:
        raise ValueError("series com formas diferentes: %s" % (sorted(shapes),))
    return np.stack(arrays)


def _member_arrays(member) -> Dict[str, np.ndarray]:
    """Normaliza um membro para {campo: array (M, nlev, nlat, nlon)}.

    Aceita {campo: [State x M]} (campos P/T/Q/U/V) ou [State x M]. Para U/V o
    vento vira complexo u+iv (um unico EOF por nivel).
    """
    if isinstance(member, Mapping):
        lists = {str(k).upper(): list(v) for k, v in member.items()}
        arrays: Dict[str, np.ndarray] = {}
        for f in ("P", "T", "Q"):
            if f in lists:
                arr = _stack_states(lists[f], _FIELD_ATTR[f])
                # P e' ps 2-D: guarda com um eixo de nivel (K=1) para que todos
                # os campos tenham a mesma forma (M, K, nlat, nlon).
                arrays[f] = arr[:, None, :, :] if f == "P" else arr
        if "U" in lists or "V" in lists:
            src = lists.get("U") or lists.get("V")
            u = _stack_states(lists.get("U", src), "u")
            v = _stack_states(lists.get("V", src), "v")
            arrays["wind"] = u + 1j * v
        return arrays

    states = list(member)
    return {
        "P": _stack_states(states, "ps")[:, None, :, :],
        "T": _stack_states(states, "t"),
        "Q": _stack_states(states, "q"),
        "wind": _stack_states(states, "u") + 1j * _stack_states(states, "v"),
    }


def _field_weights(weights, field: str):
    if weights is None:
        return None
    if isinstance(weights, Mapping):
        return weights.get(field, weights.get(field.lower(), weights.get("wind")))
    return weights


def _sigma_i_for(mode_grid: np.ndarray, field: str, mask: np.ndarray) -> np.ndarray:
    values = mode_grid[:, mask]
    if field == "Q":
        # desvio padrao por camada (Tabela 1 e' por nivel sigma)
        std = np.std(values, axis=1)
        return np.asarray(std, dtype=np.float64).reshape((-1, 1, 1))
    return np.asarray([[float(np.std(values))]], dtype=np.float64)


def _field_modes(array: np.ndarray, field: str, regions, lat, lon, weights, keep_basis):
    if np.ndim(array) != 4:
        raise ValueError("campo %r deve ter forma (M, nlev, nlat, nlon); recebido %s"
                         % (field, np.shape(array)))
    M, K, nlat, nlon = array.shape
    if (nlat, nlon) != (len(lat), len(lon)):
        raise ValueError("campo %r tem grid %dx%d; lat/lon tem %dx%d"
                         % (field, nlat, nlon, len(lat), len(lon)))
    out = {}
    w = _field_weights(weights, field)
    if w is not None:
        w = np.broadcast_to(np.asarray(w, dtype=np.float64), (K, nlat, nlon))
    for region in regions:
        mask = mb.region_mask(np.asarray(lat)[:, None], np.asarray(lon)[None, :], region)
        npts = int(mask.sum())
        if npts == 0:
            raise ValueError("regiao %r sem pontos na grade (%dx%d)" % (region, nlat, nlon))
        X = np.empty((M, K * npts), dtype=array.dtype)
        for t in range(M):
            X[t] = array[t][:, mask].reshape(-1)
        wflat = None if w is None else w[:, mask].reshape(-1)
        mode, evals, basis = _eof_first(X, wflat, keep_basis)
        mode_grid = np.zeros((K, nlat, nlon), dtype=array.dtype)
        mode_grid[:, mask] = mode.reshape(K, npts)
        out[region] = Mode(
            field=field, region=region, mode=mode_grid, eigenvalues=evals,
            sigma_i=_sigma_i_for(mode_grid, field, mask), eigenvectors=basis,
        )
    return out


def _member_modes(arrays: Mapping[str, np.ndarray], preset, lat, lon, weights, keep_basis):
    regions, fields = mb.preset(preset)
    out = {}
    for f in ("P", "T", "Q"):
        if f in fields and f in arrays:
            out[f] = _field_modes(arrays[f], f, regions, lat, lon, weights, keep_basis)
    if ("U" in fields or "V" in fields) and "wind" in arrays:
        out["wind"] = _field_modes(arrays["wind"], "wind", regions, lat, lon, weights, keep_basis)
    if not out:
        raise ValueError("nenhum campo do preset %r presente em diffs" % (preset,))
    return out


def _assemble_eofs(array_iter: Iterable[Mapping[str, np.ndarray]], preset, lat, lon,
                   weights, keep_basis, source=None) -> dict:
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    members = []
    m_times = None
    kmax = None
    for arrays in array_iter:
        if m_times is None:
            m_times = int(next(iter(arrays.values())).shape[0])
            kmax = int(max(a.shape[1] for a in arrays.values()))
        for field, arr in arrays.items():
            expected_k = 1 if field == "P" else kmax
            if (np.ndim(arr) != 4 or arr.shape[0] != m_times
                    or arr.shape[1] != expected_k):
                raise ValueError(
                    "membro com %r de forma %s; esperado (M=%d, nlev=%d, nlat, nlon)"
                    % (field, np.shape(arr), m_times, expected_k))
        members.append(_member_modes(arrays, preset, lat, lon, weights, keep_basis))
    if not members:
        raise ValueError("nenhum membro em diffs")
    if m_times is None or m_times < 2:
        raise ValueError("o MB09 precisa de M >= 2 tempos (recebido %r)" % (m_times,))
    eof_fields = [f for f in _EOF_FIELDS if f in members[0]]
    return {
        "preset": mb.DEFAULT_PRESET if preset is None else str(preset),
        "regions": list(mb.regions_of(preset)),
        "fields": list(mb.fields_of(preset)),
        "eof_fields": eof_fields,
        "n": len(members),
        "m_times": m_times,
        "kmax": kmax,
        "nlat": int(lat.size),
        "nlon": int(lon.size),
        "lat": lat,
        "lon": lon,
        "members": members,
        "meta": {} if source is None else {"source": str(source)},
    }


def _first_state(member) -> State:
    if isinstance(member, Mapping):
        for value in member.values():
            states = list(value)
            if states:
                return states[0]
        raise ValueError("membro sem estados")
    states = list(member)
    if not states:
        raise ValueError("membro sem estados")
    return states[0]


def build_eofs(diffs, preset: str = mb.DEFAULT_PRESET, weights=None,
               *, keep_basis: bool = False) -> dict:
    """Calcula o 1o EOF por membro/campo/regiao a partir das series de diferencas.

    Parameters
    ----------
    diffs:
        {membro: {campo: [State x M]}} com campo em P/T/Q/U/V (U e V formam o
        vento complexo), ou {membro: [State x M]} usando todos os campos.
        A ordem dos membros e' a ordem de insercao (k = 1..n).
    preset:
        nome do preset de 'bamci.mb09_tables' (default ETSA). Campos ausentes em
        'diffs' sao simplesmente ignorados.
    weights:
        None (default, sem pesos, D3) ou array de pesos metricos broadcastable
        para (nlev, nlat, nlon); dict por campo tambem e' aceito.
    keep_basis:
        guarda tambem a matriz de autovetores (M, M) em Mode.eigenvectors.

    Returns
    -------
    dict com preset, regions, fields, eof_fields, n, m_times, kmax, nlat, nlon,
    lat, lon e 'members' (lista de {campo: {regiao: Mode}}).
    """
    keys = list(diffs.keys())
    if not keys:
        raise ValueError("diffs vazio")
    first = _first_state(diffs[keys[0]])
    lat, lon = first.lat.copy(), first.lon.copy()
    for key in keys:
        state = _first_state(diffs[key])
        if state.lat.shape != lat.shape or not np.allclose(
                state.lat, lat, atol=1e-6, rtol=0.0):
            raise ValueError("membro %r com latitudes diferentes do primeiro" % (key,))
        if state.lon.shape != lon.shape or not np.allclose(
                state.lon, lon, atol=1e-6, rtol=0.0):
            raise ValueError("membro %r com longitudes diferentes do primeiro" % (key,))

    def gen():
        for key in keys:
            yield _member_arrays(diffs[key])

    return _assemble_eofs(gen(), preset, lat, lon, weights, keep_basis)


# ---------------------------------------------------------------------------
# aplicacao dos modos
# ---------------------------------------------------------------------------
def _targets(nlev: int, sigma=None) -> Dict[str, np.ndarray]:
    target = {
        "P": float(mb.AMPLITUDES_TARGET["P"]),
        "T": float(mb.AMPLITUDES_TARGET["T"]),
        "wind": float(mb.AMPLITUDES_TARGET["wind"]),
        "Q": np.asarray(mb.q_stddev_for_layers(nlev), dtype=np.float64),
    }
    if sigma:
        for key, value in sigma.items():
            k = _sigma_entry(key)
            if k == "wind":
                target["wind"] = float(value)
            elif k == "t":
                target["T"] = float(value)
            elif k == "ps":
                target["P"] = float(value)
            elif k in ("u", "v"):
                raise KeyError("sigma de vento deve ser informado como 'wind'")
            else:
                arr = np.asarray(value, dtype=np.float64)
                if arr.ndim == 0:
                    arr = np.full(nlev, float(arr))
                arr = np.broadcast_to(arr, (nlev,)).astype(np.float64)
                target["Q"] = arr
    return target


def _align_mode(mode: np.ndarray, eof_lat, eof_lon, lat, lon) -> np.ndarray:
    """Alinha um modo ao grid do controle (latitude invertida e/ou lon girada)."""
    eof_lat = np.asarray(eof_lat, dtype=np.float64)
    eof_lon = np.asarray(eof_lon, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    if mode.shape[-2:] != (lat.size, lon.size):
        raise ValueError(
            "modo tem grid %s; controle tem (%d, %d)" % (mode.shape[-2:], lat.size, lon.size)
        )
    out = mode
    # o .ctl do PostGrib declara ydef com 3 casas decimais (F10.3), portanto a
    # comparacao precisa tolerar ~1e-3 grau; o que importa e' a ordem dos pontos.
    if not np.allclose(eof_lat, lat, atol=1e-2, rtol=0.0):
        if np.allclose(eof_lat[::-1], lat, atol=1e-2, rtol=0.0):
            out = out[..., ::-1, :]
        else:
            raise ValueError("latitudes do EOF nao batem com o controle")
    if not np.allclose(eof_lon, lon, atol=1e-4, rtol=0.0):
        step = 360.0 / lon.size
        eof_mod = np.mod(eof_lon, 360.0)
        lon_mod = np.mod(lon, 360.0)
        shift = int(round((lon_mod[0] - eof_mod[0]) / step)) % lon.size
        rolled = np.roll(eof_mod, -shift)
        if np.allclose(rolled, lon_mod, atol=1e-4, rtol=0.0):
            out = np.roll(out, -shift, axis=-1)
        elif np.allclose(rolled[::-1], lon_mod, atol=1e-4, rtol=0.0):
            out = np.roll(out, -shift, axis=-1)[..., ::-1]
        else:
            raise ValueError("longitudes do EOF nao batem com o controle")
    return out


def apply_modes(control: State, eofs: Mapping[str, object], preset: str = mb.DEFAULT_PRESET,
                sigma=None) -> List[State]:
    """Gera os 2n estados x_c +/- p_k, na ordem +p_1, -p_1, ..., +p_n, -p_n.

    'preset' serve para validar que o eofs foi construido com o mesmo preset;
    None (ou o default ETSA) usa o preset gravado em 'eofs'. 'sigma' sobrescreve
    as amplitudes alvo (chaves wind, t, ps, q; q aceita vetor por camada).
    """
    eof_preset = eofs.get("preset", preset)
    if preset is None or str(preset) == eof_preset:
        preset = eof_preset
    elif str(preset) == mb.DEFAULT_PRESET:
        # o default da assinatura nao deve impedir usar um eofs de outro preset
        preset = eof_preset
    else:
        raise ValueError("preset %r nao e' o do eofs (%r)" % (preset, eof_preset))
    members = eofs["members"]
    targets = _targets(control.nlev, sigma)
    lat = eofs["lat"]
    lon = eofs["lon"]

    out: List[State] = []
    for member in members:
        delta = {
            "t": np.zeros_like(control.t),
            "q": np.zeros_like(control.q),
            "u": np.zeros_like(control.u),
            "v": np.zeros_like(control.v),
            "ps": np.zeros_like(control.ps),
        }
        for field, regions in member.items():
            target = targets[field]
            for region, mode in regions.items():
                scaled = mode.scaled(target)
                scaled = _align_mode(scaled, lat, lon, control.lat, control.lon)
                if field == "P":
                    delta["ps"] += scaled[0]
                elif field == "T":
                    delta["t"] += scaled
                elif field == "Q":
                    delta["q"] += scaled
                elif field == "wind":
                    delta["u"] += scaled.real
                    delta["v"] += scaled.imag
        out.append(perturbed_state(control, delta))
        out.append(perturbed_state(control, {k: -v for k, v in delta.items()}))
    return out


# ---------------------------------------------------------------------------
# escrita do ensemble
# ---------------------------------------------------------------------------
def write_ensemble(control_path: Union[str, Path], states: Sequence[State],
                   out_dir: Union[str, Path], prefix: str = "mb09", *,
                   suffix: str = "") -> List[Path]:
    """Grava controle + estados (15 arquivos no caso 2n=14) e members.json.

    Nomes: '<prefix>_000.unf<suffix>' (controle, copia do arquivo original),
    '<prefix>_001.unf<suffix>' (+p_1), '_002' (-p_1), ... cada State e' gravado
    com o cabecalho do IC de controle.
    Retorna [controle, +p_1, -p_1, ..., +p_n, -p_n].
    """
    control_path = Path(control_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths: List[Path] = []
    control_out = out_dir / ("%s_000.unf%s" % (prefix, suffix))
    if control_out.resolve() != control_path.resolve():
        shutil.copyfile(control_path, control_out)
    paths.append(control_out)

    for idx, state in enumerate(states, start=1):
        target = out_dir / ("%s_%03d.unf%s" % (prefix, idx, suffix))
        state.to_ic_file(target, template=control_path)
        paths.append(target)

    members = []
    for idx in range(1, len(states) + 1):
        k = (idx + 1) // 2
        sign = "+" if idx % 2 == 1 else "-"
        members.append({"index": idx, "k": k, "sign": sign, "path": str(paths[idx])})
    manifest = {
        "control": str(paths[0]),
        "order": "000=controle; 2k-1=+p_k; 2k=-p_k",
        "files": [str(p) for p in paths],
        "members": members,
    }
    (out_dir / "members.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return paths


def _diff_state(a: State, b: State) -> State:
    if a.shape != b.shape:
        raise ValueError("estados com formas diferentes: %s e %s" % (a.shape, b.shape))
    return State(
        t=a.t - b.t, q=a.q - b.q, u=a.u - b.u, v=a.v - b.v, ps=a.ps - b.ps,
        lat=a.lat.copy(), lon=a.lon.copy(),
    )


def _normalize_series(member_ics) -> Dict[str, list]:
    if isinstance(member_ics, Mapping):
        if "control" in member_ics and "members" in member_ics:
            return {
                "control": list(member_ics["control"]),
                "members": [list(series) for series in member_ics["members"]],
            }
        raise KeyError("member_ics mapping precisa de 'control' e 'members'")
    seq = list(member_ics)
    if not seq:
        raise ValueError("member_ics vazio")
    return {"control": list(seq[0]), "members": [list(series) for series in seq[1:]]}


def run_mb09(control_ic: Union[str, Path], member_ics, out_dir: Union[str, Path],
             preset: str = mb.DEFAULT_PRESET, seed: int = 0, *,
             prefix: str = "mb09", suffix: str = "", sigma=None) -> List[Path]:
    """Build + apply a partir das series ja integradas.

    'control_ic' e' o IC da analise de controle (estado base a perturbar).
    'member_ics' deve ser {"control": [M caminhos], "members": [[M caminhos] x n]}
    ou uma sequencia de n+1 series (controle primeiro). Cada serie tem M estados
    em formato IC (State.from_ic_file). 'seed' e' apenas proveniencia aqui: a
    integracao de 36 h e' externa.

    Retorna os 15 caminhos gravados por 'write_ensemble'.
    """
    control = State.from_ic_file(control_ic)
    series = _normalize_series(member_ics)
    n_members = len(series["members"])
    if n_members == 0:
        raise ValueError("nenhuma serie de membro em member_ics")
    m_times = len(series["control"])
    ctrl = [State.from_ic_file(p) for p in series["control"]]

    diffs = {}
    for k, paths in enumerate(series["members"], start=1):
        if len(paths) != m_times:
            raise ValueError(
                "membro %d tem %d tempos; controle tem %d" % (k, len(paths), m_times))
        states = [State.from_ic_file(p) for p in paths]
        diffs[k] = [_diff_state(s, c) for s, c in zip(states, ctrl)]
    eofs = build_eofs(diffs, preset)
    eofs["meta"]["seed"] = int(seed)
    states = apply_modes(control, eofs, preset, sigma)
    return write_ensemble(control_ic, states, out_dir, prefix, suffix=suffix)


# ---------------------------------------------------------------------------
# leitura de series GrADS/PostGrib (.ctl com M tempos no mesmo binario)
# ---------------------------------------------------------------------------
_GRADS_CANDIDATES = {
    "P": ("pslc", "ps", "psfc", "prmsl", "press", "pres", "slp", "pressao"),
    "T": ("tvir", "t", "temp", "tmp", "temperatura"),
    "Q": ("umes", "q", "sphu", "umid", "umidade", "hu"),
    "U": ("uvel", "u", "uwnd", "uwind", "vento_zonal"),
    "V": ("vvel", "v", "vwnd", "vwind", "vento_meridional"),
}


_VAR_KEY_ALIASES = {
    "PS": "P", "TVIR": "T", "TEMP": "T", "T": "T",
    "SPHU": "Q", "UMES": "Q", "Q": "Q",
    "UVEL": "U", "U": "U", "VVEL": "V", "V": "V",
    "PSLC": "P", "PRMSL": "P", "P": "P",
}


def _normalize_var_key(key) -> str:
    k = str(key).strip().upper()
    if k in _VAR_KEY_ALIASES:
        return _VAR_KEY_ALIASES[k]
    if k in _GRADS_CANDIDATES:
        return k
    raise KeyError("campo %r desconhecido (use P/T/Q/U/V)" % (key,))


def _resolve_grads_vars(ctl, variables=None) -> Dict[str, str]:
    """Resolve os 5 campos na vlist; 'variables' sobrescreve caso a caso."""
    names = {v.name.lower(): v.name for v in ctl.variables}
    out: Dict[str, str] = {}
    unresolved = []
    for key, candidates in _GRADS_CANDIDATES.items():
        for cand in candidates:
            if cand in names:
                out[key] = names[cand]
                break
        else:
            unresolved.append(key)
    if variables:
        for key, value in variables.items():
            out[_normalize_var_key(key)] = str(value)
    for key, value in list(out.items()):
        if str(value).lower() not in names:
            raise KeyError(
                "variavel %r ausente no .ctl; disponiveis: %s"
                % (value, sorted(v.name for v in ctl.variables)))
        out[key] = names[str(value).lower()]
    missing = [k for k in ("P", "T", "Q", "U", "V") if k not in out]
    if missing:
        raise KeyError(
            "nao encontrei variavel para %s no .ctl; disponiveis: %s; "
            "use variables={...} ou --vars P=...,T=..."
            % (missing, sorted(v.name for v in ctl.variables)))
    return out


# ---------------------------------------------------------------------------
# janela temporal (6..36 h, passo 3 h => M=11)
# ---------------------------------------------------------------------------
def _format_hours(values) -> str:
    return ", ".join("%g" % float(v) for v in np.asarray(values, dtype=np.float64).ravel())


def _select_hours(times_hours, wanted, label: str = "janela") -> np.ndarray:
    """Indices (0-based, em ordem cronologica) de 'wanted' dentro de 'times_hours'.

    Erro claro com a lista de tempos disponiveis quando falta algum horario.
    """
    times = np.asarray(times_hours, dtype=np.float64).ravel()
    wanted = np.asarray(wanted, dtype=np.float64).ravel()
    idx: List[int] = []
    missing: List[float] = []
    for hour in wanted:
        hit = np.flatnonzero(np.isclose(times, hour, atol=1e-6, rtol=0.0))
        if hit.size == 0:
            missing.append(float(hour))
        else:
            idx.append(int(hit[0]))
    if missing:
        raise ValueError(
            "%s pede os tempos %s h, ausentes no .ctl; tempos disponiveis (h): %s"
            % (label, _format_hours(missing), _format_hours(times)))
    order = np.argsort(times[idx]) if idx else np.asarray([], dtype=int)
    return np.asarray(idx, dtype=int)[order]


def select_window(times_hours, start: float = 6.0, stop: float = 36.0,
                  step: float = 3.0) -> np.ndarray:
    """Indices de times_hours que caem na janela start..stop (passo step).

    Default do MB09: 6, 9, ..., 36 h (M=11). Tempos fora da janela sao
    ignorados pelo chamador; horarios da janela ausentes geram ValueError com
    a lista dos tempos disponiveis.
    """
    start, stop, step = float(start), float(stop), float(step)
    if step <= 0.0:
        raise ValueError("passo da janela deve ser > 0 (recebido %r)" % (step,))
    if stop < start:
        raise ValueError("fim da janela (%r) menor que o inicio (%r)" % (stop, start))
    count = int(round((stop - start) / step))
    wanted = start + step * np.arange(count + 1, dtype=np.float64)
    return _select_hours(times_hours, wanted, label="janela %g:%g:%g" % (start, stop, step))


def select_times(times_hours, wanted_hours) -> np.ndarray:
    """Indices de uma lista explicita de horas (ex.: --times 6,9,...,36)."""
    return _select_hours(times_hours, wanted_hours, label="--times")


def _parse_analysis_time(value) -> Optional[datetime]:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value
    text = str(value).strip()
    try:
        return datetime.strptime(text, "%Y%m%d%H")
    except ValueError:
        raise ValueError(
            "--analysis-time espera YYYYMMDDHH (ex.: 2026010100); recebido %r" % (value,)
        ) from None


def _infer_analysis_time(times: Sequence[datetime], analysis_time=None) -> datetime:
    """Referencia para converter tdef em horas de previsao.

    Ordem: 'analysis_time' explicito; senao 00Z do primeiro t do .ctl (o
    PostGrib grava o tdef a partir de 00Z do dia da rodada; um .ctl que comece
    as 06Z tambem cai certo nessa referencia).
    """
    parsed = _parse_analysis_time(analysis_time)
    if parsed is not None:
        return parsed
    if not times:
        raise ValueError(".ctl sem tempos em tdef")
    first = times[0]
    return datetime(first.year, first.month, first.day)


class _GradsSeries:
    """Le um .ctl/.unf do PostGrib com os M tempos no mesmo binario.

    Implementacao local (nao usa GradsFile) porque variaveis 2-D declaradas com
    '0 99' ocupam 1 registro, enquanto GradsFile assume nz registros e
    desalinha as variaveis seguintes. Convencoes: linhas do binario sao
    invertidas quando o .ctl tem 'yrev'; o resultado fica na ordem de
    ctl.ydef.values(). 'analysis_time' (datetime ou "YYYYMMDDHH") e' a
    referencia de 'self.hours'; sem ela usa-se 00Z do primeiro tempo do .ctl.
    """

    def __init__(self, ctl_path, variables=None, ps_units: str = "hpa",
                 analysis_time=None):
        from .ctl import read_ctl  # import local: mantem o modulo leve

        self.ctl = read_ctl(ctl_path)
        self.path = Path(ctl_path)
        self.byteorder = "little" if "little_endian" in self.ctl.options else "big"
        self.dtype = np.dtype("<f4" if self.byteorder == "little" else ">f4")
        self.base = self.ctl.path.parent if self.ctl.path is not None else Path(".")
        self.nx = int(self.ctl.nx)
        self.ny = int(self.ctl.ny)
        self.nt = int(self.ctl.nt)
        self.times = list(self.ctl.tdef.times())
        self.analysis_time = _infer_analysis_time(self.times, analysis_time)
        self.hours = np.asarray(
            [(t - self.analysis_time).total_seconds() / 3600.0 for t in self.times],
            dtype=np.float64)
        self.lat = np.asarray(self.ctl.ydef.values(), dtype=np.float64)
        self.lon = np.asarray(self.ctl.xdef.values(), dtype=np.float64)
        self.vars = _resolve_grads_vars(self.ctl, variables)
        self.ps_units = str(ps_units).lower()
        if self.ps_units not in ("hpa", "pa"):
            raise ValueError("ps_units deve ser 'hpa' ou 'pa'")
        if "template" in self.ctl.options:
            raise NotImplementedError("dset com template nao e' suportado no build")

        self._layout: Dict[str, Tuple[int, int]] = {}
        cell = 4 * self.nx * self.ny
        offset = 0
        for var in self.ctl.variables:
            nlev = len(var.levels) if var.levels else 1
            self._layout[var.name] = (offset, nlev)
            offset += nlev * cell
        self.time_block = offset
        self._cell = cell

    def _file(self) -> Path:
        dset = self.ctl.dset
        return self.base / dset[1:].strip() if dset.startswith("^") else Path(dset)

    def _read(self, name: str, t: int, z: int = 1) -> np.ndarray:
        if not (1 <= t <= self.nt):
            raise IndexError("t=%d fora de 1..%d" % (t, self.nt))
        st_off, nlev = self._layout[name]
        if not (1 <= z <= nlev):
            raise IndexError("z=%d fora de 1..%d para %s" % (z, nlev, name))
        offset = (t - 1) * self.time_block + st_off + (z - 1) * self._cell
        field = self._file()
        need = offset + self._cell
        size = field.stat().st_size
        if size < need:
            raise ValueError("%s tem %d bytes; leitura pede %d" % (field, size, need))
        with field.open("rb") as fh:
            fh.seek(offset)
            data = fh.read(self._cell)
        raw = np.frombuffer(data, dtype=self.dtype).astype(np.float64)
        raw = raw.reshape(self.ny, self.nx)
        if "yrev" in self.ctl.options:
            raw = raw[::-1, :]
        bad = ~np.isfinite(raw) | (np.abs(raw) >= 1.0e30)
        if bad.any():
            raise ValueError("%s: %d valores indefinidos no campo %r (t=%d)"
                             % (field, int(bad.sum()), name, t))
        return raw

    def _volume(self, name: str, t: int) -> np.ndarray:
        _, nlev = self._layout[name]
        return np.stack([self._read(name, t, z) for z in range(1, nlev + 1)], axis=0)

    def state(self, t: int = 1) -> State:
        ps = self._read(self.vars["P"], t, 1)
        if self.ps_units == "pa":
            ps = ps / 100.0
        fields = {key: self._volume(name, t) for key, name in self.vars.items()
                  if key != "P"}
        shapes = {f.shape for f in fields.values()}
        if len(shapes) != 1:
            raise ValueError("variaveis 3-D com formas diferentes: %s" % (sorted(shapes),))
        return State(
            t=fields["T"], q=fields["Q"], u=fields["U"], v=fields["V"], ps=ps,
            lat=self.lat.copy(), lon=self.lon.copy(),
        )


def build_eofs_from_files(control_ctl: Union[str, Path],
                          member_ctls: Sequence[Union[str, Path]],
                          preset: str = mb.DEFAULT_PRESET, *,
                          variables=None, ps_units: str = "hpa",
                          weights=None, keep_basis: bool = False,
                          hours=None, window=None, analysis_time=None) -> dict:
    """Como build_eofs, mas lendo .ctl/.unf do PostGrib (grade, M tempos).

    Cada serie (controle e n membros) e' um unico .ctl com nt tempos. A janela
    do MB09 (6..36 h, passo 3 h => M=11) e' selecionada explicitamente:

      * 'window=(start, stop, step)' (default (6, 36, 3));
      * 'hours=[...]' sobrepoe a janela com uma lista de horas;
      * 'analysis_time' (datetime ou "YYYYMMDDHH") e' a referencia para
        converter o tdef em horas de previsao; sem ela usa-se 00Z do primeiro
        tempo do .ctl.

    Tempos fora da janela sao IGNORADOS (o total fica em eofs["meta"]
    ["n_ignored_times"] e as horas selecionadas em eofs["meta"]["hours"] e
    ["analysis_time"]); horarios da janela ausentes geram ValueError com a
    lista de tempos disponiveis. Os campos sao resolvidos pela vlist do .ctl
    (pslc/tvir/umes/uvel/vvel por default) ou por 'variables={"P": "ps", ...}'.
    """
    control = _GradsSeries(control_ctl, variables, ps_units, analysis_time=analysis_time)
    if control.nt < 2:
        raise ValueError("controle tem %d tempos; o MB09 precisa de M >= 2" % (control.nt))
    if hours is not None:
        idx = select_times(control.hours, hours)
    elif window is not None:
        if len(tuple(window)) != 3:
            raise ValueError("window deve ser (inicio, fim, passo)")
        idx = select_window(control.hours, *window)
    else:
        idx = select_window(control.hours)
    selected = control.hours[idx]
    n_ignored = int(control.nt - idx.size)
    if idx.size < 2:
        raise ValueError(
            "a janela selecionou %d tempo(s); o MB09 precisa de M >= 2 "
            "(horas disponiveis: %s)" % (idx.size, _format_hours(control.hours)))
    ctrl = [control.state(int(t) + 1) for t in idx]

    def gen():
        for path in member_ctls:
            series = _GradsSeries(path, variables, ps_units,
                                  analysis_time=control.analysis_time)
            if series.nt != control.nt or not np.allclose(
                    series.hours, control.hours, atol=1e-6, rtol=0.0):
                raise ValueError(
                    "%s: tempos diferentes do controle (nt=%d vs %d)"
                    % (path, series.nt, control.nt))
            if series.lat.shape != control.lat.shape or not np.allclose(
                    series.lat, control.lat, atol=1e-6, rtol=0.0):
                raise ValueError("%s: latitudes diferentes do controle" % (path,))
            if series.lon.shape != control.lon.shape or not np.allclose(
                    series.lon, control.lon, atol=1e-6, rtol=0.0):
                raise ValueError("%s: longitudes diferentes do controle" % (path,))
            diff = [_diff_state(series.state(int(t) + 1), ctrl[j])
                    for j, t in enumerate(idx)]
            yield _member_arrays(diff)

    eofs = _assemble_eofs(gen(), preset, control.lat, control.lon, weights,
                          keep_basis, source=control_ctl)
    eofs["meta"].update({
        "hours": [float(h) for h in selected],
        "analysis_time": control.analysis_time.strftime("%Y%m%d%H"),
        "n_ignored_times": n_ignored,
    })
    return eofs


# ---------------------------------------------------------------------------
# persistencia dos EOFs (.npz)
# ---------------------------------------------------------------------------
def save_eofs(path: Union[str, Path], eofs: Mapping[str, object]) -> Path:
    """Serializa o dict de build_eofs em .npz (sem pickle)."""
    meta = {
        "preset": eofs["preset"],
        "regions": list(eofs["regions"]),
        "fields": list(eofs["fields"]),
        "eof_fields": list(eofs["eof_fields"]),
        "n": int(eofs["n"]),
        "m_times": int(eofs["m_times"]),
        "kmax": int(eofs["kmax"]),
        "nlat": int(eofs["nlat"]),
        "nlon": int(eofs["nlon"]),
    }
    payload: Dict[str, np.ndarray] = {
        "meta": np.asarray(json.dumps(meta)),
        "lat": np.asarray(eofs["lat"], dtype=np.float64),
        "lon": np.asarray(eofs["lon"], dtype=np.float64),
    }
    for i, member in enumerate(eofs["members"]):
        for field, regions in member.items():
            for region, mode in regions.items():
                key = "m%d/%s/%s" % (i, field, region)
                payload[key + "::mode"] = np.asarray(mode.mode)
                payload[key + "::eval"] = np.asarray(mode.eigenvalues, dtype=np.float64)
                payload[key + "::sigi"] = np.asarray(mode.sigma_i, dtype=np.float64)
                if mode.eigenvectors is not None:
                    payload[key + "::evec"] = np.asarray(mode.eigenvectors)
    path = Path(path)
    with path.open("wb") as fh:   # file handle: nao deixa o numpy acrescentar .npz
        np.savez(fh, **payload)
    return path


def load_eofs(path: Union[str, Path]) -> dict:
    """Le o .npz gravado por 'save_eofs'."""
    with np.load(path, allow_pickle=False) as data:
        if "meta" not in data:
            raise ValueError("%s nao e' um eof.npz do mb09" % (path,))
        meta = json.loads(str(data["meta"]))
        lat = np.asarray(data["lat"], dtype=np.float64)
        lon = np.asarray(data["lon"], dtype=np.float64)
        modes: Dict[Tuple[int, str, str], Mode] = {}
        for key in data.files:
            if not key.endswith("::mode"):
                continue
            base = key[: -len("::mode")]
            try:
                idx, field, region = base[1:].split("/")
                idx = int(idx)
            except ValueError:
                raise ValueError("chave invalida em %s: %r" % (path, base))
            modes[(idx, field, region)] = Mode(
                field=field, region=region,
                mode=np.asarray(data[base + "::mode"]),
                eigenvalues=np.asarray(data[base + "::eval"], dtype=np.float64),
                sigma_i=np.asarray(data[base + "::sigi"], dtype=np.float64),
                eigenvectors=(np.asarray(data[base + "::evec"])
                              if base + "::evec" in data.files else None),
            )

    members = []
    for i in range(int(meta["n"])):
        member = {}
        for field in meta["eof_fields"]:
            regions = {}
            for region in meta["regions"]:
                mode = modes.get((i, field, region))
                if mode is not None:
                    regions[region] = mode
            if regions:
                member[field] = regions
        if not member:
            raise ValueError("%s: membro %d sem modos" % (path, i))
        members.append(member)
    return {
        "preset": meta["preset"], "regions": meta["regions"],
        "fields": meta["fields"], "eof_fields": meta["eof_fields"],
        "n": meta["n"], "m_times": meta["m_times"], "kmax": meta["kmax"],
        "nlat": meta["nlat"], "nlon": meta["nlon"], "lat": lat, "lon": lon,
        "members": members, "meta": {},
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_vars(text: Optional[str]) -> Optional[Dict[str, str]]:
    if not text:
        return None
    out = {}
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError("--vars espera campo=variavel, recebido %r" % (part,))
        key, value = part.split("=", 1)
        out[_normalize_var_key(key)] = value.strip()
    return out


def _parse_window(text) -> Optional[Tuple[float, float, float]]:
    """Interpreta --window inicio:fim:passo (aceita ':' ou ',', em horas)."""
    if not text:
        return None
    parts = str(text).replace(":", " ").replace(",", " ").split()
    if len(parts) != 3:
        raise ValueError("--window espera inicio:fim:passo (ex. 6:36:3); recebido %r"
                         % (text,))
    start, stop, step = (float(p) for p in parts)
    if step <= 0.0:
        raise ValueError("--window com passo <= 0 (%r)" % (step,))
    if stop < start:
        raise ValueError("--window com fim < inicio (%r < %r)" % (stop, start))
    return (start, stop, step)


def _parse_hours(text) -> Optional[List[float]]:
    """Interpreta --times 6,9,12,...,36 (horas de previsao)."""
    if not text:
        return None
    parts = str(text).replace(":", ",").replace(" ", ",").split(",")
    values = [float(p) for p in parts if p.strip()]
    if not values:
        raise ValueError("--times vazio")
    return values


def _parse_sigma(text: Optional[str]) -> Optional[Dict[str, object]]:
    if not text:
        return None
    out: Dict[str, object] = {}
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError("--sigma espera campo=valor, recebido %r" % (part,))
        key, value = part.split("=", 1)
        values = [float(v) for v in value.replace(" ", "").split("/") if v]
        out[key.strip()] = values[0] if len(values) == 1 else np.asarray(values)
    return out


def _cmd_random(args) -> int:
    control = State.from_ic_file(args.ic)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    suffix = (".%s" % args.resolution) if args.resolution else ""
    paths = []
    for k in range(1, int(args.n) + 1):
        state = perturbed_state(control, random_deltas(control, rng))
        path = out_dir / ("%s_%02d.unf%s" % (args.prefix, k, suffix))
        state.to_ic_file(path, template=args.ic)
        paths.append(path)
    manifest = {
        "control": str(Path(args.ic)),
        "members": [str(p) for p in paths],
        "seed": int(args.seed),
        "preset": args.preset,
        "n": int(args.n),
    }
    (out_dir / "random.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    for path in paths:
        print(path)
    return 0


def _discover_build_inputs(args):
    work = Path(args.work) if args.work else None
    control = args.control
    members = args.members
    if control is None and work is not None:
        control = str(work / "mb09" / "forecast" / "control.ctl")
    if not members and work is not None:
        members = [str(p) for p in sorted((work / "mb09" / "forecast").glob("member_*.ctl"))]
    if control is None:
        raise SystemExit("build: informe --control <controle.ctl> ou --work")
    if not members:
        raise SystemExit("build: informe --members <m1.ctl> ... <m7.ctl>")
    return control, list(members)


def _cmd_build(args) -> int:
    control, members = _discover_build_inputs(args)
    hours = _parse_hours(args.times)
    window = None if hours is not None else _parse_window(args.window)
    eofs = build_eofs_from_files(
        control, members, preset=args.preset,
        variables=_parse_vars(args.vars), ps_units=args.ps_units,
        hours=hours, window=window, analysis_time=args.analysis_time,
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    save_eofs(out, eofs)
    print("eof: %s (n=%d, M=%d, preset=%s, campos=%s, horas=%s)"
          % (out, eofs["n"], eofs["m_times"], eofs["preset"],
             ",".join(eofs["eof_fields"]),
             ",".join("%g" % h for h in eofs["meta"].get("hours", []))))
    ignored = int(eofs["meta"].get("n_ignored_times", 0))
    if ignored:
        print("aviso: %d tempo(s) fora da janela ignorado(s) (use --window/--times)"
              % (ignored,))
    return 0


def _cmd_apply(args) -> int:
    eofs = load_eofs(args.eof)
    control = State.from_ic_file(args.ic)
    preset = args.preset or eofs["preset"]
    states = apply_modes(control, eofs, preset, _parse_sigma(args.sigma))
    suffix = (".%s" % args.resolution) if args.resolution else ""
    paths = write_ensemble(args.ic, states, args.out_dir, args.prefix, suffix=suffix)
    for path in paths:
        print(path)
    return 0


def build_parser():
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m bamci.mb09", description="Metodo MB09 (EOF) para o ensemble do BAM")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("random", help="gera n ICs com perturbacoes aleatorias")
    p.add_argument("--ic", required=True, help="IC da analise de controle")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n", type=int, default=mb.N_PERTURBATIONS)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--preset", default=mb.DEFAULT_PRESET)
    p.add_argument("--resolution", default=None, help="sufixo do nome, ex. TQ0126L042")
    p.add_argument("--prefix", default="mb09rand")
    p.set_defaults(func=_cmd_random)

    p = sub.add_parser("build", help="calcula EOFs das series controle + membros")
    p.add_argument("--work", default=None, help="usa <work>/mb09/forecast/*.ctl")
    p.add_argument("--control", default=None, help=".ctl do controle (PostGrib)")
    p.add_argument("--members", nargs="+", default=None, help=".ctl dos n membros")
    p.add_argument("--preset", default=mb.DEFAULT_PRESET)
    p.add_argument("--vars", default=None, help="P=pslc,T=tvir,Q=umes,U=uvel,V=vvel")
    p.add_argument("--ps-units", default="hpa", choices=("hpa", "pa"))
    p.add_argument("--window", default="6:36:3",
                   help="janela inicio:fim:passo em horas (default 6:36:3 => M=11)")
    p.add_argument("--times", default=None,
                   help="horas explicitas, ex. 6,9,...,36 (sobrepoe --window)")
    p.add_argument("--analysis-time", default=None,
                   help="YYYYMMDDHH da analise (default: 00Z do 1o tempo do .ctl)")
    p.add_argument("--out", required=True)
    p.set_defaults(func=_cmd_build)

    p = sub.add_parser("apply", help="aplica os EOFs e grava os 15 ICs")
    p.add_argument("--eof", required=True)
    p.add_argument("--ic", required=True, help="IC da analise de controle")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--preset", default=None)
    p.add_argument("--resolution", default=None)
    p.add_argument("--prefix", default="mb09")
    p.add_argument("--sigma", default=None, help="ex. T=1.5,wind=5.0,p=1.0")
    p.set_defaults(func=_cmd_apply)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError, KeyError, IndexError,
            NotImplementedError) as exc:
        import sys
        print("bamci.mb09: erro: %s" % (exc,), file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
