"""Interpolacao vertical das analises para os niveis do modelo.

Porte de ICnAGCM/VerticalInterpolation.f90 (VertSigmaInter) e de
ICnAGCM/Utils.f90 (NewSigma, NewPs). Regras (docs/spec_fortran.md §6-7):

* coordenada de interpolacao z = -ln(p), crescente porque p decresce
  (niveis ordenados superficie -> topo);
* Lagrange cubico de 4 nos no interior, linear nas duas camadas das pontas,
  valor da borda fora do dominio;
* clamp de monotonicidade entre os dois nos internos do intervalo (evita
  overshoot do cubico);
* abaixo do nivel mais baixo (p_out > p_in[0]): temperatura com gradiente
  constante de -6.5 K/km e umidade com UR constante.
"""

from __future__ import annotations

import numpy as np

GRAV = 9.80665
RD = 287.05
RV = 461.50
CP = 1004.6
LC = 2.5e6
GAMMA = -6.5e-3
GORd = GRAV / RD
ROCP = RD / CP
ROCP1 = ROCP + 1.0
DLTDZ = GAMMA * RD / GRAV          # -0.19029 (K/m normalizado)
DLPVDRT = -LC / RV                 # -5417.1

__all__ = ["new_sigma", "new_ps", "interp_logp", "vert_sigma_inter", "surface_fix"]


def new_sigma(delta_sigma) -> tuple:
    """Interfaces e camadas sigma (Utils.f90:18-53)."""
    d = np.asarray(delta_sigma, dtype=np.float64)
    n = d.size
    inter = np.empty(n + 1)
    inter[0] = 1.0
    for k in range(n):
        inter[k + 1] = inter[k] - d[k]
    inter[-1] = 0.0
    s0, s1 = inter[:-1], inter[1:]
    diff = s0 - s1
    safe = np.where(np.abs(diff) > 1e-12, diff, 1.0)
    layer = ((s0 ** ROCP1 - s1 ** ROCP1) / (ROCP1 * safe)) ** (1.0 / ROCP)
    layer = np.where(np.abs(diff) > 1e-12, layer, s0)
    return inter, layer


def new_ps(ps_in, p_levels, tvir, topo_out, g_topo: float = 0.0, g_eps: float = 1e-9) -> np.ndarray:
    """Pressao de superficie consistente com a topografia (Utils.f90:56-158)."""
    ps = np.asarray(ps_in, dtype=np.float64)
    tv = np.asarray(tvir, dtype=np.float64)
    topo = np.asarray(topo_out, dtype=np.float64)
    p = np.asarray(p_levels, dtype=np.float64).ravel()
    if ps.shape != topo.shape or tv.shape[1:] != ps.shape:
        raise ValueError("formas incompativeis: ps%s tv%s topo%s" % (ps.shape, tv.shape, topo.shape))

    g = GAMMA
    zu = g_topo - tv[0] / g * ((ps / p[0]) ** (-g / GORd) - 1.0)
    base0 = np.maximum(1.0 + g / tv[0] * (topo - zu), 1e-6)
    ps_out = np.where(topo <= zu, p[0] * base0 ** (-GORd / g), 0.0)

    for k in range(1, p.size):
        gk = -GORd * np.log(tv[k - 1] / tv[k]) / np.log(p[k - 1] / p[k])
        strong = np.abs(gk) > g_eps
        safe = np.where(strong, gk, 1.0)
        zu = zu + np.where(strong, -tv[k] / safe * ((p[k - 1] / p[k]) ** (-safe / GORd) - 1.0),
                           tv[k] / GORd * np.log(p[k - 1] / p[k]))
        fill = (ps_out == 0.0) & (topo <= zu)
        if np.any(fill):
            safe_all = np.where(strong, gk, 1.0)
            basek = np.maximum(1.0 + safe_all / tv[k] * (topo - zu), 1e-6)
            with np.errstate(over="ignore", invalid="ignore"):
                new = np.where(
                    strong,
                    p[k] * basek ** (-GORd / safe_all),
                    p[k] * np.exp(-GORd / tv[k] * (topo - zu)),
                )
            ps_out = np.where(fill, new, ps_out)

    still = ps_out == 0.0
    if np.any(still):
        ps_out = np.where(still, p[-1] * np.exp(-GORd / tv[-1] * (topo - zu)), ps_out)
    return ps_out


def _gather(a: np.ndarray, i: np.ndarray) -> np.ndarray:
    i = np.clip(i, 0, a.shape[0] - 1)
    return np.take_along_axis(a, i[None, ...], axis=0)[0]


def interp_logp(field, p_in, p_out) -> np.ndarray:
    """Interpola field(K, ...) de p_in(K,) para p_out(K2, ...) em -ln(p)."""
    f = np.asarray(field, dtype=np.float64)
    p1 = np.asarray(p_in, dtype=np.float64).ravel()
    p2 = np.asarray(p_out, dtype=np.float64)
    if f.shape[0] != p1.size:
        raise ValueError("field tem %d niveis, p_in tem %d" % (f.shape[0], p1.size))
    if np.any(p1 <= 0) or np.any(p2 <= 0):
        raise ValueError("pressoes devem ser positivas")
    z1 = -np.log(p1)
    if np.any(np.diff(z1) <= 0):
        raise ValueError("p_in deve decrescer com o nivel (superficie -> topo)")
    z2 = -np.log(p2)
    n = p1.size
    out = np.empty_like(z2)
    for k in range(z2.shape[0]):
        out[k] = _interp_level(f, z1, z2[k])
    return out


def _interp_level(f: np.ndarray, z1: np.ndarray, z: np.ndarray) -> np.ndarray:
    n = z1.size
    idx = np.searchsorted(z1, z) - 1
    i = np.clip(idx, 1, n - 3)

    a, b, c, d = (_gather(f, i - 1), _gather(f, i), _gather(f, i + 1), _gather(f, i + 2))
    za, zb, zc, zd = z1[i - 1], z1[i], z1[i + 1], z1[i + 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        wa = (z - zb) * (z - zc) * (z - zd) / ((za - zb) * (za - zc) * (za - zd))
        wb = (z - za) * (z - zc) * (z - zd) / ((zb - za) * (zb - zc) * (zb - zd))
        wc = (z - za) * (z - zb) * (z - zd) / ((zc - za) * (zc - zb) * (zc - zd))
        wd = (z - za) * (z - zb) * (z - zc) / ((zd - za) * (zd - zb) * (zd - zc))
        out = wa * a + wb * b + wc * c + wd * d
    out = np.minimum(np.maximum(out, np.minimum(b, c)), np.maximum(b, c))

    lin = (idx == 0) | (idx == n - 2)
    if np.any(lin):
        j = np.clip(idx, 0, n - 2)
        f0, f1 = _gather(f, j), _gather(f, j + 1)
        w = (z - z1[j]) / (z1[j + 1] - z1[j])
        out = np.where(lin, f0 * (1.0 - w) + f1 * w, out)

    out = np.where(idx < 0, f[0], out)
    out = np.where(idx >= n - 1, f[n - 1], out)
    return out


def surface_fix(t2, q2, p_in, p_out, t1, q1) -> tuple:
    """Corrige t/q abaixo do nivel mais baixo (VerticalInterpolation.f90:126-176).

    p_out > p_in[0] significa alvo abaixo da superficie da analise:
      * primeira camada: usa os valores de entrada do nivel 1;
      * camadas internas: suavizacao 1-2-1 nos valores ja interpolados;
      * ultima camada: gradiente constante t = t1(1) exp(dltdz*dz) e
        q = q1(1) exp(dlpvdrt (1/t2 - 1/t1(1)) - dz)  (UR constante).
    """
    p_in = np.asarray(p_in, dtype=np.float64).ravel()
    p_out = np.asarray(p_out, dtype=np.float64)
    z1 = -np.log(p_in[0])
    z2 = -np.log(p_out)
    dz = z2 - z1

    below = dz < 0
    n2 = p_out.shape[0]
    t_ref = np.array(t2, dtype=np.float64, copy=True)   # valores interpolados originais
    q_ref = np.array(q2, dtype=np.float64, copy=True)
    t2 = np.array(t2, dtype=np.float64, copy=True)
    q2 = np.array(q2, dtype=np.float64, copy=True)

    if np.any(below):
        t1s = t1[0]
        q1s = q1[0]
        for k in range(n2):
            mask = below[k]
            if not np.any(mask):
                continue
            if k == 0:
                t2[k] = np.where(mask, t1s, t2[k])
                q2[k] = np.where(mask, q1s, q2[k])
            elif k == n2 - 1:
                t_new = t1s * np.exp(DLTDZ * dz[k])
                q_new = q1s * np.exp(DLPVDRT * (1.0 / t_new - 1.0 / t1s) - dz[k])
                t2[k] = np.where(mask, t_new, t2[k])
                q2[k] = np.where(mask, q_new, q2[k])
            else:
                t_s = 0.25 * t_ref[k - 1] + 0.5 * t_ref[k] + 0.25 * t_ref[k + 1]
                q_s = 0.25 * q_ref[k - 1] + 0.5 * q_ref[k] + 0.25 * q_ref[k + 1]
                t2[k] = np.where(mask, t_s, t2[k])
                q2[k] = np.where(mask, q_s, q2[k])
    return t2, q2


def vert_sigma_inter(fields, p_in, p_out) -> dict:
    """Interpola u, v, t, q para os niveis de saida e aplica a correcao de superficie.

    fields: {"u": (K, ...), "v": (K, ...), "t": (K, ...), "q": (K, ...)}
    p_out: (K2, ...) pressao dos niveis de saida em cada coluna.
    """
    out = {name: interp_logp(fields[name], p_in, p_out) for name in ("u", "v")}
    t2 = interp_logp(fields["t"], p_in, p_out)
    q2 = interp_logp(fields["q"], p_in, p_out)
    t2, q2 = surface_fix(t2, q2, p_in, p_out, fields["t"], fields["q"])
    out["t"] = t2
    out["q"] = q2
    return out
