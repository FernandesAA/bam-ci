"""Leitura do historico do BAM 1.2 (.dir + .fct), no lugar do PostGrib.

O PostGrib (pos/source/post.f90 do modelo) le o .fct e reescreve um .unf + .ctl
para o GrADS. Aqui o .fct e' lido direto, sem binario intermediario (o .ctl
equivalente continua disponivel para o resto do pacote, ver write_grads).

Formato confirmado em 2026-11-24 contra as 8 rodadas de 36 h em
/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126 (jaci):

* o tamanho do arquivo fecha exatamente:
  total = 48 + soma_sobre_variaveis nlev*(8 + 4*npts), 8 = dois marcadores de 4 B;
* o primeiro registro tem 40 B = 10 int32, casando com ReadHeader de
  pos/source/FileAccess.f90: ifday, nsec, idate(4), idatec(4);
* a topografia espectral lida do .fct e' BIT-IDENTICA (max|diff| = 0.0) aos
  coeficientes do IC de controle GANLNMC2026111900S.unf.TQ0126L042, o que fixa
  a ordem e a normalizacao dos coeficientes;
* T(k=1) ~ 288 K e q(k=1) ~ 9e-3 (camada 1 = fundo), T(k=42) ~ 255 K (topo).

.dir (texto)
-----------
linha 1: CPTEC SIGMA VERS 2.0
linha 2: 0001 SEQU T0126   127   42   42    3   19   11 2026    0   11   19 2026 TAPE
         [0]=versao do arquivo (0001)  [1]=SEQU  [2]=truncacao (T0126)
         [3]=mwv=Mend+1 (127)  [4]=nz  [5]=nz
         [6..9] = data corrente  (hh, dd, mm, yyyy)   <- dia ANTES do mes
         [10..13] = data inicial (hh, mm, dd, yyyy)   <- mes ANTES do dia
         [14]=suporte (TAPE). A ordem dd/mm difere entre os dois grupos (quirk
         do legado); as datas do .fct sao a referencia (ver FctHeader).
linha 3: CPTEC AGCM REVIS 1.0 2000  T012642   COLD
depois : K valores de DeltaSigma, camada 1 = FUNDO, soma = 1
depois : tabela de variaveis, uma por linha,
         NOME (40 colunas) TIPO(FIXD|PROG|DIAG|SEQU) npts nlevels codigo

.fct (binario)
--------------
Fortran UNFORMATTED SEQUENCIAL, big endian, marcadores de 4 bytes:
  registro 1: 10 int32 (40 B) = ifday, nsec, idate(4), idatec(4), com
              idate/idatec = (hh, mm, dd, yyyy);
  em seguida, na ordem da tabela do .dir e, dentro de cada variavel, nos
  niveis k=1..nlev, um registro de npts float32:
    npts = (M+1)(M+2)  -> coeficientes espectrais na ordem do IC (m externo,
                          n interno, par cos/sen);
    npts = nlon*nlat   -> grade: i (longitude) mais rapido, linha 1 = polo
                          norte (verificado contra T2m/q2m e contra a sintese
                          espectral).
  Nao ha indice nem cabecalho por variavel: o deslocamento de cada registro sai
  da tabela + dos tamanhos. O campo nsec e' escrito como INTEIRO (segundos),
  mas o ReadHeader do PostGrib o le como REAL4 (FileAccess.f90:456) - aqui
  lemos o int32 e expomos tambem a leitura "real" por compatibilidade.

Uso tipico
----------
    info = parse_dir("GFCT...F.dir.TQ0126L042")
    fct = read_fct("GFCT...F.fct.TQ0126L042")     # .dir inferido do nome
    campos = fct.fields()                         # ps/t/q/topo + div/vort
    write_grads(...) / fct_to_grads(...)          # .ctl equivalente
"""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .derivatives import _fourier_to_grid as fourier_to_grid
from .derivatives import dpln_columns, wind_from_div_vort
from .ic_format import n_coeffs, pair_index
from .levels import gaussian_latitudes, grid_shape
from .spectral import legendre_columns
from .vertical import new_sigma

__all__ = [
    "FctError", "DirVariable", "DirHeader", "DirInfo", "FctHeader", "FctFile",
    "ALIASES", "parse_dir", "read_fct", "dir_for", "fct_files", "read_fct_series",
    "write_grads", "fct_to_grads",
]

HEADER_INTS = 10                 # ifday, nsec, idate(4), idatec(4)
HEADER_BYTES = 4 * HEADER_INTS   # 40
FCT_PATTERN = "*.fct*"
GRADS_UNDEF = 1e20

#: chave normalizada -> nomes aceitos na tabela do .dir
ALIASES: Dict[str, Tuple[str, ...]] = {
    "PS": ("LN SURFACE PRESSURE", "LOG SURFACE PRESSURE", "SURFACE PRESSURE"),
    "T": ("VIRTUAL TEMPERATURE", "TEMPERATURE"),
    "Q": ("SPECIFIC HUMIDITY",),
    "DIV": ("DIVERGENCE",),
    "VORT": ("VORTICITY",),
    "TOPO": ("TOPOGRAPHY",),
    "TSFC": ("SURFACE TEMPERATURE", "GROUND/SURFACE COVER TEMPERATURE"),
    "T2M": ("TEMPERATURE AT 2-M FROM SURFACE",),
    "Q2M": ("SPECIFIC HUMIDITY AT 2-M FROM SURFACE",),
    "LSM": ("LAND SEA MASK",),
    "O3": ("OZONE MIXING RATIO",),
}

#: chaves que o leitor sabe montar em grade a partir do espectro
SPECTRAL_KEYS = ("TOPO", "PS", "T", "Q", "DIV", "VORT")

#: nome da variavel no .gr/.ctl para cada chave
GRADS_NAMES = {"PS": "pslc", "T": "tvir", "Q": "umes", "U": "uvel", "V": "vvel",
               "TOPO": "topo"}
GRADS_UNITS = {"PS": "hPa", "T": "K", "Q": "kg/kg", "U": "m/s", "V": "m/s",
               "TOPO": "m"}

_VAR_RE = re.compile(
    r"^(?P<name>\S.*?)\s+(?P<kind>FIXD|PROG|DIAG|SEQU)\s+"
    r"(?P<npts>\d+)\s+(?P<nlev>\d+)\s+(?P<code>-?\d+)\s*$")

_MONTHS = ("jan", "feb", "mar", "apr", "may", "jun",
           "jul", "aug", "sep", "oct", "nov", "dec")


class FctError(ValueError):
    """Arquivo .dir/.fct invalido ou incompativel."""


# ---------------------------------------------------------------------------
# .dir
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DirVariable:
    """Uma linha da tabela de variaveis do .dir."""

    name: str
    kind: str          # FIXD | PROG | DIAG | SEQU
    npts: int          # palavras por registro (coeficientes ou pontos de grade)
    nlev: int          # registros (niveis) no .fct
    code: int          # codigo WMO-ish do legado (informativo)
    index: int         # ordem no arquivo (0-based)

    @property
    def nbytes(self) -> int:
        return 4 * self.npts


@dataclass(frozen=True)
class DirHeader:
    """Segunda linha do .dir (a primeira e a terceira sao fixas/informativas)."""

    file_version: str      # '0001'
    layout: str            # 'SEQU'
    trunc: str             # 'T0126'
    mwv: int               # 127 = Mend+1
    nz: int                # 42
    nz_layers: int         # 42 (repetido no arquivo)
    date_current: Tuple[int, int, int, int]   # (hh, dd, mm, yyyy)
    date_initial: Tuple[int, int, int, int]   # (hh, mm, dd, yyyy)
    support: str           # 'TAPE'

    @property
    def mend(self) -> int:
        return int(self.mwv) - 1

    @property
    def current_time(self) -> datetime:
        hh, dd, mm, yyyy = self.date_current
        return datetime(int(yyyy), int(mm), int(dd), int(hh))

    @property
    def initial_time(self) -> datetime:
        hh, mm, dd, yyyy = self.date_initial
        return datetime(int(yyyy), int(mm), int(dd), int(hh))


@dataclass(frozen=True)
class DirInfo:
    """Conteudo do .dir: cabecalho, sigma e tabela de variaveis."""

    path: Optional[Path]
    version: str                        # 'CPTEC SIGMA VERS 2.0'
    header: DirHeader
    revision: str                       # 'CPTEC AGCM REVIS 1.0 ...'
    sigma: np.ndarray                   # (K,) DeltaSigma, camada 1 = fundo
    sigma_interfaces: np.ndarray        # (K+1,) 1 -> 0 (fundo -> topo)
    sigma_layers: np.ndarray            # (K,) centros (formula de Phillips)
    variables: Tuple[DirVariable, ...]

    @property
    def mend(self) -> int:
        return self.header.mend

    @property
    def ncoeff(self) -> int:
        return n_coeffs(self.mend)

    @property
    def kmax(self) -> int:
        return len(self.sigma)

    @property
    def grid_points(self) -> int:
        nlon, nlat = grid_shape(self.mend)
        return nlon * nlat

    @property
    def names(self) -> Tuple[str, ...]:
        return tuple(v.name for v in self.variables)

    def variable(self, name: str) -> DirVariable:
        key = str(name).strip().upper()
        for var in self.variables:
            if var.name.upper() == key:
                return var
        for cand in ALIASES.get(key, ()):
            for var in self.variables:
                if var.name.upper() == cand:
                    return var
        raise KeyError("variavel %r ausente no .dir; disponiveis: %s"
                       % (name, list(self.names)))

    def to_dict(self) -> dict:
        return dict(version=self.version, header=self.header, sigma=self.sigma,
                    sigma_interfaces=self.sigma_interfaces,
                    sigma_layers=self.sigma_layers,
                    variables=list(self.variables), path=self.path)

    def __len__(self) -> int:
        return len(self.variables)


def _split_header_line(line: str, path: Path) -> DirHeader:
    tokens = line.split()
    if len(tokens) < 15:
        raise FctError("%s: linha 2 com %d campos (esperado >= 15): %r"
                       % (path, len(tokens), line))
    try:
        nums = [int(t) for t in tokens[3:14]]
    except ValueError:
        raise FctError("%s: campos numericos invalidos na linha 2: %r"
                       % (path, line)) from None
    return DirHeader(file_version=tokens[0], layout=tokens[1], trunc=tokens[2],
                     mwv=nums[0], nz=nums[1], nz_layers=nums[2],
                     date_current=tuple(nums[3:7]), date_initial=tuple(nums[7:11]),
                     support=tokens[14])


def parse_dir(path: str | Path) -> DirInfo:
    """Le o .dir e devolve cabecalho, sigma e tabela de variaveis.

    Equivale a opnpos/InitFiles do PostGrib. Levanta FctError com mensagem
    explicita quando o arquivo nao casa com o formato esperado.
    """
    p = Path(path)
    text = p.read_text(encoding="latin-1")
    lines = text.splitlines()
    if len(lines) < 4:
        raise FctError("%s: .dir com %d linhas (esperado >= 4)" % (p, len(lines)))
    version = lines[0].strip()
    if "CPTEC" not in version.upper() or "SIGMA" not in version.upper():
        raise FctError("%s: primeira linha inesperada: %r" % (p, version))
    header = _split_header_line(lines[1], p)
    revision = lines[2].strip()

    variables: List[DirVariable] = []
    sigma_vals: List[float] = []
    for raw in lines[3:]:
        line = raw.rstrip()
        if not line.strip():
            continue
        match = _VAR_RE.match(line)
        if match:
            variables.append(DirVariable(
                name=match.group("name").strip(),
                kind=match.group("kind"),
                npts=int(match.group("npts")),
                nlev=int(match.group("nlev")),
                code=int(match.group("code")),
                index=len(variables)))
            continue
        if variables:
            raise FctError("%s: linha invalida depois da tabela: %r" % (p, line))
        try:
            sigma_vals.extend(float(t) for t in line.split())
        except ValueError:
            raise FctError("%s: valor de sigma invalido: %r" % (p, line)) from None
    if not variables:
        raise FctError("%s: tabela de variaveis vazia" % (p,))
    sigma = np.asarray(sigma_vals, dtype=np.float64)
    if sigma.size != header.nz:
        raise FctError("%s: %d valores de sigma para nz=%d"
                       % (p, sigma.size, header.nz))
    total = float(sigma.sum())
    if abs(total - 1.0) > 1e-4:
        raise FctError("%s: sigma nao soma 1 (soma = %.8f)" % (p, total))
    inter, layers = new_sigma(sigma)
    return DirInfo(path=p, version=version, header=header, revision=revision,
                   sigma=sigma, sigma_interfaces=np.asarray(inter),
                   sigma_layers=np.asarray(layers),
                   variables=tuple(variables))


# ---------------------------------------------------------------------------
# .fct
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FctHeader:
    """Primeiro registro do .fct (40 bytes, 10 int32)."""

    forecast_day: int
    seconds: int                                  # segundos do dia de previsao
    date_inp: Tuple[int, int, int, int]           # (hh, mm, dd, yyyy)
    date_cur: Tuple[int, int, int, int]           # (hh, mm, dd, yyyy)

    @property
    def time_of_day_real(self) -> float:
        """Leitura do legado (REAL4): 1.4e-41, o PostGrib le assim mesmo."""
        return struct.unpack(">f", struct.pack(">i", int(self.seconds)))[0]

    @property
    def analysis_time(self) -> datetime:
        hh, mm, dd, yyyy = self.date_inp
        return datetime(int(yyyy), int(mm), int(dd), int(hh))

    @property
    def valid_time(self) -> datetime:
        hh, mm, dd, yyyy = self.date_cur
        return datetime(int(yyyy), int(mm), int(dd), int(hh))

    @property
    def lead_hours(self) -> float:
        """Horas de previsao (valid_time - analysis_time)."""
        delta = self.valid_time - self.analysis_time
        return delta.total_seconds() / 3600.0


def dir_for(fct_path: str | Path) -> Path:
    """Caminho do .dir correspondente a um .fct (mesmo prefixo)."""
    p = Path(fct_path)
    name = p.name.replace(".fct.", ".dir.")
    cand = p.with_name(name)
    if cand.exists():
        return cand
    if p.parent.is_dir():
        hits = sorted(p.parent.glob("*.dir.*"))
        if len(hits) == 1:
            return hits[0]
    raise FileNotFoundError("nao encontrei o .dir de %s (tentei %s)" % (p, cand))


def _read_record_marker(fh) -> int:
    raw = fh.read(4)
    if len(raw) != 4:
        raise FctError("fim de arquivo inesperado ao ler o marcador de registro")
    marker = struct.unpack(">i", raw)[0]
    if marker < 0:
        raise FctError("marcador de registro negativo (%d): arquivo corrompido"
                       % marker)
    return marker


@dataclass
class FctFile:
    """Um arquivo .fct (um tempo de previsao) com acesso aleatorio aos campos."""

    path: Path
    dir: DirInfo
    header: FctHeader
    offsets: Dict[str, Tuple[int, int]]      # nome -> (offset do 1o payload, nlev)

    # -- construcao ---------------------------------------------------------
    @property
    def variables(self) -> Tuple[DirVariable, ...]:
        return self.dir.variables

    @property
    def names(self) -> Tuple[str, ...]:
        return self.dir.names

    def __len__(self) -> int:
        return len(self.variables)

    @property
    def mend(self) -> int:
        return self.dir.mend

    @property
    def nlon(self) -> int:
        return grid_shape(self.mend)[0]

    @property
    def nlat(self) -> int:
        return grid_shape(self.mend)[1]

    @property
    def lead_hours(self) -> float:
        return self.header.lead_hours

    # -- leitura ------------------------------------------------------------
    def _offset(self, name: str, level: int) -> Tuple[int, int]:
        var = self.dir.variable(name)
        base, nlev = self.offsets[var.name]
        if not (0 <= int(level) < nlev):
            raise IndexError("nivel %d fora de 1..%d para %r"
                             % (int(level) + 1, nlev, var.name))
        return base + int(level) * (var.nbytes + 8), var.npts

    def read(self, name: str, level: Optional[int] = None) -> np.ndarray:
        """Le uma variavel (niveis x npts) ou um nivel (npts), em float64."""
        var = self.dir.variable(name)
        npts = var.npts
        levels = range(var.nlev) if level is None else (int(level),)
        out = np.empty((var.nlev if level is None else 1, npts))
        for k, lev in enumerate(levels):
            off, _ = self._offset(var.name, lev)
            out[k] = np.fromfile(self.path, dtype=">f4", count=npts,
                                 offset=off + 4).astype(np.float64)
        return out if level is None else out[0]

    def spectral(self, name: str) -> np.ndarray:
        """Coeficientes espectrais (nlev, (M+1)(M+2)) de uma variavel."""
        var = self.dir.variable(name)
        if var.npts != self.dir.ncoeff:
            raise FctError("%r tem npts=%d; esperado %d coeficientes espectrais"
                           % (var.name, var.npts, self.dir.ncoeff))
        return self.read(var.name)

    def grid(self, name: str) -> np.ndarray:
        """Campo em grade (nlev, nlat, nlon) de uma variavel em grade."""
        var = self.dir.variable(name)
        if var.npts != self.nlon * self.nlat:
            raise FctError("%r tem npts=%d; esperado %d pontos de grade (%dx%d)"
                           % (var.name, var.npts, self.nlon * self.nlat,
                              self.nlat, self.nlon))
        return self.read(var.name).reshape(var.nlev, self.nlat, self.nlon)

    def is_spectral(self, name: str) -> bool:
        return self.dir.variable(name).npts == self.dir.ncoeff

    # -- campos fisicos -----------------------------------------------------
    def fields(self, *, keys: Sequence[str] = SPECTRAL_KEYS,
               clip_q: bool = True) -> dict:
        """Campos em grade a partir do espectro (+ div/vort complexos).

        Devolve dict com:
          ps (nlat,nlon) hPa; t (K), q (kg/kg) (nlev,nlat,nlon); topo (m);
          div_sh, vort_sh (nlev, m+1, n+1) complexos; lat, lon, nlat, nlon,
          mend, header, sigma_layers.
        'q' e' clipada em 0 (ringing espectral deixa valores ~-1e-6).
        """
        want = tuple(k.upper() for k in keys)
        mend, nlat, nlon = self.mend, self.nlat, self.nlon
        out: dict = dict(header=self.header, mend=mend, nlat=nlat, nlon=nlon,
                         lat=gaussian_latitudes(nlat),
                         lon=np.arange(nlon, dtype=np.float64) * 360.0 / nlon,
                         sigma=self.dir.sigma,
                         sigma_interfaces=self.dir.sigma_interfaces,
                         sigma_layers=self.dir.sigma_layers)
        cache: Dict[str, np.ndarray] = {}
        mu = np.sin(np.radians(gaussian_latitudes(nlat)))
        columns = tuple(legendre_columns(mend, mu))

        def coefs(key: str) -> np.ndarray:
            return _split_pairs(self.spectral(key), mend)

        if "TOPO" in want:
            out["topo"] = _synth_stack(coefs("TOPO"), mend, nlat, nlon, columns)[0]
        if "PS" in want:
            lnps = _synth_stack(coefs("PS"), mend, nlat, nlon, columns)[0]
            out["ps"] = 10.0 * np.exp(lnps)
            out["lnps"] = lnps
        for key, target in (("T", "t"), ("Q", "q")):
            if key in want:
                cache[target] = _synth_stack(coefs(key), mend, nlat, nlon, columns)
        if "q" in cache and clip_q:
            np.maximum(cache["q"], 0.0, out=cache["q"])
        out.update({k: v for k, v in cache.items()})
        if "DIV" in want:
            out["div_sh"] = coefs("DIV")
        if "VORT" in want:
            out["vort_sh"] = coefs("VORT")
        if "div_sh" in out and "vort_sh" in out:
            dcolumns = tuple(dpln_columns(mend, mu))
            out["u"], out["v"] = wind_from_div_vort(
                out["div_sh"], out["vort_sh"], mend, nlat, nlon,
                columns=columns, dcolumns=dcolumns)
        return out

    # -- verificacao --------------------------------------------------------
    def verify(self) -> List[str]:
        """Checagens objetivas de estrutura (lista vazia = ok)."""
        problems: List[str] = []
        expected = fct_size(self.dir)
        size = self.path.stat().st_size
        if size != expected:
            problems.append("tamanho %d != esperado %d (diferenca %d bytes)"
                            % (size, expected, size - expected))
        try:
            with self.path.open("rb") as fh:
                marker = _read_record_marker(fh)
                if marker != HEADER_BYTES:
                    problems.append("marcador do cabecalho = %d (esperado %d)"
                                    % (marker, HEADER_BYTES))
                fh.seek(marker + 4, 1)      # payload do cabecalho + marcador final
                for var in self.variables:
                    for k in range(var.nlev):
                        pos = fh.tell()
                        marker = _read_record_marker(fh)
                        if marker != var.nbytes:
                            problems.append("%s k=%d: marcador %d != %d (offset %d)"
                                            % (var.name, k + 1, marker, var.nbytes, pos))
                        fh.seek(var.nbytes, 1)
                        tail = _read_record_marker(fh)
                        if tail != marker:
                            problems.append("%s k=%d: marcador final %d != inicial %d"
                                            % (var.name, k + 1, tail, marker))
        except (FctError, OSError) as exc:
            problems.append("leitura interrompida: %s" % (exc,))
        if abs(self.header.lead_hours - self.header.seconds / 3600.0) > 1e-6:
            problems.append("horas: datas dao %g h mas nsec/3600 = %g h"
                            % (self.header.lead_hours, self.header.seconds / 3600.0))
        return problems

    def summary(self) -> str:
        return ("%s: %s %02dz%02d%s%04d, %d variaveis, M=%d, K=%d, %.0f h"
                % (self.path.name, self.dir.header.trunc,
                   self.header.date_cur[0], self.header.date_cur[2],
                   _MONTHS[int(self.header.date_cur[1]) - 1],
                   self.header.date_cur[3], len(self.variables), self.mend,
                   self.dir.kmax, self.lead_hours))


def _split_pairs(vectors, mend: int) -> np.ndarray:
    """(nlev, ncoeff) na ordem do arquivo -> (nlev, m+1, n+1) complexo."""
    v = np.atleast_2d(np.asarray(vectors, dtype=np.float64))
    idx = pair_index(mend)
    mask = idx >= 0
    m_pos, n_pos = np.nonzero(mask)
    pair = idx[mask]
    coef = np.zeros((v.shape[0], mend + 1, mend + 1), dtype=np.complex128)
    coef[:, m_pos, n_pos] = v[:, 2 * pair] + 1j * v[:, 2 * pair + 1]
    return coef


def _synth_stack(coefs, mend: int, nlat: int, nlon: int, columns=None) -> np.ndarray:
    """Sintese vetorizada de todos os niveis de uma vez -> (nlev, nlat, nlon).

    Mesma convencao de bamci.spectral.synthesis (Legendre normalizado, fator 2
    para m>0), mas com uma unica passada de legendre_columns por campo.
    'columns' reaproveita as colunas de legendre_columns entre campos.
    """
    c = np.asarray(coefs, dtype=np.complex128)
    if c.ndim == 2:
        c = c[None]
    lat = gaussian_latitudes(nlat)
    mu = np.sin(np.radians(lat))
    if columns is None:
        columns = tuple(legendre_columns(mend, mu))
    y = np.zeros((c.shape[0], nlat, mend + 1), dtype=np.complex128)
    for m, pln_m in columns:
        y[:, :, m] = c[:, m, m:] @ pln_m              # (nlev, nlat) complexo
    return fourier_to_grid(y, nlon)


def fct_size(info: DirInfo) -> int:
    """Tamanho esperado do .fct (bytes) para o .dir dado."""
    total = 8 + HEADER_BYTES
    for var in info.variables:
        total += var.nlev * (8 + var.nbytes)
    return total


def read_fct(fct_path: str | Path, dir_path: Optional[str | Path] = None,
             *, check: bool = True) -> FctFile:
    """Abre um .fct (walk dos marcadores) e devolve um FctFile.

    dir_path=None usa o .dir de mesmo prefixo (dir_for). check=True valida o
    marcador do cabecalho e o tamanho total do arquivo.
    """
    p = Path(fct_path)
    info = parse_dir(dir_path if dir_path is not None else dir_for(p))
    size = p.stat().st_size
    with p.open("rb") as fh:
        marker = _read_record_marker(fh)
        if check and marker != HEADER_BYTES:
            raise FctError("%s: marcador do cabecalho = %d (esperado %d)"
                           % (p, marker, HEADER_BYTES))
        head_raw = fh.read(marker)
        if len(head_raw) != marker:
            raise FctError("%s: cabecalho truncado" % (p,))
    words = np.frombuffer(head_raw[:HEADER_BYTES], dtype=">i4")
    if words.size != HEADER_INTS:
        raise FctError("%s: cabecalho com %d palavras (esperado %d)"
                       % (p, words.size, HEADER_INTS))
    header = FctHeader(forecast_day=int(words[0]), seconds=int(words[1]),
                       date_inp=tuple(int(v) for v in words[2:6]),
                       date_cur=tuple(int(v) for v in words[6:10]))

    offsets: Dict[str, Tuple[int, int]] = {}
    with p.open("rb") as fh:
        fh.seek(8 + marker)                    # pula cabecalho + marcador final
        for var in info.variables:
            pos = fh.tell()
            marker = _read_record_marker(fh)
            if check and marker != var.nbytes:
                raise FctError("%s: %s: marcador %d != %d (offset %d)"
                               % (p, var.name, marker, var.nbytes, pos))
            offsets[var.name] = (pos, var.nlev)
            fh.seek(var.nlev * (var.nbytes + 8) - 4, 1)
    if check and size != fct_size(info):
        raise FctError("%s: tamanho %d != esperado %d (use o .dir correto)"
                       % (p, size, fct_size(info)))
    return FctFile(path=p, dir=info, header=header, offsets=offsets)


def fct_files(path: str | Path, pattern: str = FCT_PATTERN) -> List[Path]:
    """Lista .fct de um diretorio/glob/arquivo, ordenado pela hora de previsao."""
    p = Path(path)
    if p.is_file():
        return [p]
    if p.is_dir():
        found = sorted(p.glob(pattern))
    elif any(ch in str(path) for ch in "*?["):
        found = sorted(Path().glob(str(path)))
    else:
        raise FileNotFoundError("nao existe: %s" % (path,))
    found = [f for f in found if f.is_file() and ".fct" in f.name]
    if not found:
        raise FileNotFoundError("nenhum %r em %s" % (pattern, path))
    def lead(f: Path) -> float:
        try:
            return read_fct(f, check=False).lead_hours
        except (FctError, OSError):
            return float("inf")
    return sorted(found, key=lead)


def read_fct_series(path: str | Path, pattern: str = FCT_PATTERN) -> List[FctFile]:
    """Abre todos os .fct de um diretorio/glob, em ordem de hora de previsao."""
    return [read_fct(f) for f in fct_files(path, pattern)]


# ---------------------------------------------------------------------------
# .ctl equivalente (GrADS) para reuso pelo resto do pacote
# ---------------------------------------------------------------------------
def write_grads(fields: Mapping[str, np.ndarray], prefix: str | Path, *,
                lat, lon, levels=None, time: Optional[datetime] = None,
                title: str = "", units: Optional[Mapping[str, str]] = None,
                level_fmt: str = "%.6f") -> Tuple[Path, Path]:
    """Escreve um par .dat + .ctl no layout que o GrADS (e bamci.grads) le.

    fields: nome -> (nlat,nlon) ou (nlev,nlat,nlon), na ordem em que devem
    aparecer no .ctl/.dat. lat/lon em graus (lat norte -> sul, como no modelo).
    O binario e' float32 big endian, linhas norte -> sul, x mais rapido, com
    'options yrev big_endian' e ydef de sul para norte (igual ao PostGrib).
    """
    base = Path(prefix)
    dat = base.with_suffix(".dat")
    ctl = base.with_suffix(".ctl")
    lat = np.asarray(lat, dtype=np.float64).ravel()
    lon = np.asarray(lon, dtype=np.float64).ravel()
    nlat, nlon = lat.size, lon.size
    entries: List[Tuple[str, np.ndarray]] = []
    for key, value in fields.items():
        arr = np.asarray(value, dtype=np.float64)
        if arr.ndim == 2:
            arr = arr[None]
        if arr.ndim != 3 or arr.shape[1:] != (nlat, nlon):
            raise ValueError("campo %r com forma %s; esperado (nlev,%d,%d)"
                             % (key, np.shape(value), nlat, nlon))
        entries.append((str(key), arr))
    if not entries:
        raise ValueError("nenhum campo para escrever")
    nlev = max(arr.shape[0] for _, arr in entries)
    if levels is None:
        levels = np.arange(1, nlev + 1, dtype=np.float64)
    levels = np.asarray(levels, dtype=np.float64).ravel()
    if levels.size != nlev:
        raise ValueError("%d niveis declarados para %d camadas"
                         % (levels.size, nlev))
    with dat.open("wb") as fh:
        for _, arr in entries:
            for k in range(arr.shape[0]):
                fh.write(np.ascontiguousarray(arr[k], dtype=">f4").tobytes())
    when = time if time is not None else datetime(2000, 1, 1)
    stamp = "%02dz%02d%s%04d" % (when.hour, when.day, _MONTHS[when.month - 1],
                                when.year)
    units = dict(units or {})
    lines = ["dset ^%s" % dat.name,
             "title %s" % (title or base.name),
             "options yrev big_endian",
             "undef %g" % GRADS_UNDEF,
             "xdef %d linear 0.0 %g" % (nlon, 360.0 / nlon),
             # ydef/zdef numa unica linha: e' o que bamci.ctl._parse_axis le
             # (o GrADS tambem aceita a lista numa linha so')
             "ydef %d levels %s" % (nlat, " ".join("%.3f" % float(v)
                                                   for v in lat[::-1])),
             "zdef %d levels %s" % (nlev, " ".join(level_fmt % float(v)
                                                   for v in levels))]
    lines += ["tdef 1 linear %s 1hr" % stamp,
              "vars %d" % len(entries)]
    for key, arr in entries:
        name = GRADS_NAMES.get(key.upper(), key.lower())
        unit = units.get(key, GRADS_UNITS.get(key.upper(), "")) or "-"
        desc = key.upper()
        if arr.shape[0] == 1:
            # 2-D: mesma forma do PostGrib ('0 99' = um nivel), que o parser
            # de bamci.ctl e o GrADS entendem
            lines.append("%-8s 0 99 %s %s" % (name, unit, desc))
        else:
            lev_txt = " ".join(level_fmt % float(v) for v in levels[:arr.shape[0]])
            lines.append("%-8s %d %s %s %s"
                         % (name, arr.shape[0], lev_txt, unit, desc))
    lines.append("endvars")
    ctl.write_text("\n".join(lines) + "\n", encoding="ascii")
    return dat, ctl


def fct_to_grads(fct: FctFile, prefix: str | Path, *, keys=("PS", "T", "Q", "U", "V"),
                 title: str = "") -> Tuple[Path, Path]:
    """Grava .dat + .ctl de um .fct com ps/t/q/u/v nos niveis sigma do modelo."""
    data = fct.fields(keys=("PS", "T", "Q", "DIV", "VORT"))
    fields: Dict[str, np.ndarray] = {}
    for key in keys:
        up = str(key).upper()
        if up == "PS":
            fields["PS"] = data["ps"]
        elif up in ("T", "Q", "U", "V"):
            fields[up] = data[up.lower()]
        else:
            raise KeyError("chave %r nao suportada (PS/T/Q/U/V)" % (key,))
    return write_grads(fields, prefix, lat=data["lat"], lon=data["lon"],
                       levels=data["sigma_layers"], time=fct.header.valid_time,
                       title=title or fct.summary())
