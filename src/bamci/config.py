"""Configuracao da geracao de condicoes iniciais do BAM.

Reproduz o namelist ICnAGCMNML do Fortran legado (Sources/ICnAGCM/InputParameters.f90)
com os presets por fonte de dados, estendidos com ERA5.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple


#: parametros por tipo de dado de entrada (ImaxInp, JmaxInp, KmaxInp, KmaxInq, DataPrs)
CASE_PRESETS: Dict[str, Dict[str, object]] = {
    "E40": dict(imax=320, jmax=161, kmax=23, kmaxq=23, levels="PressureLevels.ERA40C", grid_deg=1.125),
    "E4W": dict(imax=144, jmax=73, kmax=23, kmaxq=23, levels="PressureLevels.ERA40W", grid_deg=2.5),
    "EIT": dict(imax=240, jmax=121, kmax=37, kmaxq=37, levels="PressureLevels.ERAINT", grid_deg=1.5),
    "EIH": dict(imax=480, jmax=241, kmax=37, kmaxq=37, levels="PressureLevels.ERAINH", grid_deg=0.75),
    "FNL": dict(imax=360, jmax=181, kmax=26, kmaxq=21, levels="PressureLevels.NMCFNL", grid_deg=1.0),
    "GFS": dict(imax=1440, jmax=721, kmax=26, kmaxq=21, levels="PressureLevels.NMCGFS", grid_deg=0.25),
    "ERA5": dict(imax=1440, jmax=721, kmax=37, kmaxq=37, levels="PressureLevels.ERA5", grid_deg=0.25),
}


@dataclass
class Config:
    """Todos os parametros de uma execucao de pre-processamento."""

    case: str = "ERA5"
    date: str = "2026100100"          # YYYYMMDDHH
    mend_out: int = 126               # truncacao espectral de saida (T126 para o subsazonal)
    kmax_out: int = 42                # numero de camadas do modelo
    grid: str = "Q"                   # Q = quadratica, L = linear
    sigma_level: bool = True          # BAM 1.2 subsazonal le SIGMA (decidido em 2026-11-07)

    # grade de entrada (sobrescreve o preset quando informado)
    imax_in: int | None = None
    jmax_in: int | None = None
    kmax_in: int | None = None
    kmaxq_in: int | None = None

    # diretorios
    data_in: Path = Path(".")         # .../<CASE>/*.dat (equivalente a DataAnl)
    data_out: Path = Path(".")        # diretorio de saida
    levels_dir: Path = Path(".")      # PressureLevels.*, HybridLevels.Lxxx, DeltaSigma.Lxxx
    topo_dir: Path = Path(".")        # Topography.TQxxxx / Topography.Gxxxxx

    # arquivos de niveis
    levels_file: str | None = None
    hybrid_file: str | None = None
    sigma_file: str | None = None
    topo_file: str | None = None
    topo_grid_file: str | None = None

    # suavizacao da topografia
    iter_smooth: int = 10
    smooth_percent: float = 0.12
    smooth_fields: bool = True      # filtro de Hoskins nos coeficientes (producao usa)
    lon_shift: int = 0              # deslocamento inteiro em indices (medido: 0)
    lon_reverse: bool = True        # grade do modelo com longitudes decrescentes (medido contra a producao)

    # opcoes de saida
    grads: bool = False               # escreve tambem o .GrADS de diagnostico
    output_template: str = "GANL{case}{date}S.unf.T{grid}{mend:04d}L{kmax:03d}"
    sigma_suffix: bool = False        # legado: _sigma/_hibrid no nome

    # constantes fisicas (identicas ao Fortran)
    earth_radius: float = 6.37e6
    grav: float = 9.80665
    rd: float = 287.05
    rv: float = 461.50
    cp: float = 1004.6

    extra: Dict[str, object] = field(default_factory=dict)

    # -- derivados ---------------------------------------------------------
    def preset(self) -> Dict[str, object]:
        if self.case not in CASE_PRESETS:
            raise KeyError("case desconhecido: %r (use %s)" % (self.case, sorted(CASE_PRESETS)))
        return CASE_PRESETS[self.case]

    @property
    def nx(self) -> int:
        return int(self.imax_in or self.preset()["imax"])

    @property
    def ny(self) -> int:
        return int(self.jmax_in or self.preset()["jmax"])

    @property
    def nz(self) -> int:
        return int(self.kmax_in or self.preset()["kmax"])

    @property
    def nzq(self) -> int:
        return int(self.kmaxq_in or self.preset()["kmaxq"])

    @property
    def levels_name(self) -> str:
        return str(self.levels_file or self.preset()["levels"])

    @property
    def ctvr(self) -> float:
        """cTv = Rv/Rd - 1 (fator da temperatura virtual)."""
        return self.rv / self.rd - 1.0

    @property
    def rocp(self) -> float:
        return self.rd / self.cp

    def output_path(self) -> Path:
        name = self.output_template.format(
            case=self.case, date=self.date, grid=self.grid, mend=self.mend_out, kmax=self.kmax_out
        )
        if self.sigma_suffix:
            name += "_sigma" if self.sigma_level else "_hibrid"
        return Path(self.data_out) / name

    def check(self) -> List[str]:
        problems: List[str] = []
        if len(self.date) != 10 or not self.date.isdigit():
            problems.append("date deve ser YYYYMMDDHH")
        if self.grid not in ("Q", "L"):
            problems.append("grid deve ser Q ou L")
        if self.kmax_out <= 0 or self.mend_out <= 0:
            problems.append("mend_out/kmax_out devem ser positivos")
        return problems
