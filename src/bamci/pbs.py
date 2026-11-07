"""Submissao dos jobs na jaci via PBS (fila aux).

O pre-processamento roda como um unico programa Python, mas na jaci o padrao da
casa e submeter por PBS. Este modulo gera o script de submissao e (opcionalmente)
chama o qsub.

Convencoes medidas na jaci (2026-11):
  * PBS Pro/OpenPBS em /opt/pbs/bin (qsub, qstat, qdel);
  * fila aux habilitada, walltime maximo 24:00:00;
  * nos de 512 cpus (cn-0001..), memoria ~756 GB.
"""

from __future__ import annotations

import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

QUEUE = "aux"
MAX_WALLTIME = 24
DEFAULT_ENV = ["module load anaconda"]
#: interpretador que existe nos nos de calculo da jaci (o 'python' do PATH nao existe la)
DEFAULT_PYTHON = "/p/projetos/ioper/models/BAM-1.2-subsaz/bam-ci/.venv/bin/python"


@dataclass
class PBSSpec:
    """Parametros de um job PBS."""

    job_name: str
    command: Sequence[str]
    queue: str = QUEUE
    ncpus: int = 1
    walltime: str = "01:00:00"
    workdir: str = "."
    log_name: Optional[str] = None
    env_setup: List[str] = field(default_factory=lambda: list(DEFAULT_ENV))
    export_env: bool = True

    def log_path(self) -> str:
        if self.log_name:
            return self.log_name
        return "%s.log" % self.job_name

    def validate(self) -> List[str]:
        problems: List[str] = []
        if self.ncpus < 1:
            problems.append("ncpus deve ser >= 1")
        hours = int(str(self.walltime).split(":")[0])
        if hours > MAX_WALLTIME:
            problems.append("walltime acima do maximo da fila %s (%dh)" % (self.queue, MAX_WALLTIME))
        if not self.command:
            problems.append("comando vazio")
        return problems


def render(spec: PBSSpec) -> str:
    """Devolve o texto do script PBS."""
    problems = spec.validate()
    if problems:
        raise ValueError("; ".join(problems))
    cmd = " ".join(shlex.quote(str(c)) for c in spec.command)
    lines = [
        "#!/bin/bash",
        "#PBS -N %s" % spec.job_name,
        "#PBS -q %s" % spec.queue,
        "#PBS -l select=1:ncpus=%d" % spec.ncpus,
        "#PBS -l walltime=%s" % spec.walltime,
        "#PBS -o %s" % spec.log_path(),
        "#PBS -j oe",
    ]
    if spec.export_env:
        lines.append("#PBS -V")
    if str(spec.workdir) in (".", ""):
        # o PBS inicia o job em $HOME; o diretorio de submissao vem em $PBS_O_WORKDIR
        cd_line = 'cd "$PBS_O_WORKDIR" || exit 1'
    else:
        cd_line = "cd %s || exit 1" % shlex.quote(str(spec.workdir))
    lines += ["", cd_line, ""]
    lines += [str(e) for e in spec.env_setup]
    lines += ["", cmd, ""]
    return "\n".join(lines)


def write_script(spec: PBSSpec, path: str | Path) -> Path:
    p = Path(path)
    p.write_text(render(spec), encoding="utf-8")
    p.chmod(0o755)
    return p


def submit(path: str | Path, dry_run: bool = False) -> str:
    """Submete o script com qsub; dry_run devolve o comando sem executar."""
    cmd = ["qsub", str(path)]
    if dry_run:
        return " ".join(cmd)
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("qsub falhou (%d): %s" % (out.returncode, out.stderr.strip()))
    return out.stdout.strip()
