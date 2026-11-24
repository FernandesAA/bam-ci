"""Driver do modelo BAM 1.2 subsazonal e do ensaio MB09.

Este modulo faz a ponte entre o pre-processamento (bamci run / bamci sst), o
nucleo MB09 (bamci.mb09) e a execucao do modelo na JACI. Ele NAO reimplementa o
modelo nem o MB09: apenas

* renderiza o namelist MODELIN (&MODEL_RES) preservando o restante do arquivo;
* renderiza o script PBS do modelo (aprun + OpenMP);
* prepara o diretorio de execucao e os arquivos de entrada que o modelo le
  (GANLNMC<data>S.unf.<res> e SSTMonthlyDirec<data>.G00<lat>);
* submete/espera jobs PBS;
* coleta as saidas e roda o PostGrib.

Convencoes medidas na JACI (2026-11) estao nos defaults das constantes abaixo.
As funcoes sao puras ou aceitam dry_run=True para serem testadas sem a JACI e
sem disparar subprocess.
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable, Mapping, Sequence

#: suite operacional do BAM 1.2 subsazonal na JACI
DEFAULT_SUITE = "/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper"
#: executavel do modelo (MPI + OpenMP, namelist na stdin)
DEFAULT_MODEL_EXE = (
    "/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper/model/exec/ParModel_MPI"
)
#: carrega intel/2022.0.2, mpi/2021.5.1 e NETCDF/PNETCDF/PIO
DEFAULT_ENV_SOURCE = (
    "/p/projetos/ioper/fix/sazonal/load_modulei_intel_2022_oneapi_mpi.bash"
)
#: qsub do PBS Pro na JACI (o do PATH nem sempre e o certo)
DEFAULT_QSUB = "/opt/pbs/default/bin/qsub"
DEFAULT_QSTAT = "/opt/pbs/bin/qstat"

DEFAULT_QUEUE = "workq"
DEFAULT_ACCOUNT = "CPTEC"
DEFAULT_WALLTIME = "23:00:00"
DEFAULT_NODES = 3
DEFAULT_PPN = 40
DEFAULT_OMP_THREADS = 10
#: ranks MPI usados para TQ0126L042 (12 ranks, 4 por no, 10 threads cada)
DEFAULT_RANKS = 12
DEFAULT_RANKS_PER_NODE = 4
#: resolucao alvo do ensaio MB09
DEFAULT_RESOLUTION = "TQ0126L042"
MODELIN_NAME = "MODELIN"
POSTGRIB_LOG = "postgrib.log"
OMP_STACKSIZE = "256M"

#: nomes que o modelo efetivamente abre em <suite>/model/datain
IC_TEMPLATE = "GANLNMC{date}S.unf.{resolution}"
#: nomenclatura opcional do ensaio operacional (CTR, 01P, 01N, ...)
IC_OPER_TEMPLATE = "GANL{tag}{date}S.unf.{resolution}"
SST_TEMPLATE = "SSTMonthlyDirec{date}.G{latitude:05d}"

#: chaves cujo valor textual precisa de aspas no namelist
_STRING_KEYS = frozenset(
    {
        "path_in",
        "dirfnameoutput",
        "nmsst",
        "prefx",
        "prefy",
        "table",
        "iccon",
        "iscon",
        "ilcon",
        "igwd",
        "ocflux",
        "slabocean",
        "icemodel",
    }
)

#: chaves de data no formato hh,dd,mm,yyyy
_DATE_KEYS = frozenset({"idatei", "idatef", "idatew"})


# ---------------------------------------------------------------------------
# datas e nomes
# ---------------------------------------------------------------------------
def parse_date(value: str | datetime) -> str:
    """Normaliza YYYYMMDDHH (aceita YYYY-MM-DDTHH, YYYYMMDD etc.)."""
    if isinstance(value, datetime):
        return value.strftime("%Y%m%d%H")
    raw = str(value).strip()
    if raw.endswith("Z"):
        raw = raw[:-1]
    digits = re.sub(r"\D", "", raw)
    if len(digits) == 8:
        digits += "00"
    elif len(digits) == 12:  # minutos/segundos: a hora basta
        digits = digits[:10]
    if len(digits) != 10:
        raise ValueError("data invalida: %r (use YYYYMMDDHH)" % (value,))
    try:
        datetime.strptime(digits, "%Y%m%d%H")
    except ValueError:
        raise ValueError("data inexistente: %r" % (value,)) from None
    return digits


def format_idate(value: str | datetime) -> str:
    """Converte 2026110500 em 00,05,11,2026 (ordem hh,dd,mm,yyyy)."""
    moment = datetime.strptime(parse_date(value), "%Y%m%d%H")
    return "%02d,%02d,%02d,%04d" % (
        moment.hour,
        moment.day,
        moment.month,
        moment.year,
    )


def shift_date(value: str | datetime, hours: int) -> str:
    """Soma horas a uma data YYYYMMDDHH e devolve o mesmo formato."""
    moment = datetime.strptime(parse_date(value), "%Y%m%d%H")
    return (moment + timedelta(hours=int(hours))).strftime("%Y%m%d%H")


def run_id(labeli: str | datetime, labelf: str | datetime) -> str:
    """Rotulo do diretorio de saida: <LABELI>-<LABELF>."""
    return "%s-%s" % (parse_date(labeli), parse_date(labelf))


def resolution_from_path(path: str | Path) -> str | None:
    """Extrai TQ0126L042 do sufixo ....unf.TQ0126L042."""
    name = Path(path).name
    match = re.search(r"\.unf\.([A-Za-z][A-Za-z0-9_]*)$", name)
    return match.group(1) if match else None


def nlat_from_path(path: str | Path) -> int | None:
    """Extrai 192 de ...G00192."""
    match = re.search(r"\.G(\d{5})(?:\.[^.]+)?$", Path(path).name)
    return int(match.group(1)) if match else None


# ---------------------------------------------------------------------------
# namelist MODELIN
# ---------------------------------------------------------------------------
def _modelin_value(key: str, value) -> str:
    """Formata um valor de override conforme o tipo da chave."""
    low = key.lower()
    if low in _DATE_KEYS:
        text = str(value)
        if isinstance(value, datetime) or re.fullmatch(r"\d{10}", text.strip()):
            return format_idate(value)
        return text
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if low in _STRING_KEYS:
        text = str(value)
        if text[:1] not in ("'", '"'):
            text = "'" + text.replace("'", "''") + "'"
        return text
    return str(value)


#: um valor escalar do namelist vai ate a virgula, o comentario ou o fim da linha
_SCALAR_VALUE = r"[^,!\n]+"
#: string Fortran (aspas simples ou duplas, com aspas duplicadas internas)
_QUOTED_VALUE = r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\""
#: data hh,dd,mm,yyyy (as virgulas internas nao terminam o valor)
_DATE_VALUE = r"[-+]?\d+\s*,\s*[-+]?\d+\s*,\s*[-+]?\d+\s*,\s*[-+]?\d+"


def _commented(text: str, position: int) -> bool:
    """True se position esta depois de um '!' de comentario na mesma linha."""
    line_start = text.rfind("\n", 0, position) + 1
    return "!" in text[line_start:position]


def _namelist_pattern(key: str, value_pattern: str) -> re.Pattern:
    return re.compile(
        r"(?i)(?<![A-Za-z0-9_])(?P<pre>" + re.escape(key) + r"[ \t]*=[ \t]*)"
        r"(?P<val>" + value_pattern + r")(?P<sep>[ \t]*,)?"
    )


def _replace_namelist_key(text: str, key: str, value) -> str:
    """Troca so o valor de key, preservando virgula, espacos e comentario.

    Aceita varias chaves por linha (o MODELIN real agrupa algumas), datas com
    virgulas internas e chaves repetidas (erro). Chaves citadas em comentarios
    sao ignoradas.
    """
    low = key.lower()
    if low in _DATE_KEYS:
        candidates = [_DATE_VALUE]
    elif low in _STRING_KEYS:
        candidates = [_QUOTED_VALUE, _SCALAR_VALUE]
    else:
        candidates = [_SCALAR_VALUE]

    rendered = _modelin_value(key, value)
    for value_pattern in candidates:
        pattern = _namelist_pattern(key, value_pattern)
        matches = [m for m in pattern.finditer(text) if not _commented(text, m.start())]
        if not matches:
            continue
        if len(matches) > 1:
            raise ValueError("chave %r aparece %d vezes no MODELIN" % (key, len(matches)))
        match = matches[0]
        replacement = match.group("pre") + rendered + (match.group("sep") or "")
        return text[: match.start()] + replacement + text[match.end():]
    raise KeyError("chave %r nao encontrada no MODELIN" % (key,))


def render_modelin_text(template_text: str, overrides: Mapping | None = None, **kw) -> str:
    """Renderiza o namelist em memoria (funcao pura, usada nos testes).

    Aceita date=... como atalho para idatei e idatef (sem sobrescrever esses
    dois quando passados explicitamente). As demais chaves sao repassadas como
    chave=valor (case-insensitive, como no Fortran).
    """
    values = dict(overrides or {})
    values.update(kw)
    if "date" in values:
        date = values.pop("date")
        values.setdefault("idatei", date)
        values.setdefault("idatef", date)
    text = template_text
    for key, value in values.items():
        if value is None:
            continue
        text = _replace_namelist_key(text, str(key), value)
    return text


def render_modelin(template_path: str | Path, out_path: str | Path, **overrides) -> Path:
    """Le o template, troca as chaves pedidas e grava out_path.

    Apenas as chaves passadas em overrides mudam; ordem, comentarios e demais
    linhas do namelist ficam intactos.
    """
    text = Path(template_path).read_text(encoding="utf-8")
    rendered = render_modelin_text(text, overrides)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(rendered, encoding="utf-8")
    return out


# ---------------------------------------------------------------------------
# script PBS do modelo
# ---------------------------------------------------------------------------
def _q(value) -> str:
    return shlex.quote(str(value))


def render_pbs(
    job_name: str,
    workdir: str | Path,
    exec_dir: str | Path,
    modelin: str | Path = MODELIN_NAME,
    aprun_args: str | Sequence[str] | None = None,
    queue: str = DEFAULT_QUEUE,
    account: str = DEFAULT_ACCOUNT,
    walltime: str = DEFAULT_WALLTIME,
    nodes: int | None = None,
    ppn: int | None = None,
    omp_threads: int = DEFAULT_OMP_THREADS,
    env_source: str | Path | Sequence[str] | None = DEFAULT_ENV_SOURCE,
    *,
    ranks: int = DEFAULT_RANKS,
    ranks_per_node: int = DEFAULT_RANKS_PER_NODE,
    pbs_style: str = "mpp",
    mppwidth: int | None = None,
    mppnppn: int | None = None,
    mppdepth: int | None = None,
    exe_name: str = "ParModel_MPI",
    stdout: str = "print.out",
    log_name: str | None = None,
    omp_stacksize: str = OMP_STACKSIZE,
) -> str:
    """Texto do script PBS que roda o modelo com aprun.

    pbs_style escolhe a sintaxe do recurso:
      * mpp:    -l mppwidth=N,mppnppn=N,mppdepth=N (Moab/Torque do oensMB09)
      * nodes:  -l nodes=N:ppn=N (bam_SeasonalClima/Torque classico)
      * select: -l select=N:ncpus=N:mpiprocs=N (PBS Pro)
    Em nodes/select o numero de chunks vem de nodes (ou de
    ceil(ranks/ranks_per_node)). aprun_args pode ser uma string pronta ou uma
    sequencia; quando omitido usa -m1000h -n <ranks> -N <ranks_per_node> -d
    <threads>, que precisa casar com o recurso pedido.
    """
    style = str(pbs_style).lower()
    if style not in ("mpp", "nodes", "select"):
        raise ValueError("pbs_style deve ser 'mpp', 'nodes' ou 'select'")
    ranks = int(ranks)
    ranks_per_node = int(ranks_per_node)
    omp_threads = int(omp_threads)
    if ranks < 1 or ranks_per_node < 1 or omp_threads < 1:
        raise ValueError("ranks, ranks_per_node e omp_threads devem ser >= 1")
    try:
        hours = int(str(walltime).split(":")[0])
    except (TypeError, ValueError):
        raise ValueError("walltime invalido: %r (use HH:MM:SS)" % (walltime,)) from None
    if hours > 24:
        raise ValueError("walltime acima do maximo de 24 h: %s" % (walltime,))

    if aprun_args is None:
        aprun = "-m1000h -n %d -N %d -d %d" % (ranks, ranks_per_node, omp_threads)
    elif isinstance(aprun_args, str):
        aprun = aprun_args
    else:
        aprun = " ".join(_q(arg) for arg in aprun_args)

    modelin_text = str(modelin) if modelin else MODELIN_NAME
    exec_path = Path(exec_dir) if exec_dir else Path(workdir)
    command = "aprun %s ./%s < %s > %s" % (
        aprun,
        _q(exe_name),
        _q(modelin_text),
        _q(stdout),
    )

    chunks = (
        int(nodes) if nodes is not None
        else (ranks + ranks_per_node - 1) // ranks_per_node
    )
    if int(chunks) < 1:
        raise ValueError("nodes/chunks deve ser >= 1")
    if style == "mpp":
        width = ranks if mppwidth is None else int(mppwidth)
        nppn = ranks_per_node if mppnppn is None else int(mppnppn)
        depth = omp_threads if mppdepth is None else int(mppdepth)
        if width < 1 or nppn < 1 or depth < 1:
            raise ValueError("mppwidth/mppnppn/mppdepth devem ser >= 1")
        resource = "#PBS -l mppwidth=%d,mppnppn=%d,mppdepth=%d" % (
            width, nppn, depth)
    elif style == "nodes":
        ppn_use = ranks_per_node if ppn is None else int(ppn)
        if ppn_use < 1:
            raise ValueError("ppn deve ser >= 1")
        resource = "#PBS -l nodes=%d:ppn=%d" % (int(chunks), ppn_use)
    else:  # select (PBS Pro)
        resource = "#PBS -l select=%d:ncpus=%d:mpiprocs=%d" % (
            int(chunks), omp_threads, ranks_per_node)

    log = log_name or "%s.log" % job_name
    lines = [
        "#!/bin/bash",
        "#PBS -N %s" % job_name,
        "#PBS -q %s" % queue,
        "#PBS -A %s" % account,
        resource,
        "#PBS -l walltime=%s" % walltime,
        "#PBS -j oe",
        "#PBS -V",
        "#PBS -o %s" % log,
        "",
        "cd %s || exit 1" % _q(exec_path),
        "",
    ]
    if env_source:
        sources = [env_source] if isinstance(env_source, (str, Path)) else list(env_source)
        for source in sources:
            lines.append("source %s" % _q(source))
    lines += [
        "ulimit -s unlimited",
        "export OMP_NUM_THREADS=%d" % omp_threads,
        "export OMP_STACKSIZE=%s" % omp_stacksize,
        "",
        command,
        "",
    ]
    return "\n".join(lines)


def write_pbs(path: str | Path, text: str) -> Path:
    """Grava o script PBS com permissao de execucao."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    out.chmod(0o755)
    return out


# ---------------------------------------------------------------------------
# diretorio de execucao e entradas do modelo
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RunDir:
    """Caminhos de uma execucao do modelo."""

    label: str
    run_dir: Path
    executable: Path
    modelin: Path


def _safe_label(label: str | Path) -> str:
    text = str(label).strip().replace("/", "_").replace(os.sep, "_")
    if text in ("", ".", ".."):
        raise ValueError("rotulo invalido: %r" % (label,))
    return text


def prepare_run_dir(
    work_dir: str | Path,
    label: str | Path,
    executable: str | Path,
    *,
    modelin_name: str = MODELIN_NAME,
    dry_run: bool = False,
) -> RunDir:
    """Cria exec_NMC.<label>/ e copia o executavel do modelo."""
    safe = _safe_label(label)
    run_dir = Path(work_dir) / ("exec_NMC.%s" % safe)
    exe = run_dir / Path(executable).name
    modelin = run_dir / modelin_name
    if not dry_run:
        source = Path(executable)
        if not source.is_file():
            raise FileNotFoundError("executavel do modelo nao encontrado: %s" % source)
        run_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, exe)
        exe.chmod(0o755)
    return RunDir(label=safe, run_dir=run_dir, executable=exe, modelin=modelin)


@dataclass(frozen=True)
class LinkedInputs:
    """Arquivos de entrada com os nomes que o modelo abre."""

    ic: Path
    sst: Path
    mode: str


def _place(src: str | Path, dst: str | Path, mode: str = "link") -> str:
    """Coloca src em dst (link simbolico, com fallback para copia)."""
    source = Path(src)
    target = Path(dst)
    if not source.is_file():
        raise FileNotFoundError("arquivo de entrada nao encontrado: %s" % source)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink():
        try:
            if Path(os.readlink(target)) == source.resolve():
                return "link"
        except OSError:
            pass
        target.unlink()
    elif target.exists():
        target.unlink()
    if str(mode).lower() == "copy":
        shutil.copy2(source, target)
        return "copy"
    try:
        os.symlink(source.resolve(), target)
        return "link"
    except OSError:
        shutil.copy2(source, target)
        return "copy"


def _ic_target_name(
    label: str,
    resolution: str,
    ic_naming: str = "simple",
    *,
    control: bool = False,
    member: int | None = None,
    sign: str | None = None,
) -> str:
    """Nome do IC no datain: GANLNMC... (simple) ou GANL{CTR|NNP|NNN}... (operacional)."""
    naming = str(ic_naming).lower()
    if naming not in ("simple", "operacional"):
        raise ValueError("ic_naming deve ser 'simple' ou 'operacional'")
    if naming == "simple":
        return IC_TEMPLATE.format(date=label, resolution=resolution)
    if control:
        tag = "CTR"
    else:
        if member is None or sign is None:
            raise ValueError(
                "ic_naming='operacional' exige member (1..N) e sign (+/-) para membros"
            )
        suffix = {"+": "P", "-": "N", "P": "P", "N": "N"}.get(str(sign).strip().upper())
        if suffix is None:
            raise ValueError("sign deve ser '+' ou '-' (ou P/N), recebido %r" % (sign,))
        index = int(member)
        if index < 1:
            raise ValueError("member deve ser >= 1")
        tag = "%02d%s" % (index, suffix)
    return IC_OPER_TEMPLATE.format(tag=tag, date=label, resolution=resolution)


def link_inputs(
    suite: str | Path,
    ic_path: str | Path,
    sst_path: str | Path,
    labeli: str | datetime,
    *,
    resolution: str | None = None,
    mode: str = "link",
    datain: str | Path | None = None,
    nlat: int | None = None,
    ic_naming: str = "simple",
    control: bool = False,
    member: int | None = None,
    sign: str | None = None,
    dry_run: bool = False,
) -> LinkedInputs:
    """Cria o IC e a SST com os nomes que o modelo abre.

    O modelo le GANLNMC<data>S.unf.<res> e SSTMonthlyDirec<data>.G00<lat>. Com
    ic_naming='operacional' o IC tambem recebe o nome do ensaio
    (GANLCTR..., GANL01P..., GANL01N...) e um alias GANLNMC aponta para ele,
    para o modelo continuar lendo o nome fixo. datain permite apontar para um
    diretorio por membro do ensemble.
    """
    label = parse_date(labeli)
    root = Path(datain) if datain is not None else Path(suite) / "model" / "datain"
    ic_source = Path(ic_path)
    sst_source = Path(sst_path)
    res = resolution or resolution_from_path(ic_source) or DEFAULT_RESOLUTION
    lat = int(nlat or nlat_from_path(sst_source) or 192)
    ic_target = root / _ic_target_name(
        label, res, ic_naming, control=control, member=member, sign=sign
    )
    simple_name = root / IC_TEMPLATE.format(date=label, resolution=res)
    sst_target = root / SST_TEMPLATE.format(date=label[:8], latitude=lat)
    if dry_run:
        return LinkedInputs(ic=ic_target, sst=sst_target, mode=str(mode))
    root.mkdir(parents=True, exist_ok=True)
    ic_mode = _place(ic_source, ic_target, mode)
    if ic_target != simple_name:
        # o modelo le GANLNMC...; o nome operacional e o que fica documentado
        _place(ic_target, simple_name, "link" if str(mode).lower() == "link" else mode)
    _place(sst_source, sst_target, mode)
    return LinkedInputs(ic=ic_target, sst=sst_target, mode=ic_mode)


def stage_datain(
    work_dir: str | Path,
    label: str | Path,
    suite: str | Path,
    ic_path: str | Path,
    sst_path: str | Path,
    *,
    labeli: str | datetime | None = None,
    resolution: str | None = None,
    mode: str = "link",
    ic_naming: str = "simple",
    control: bool = False,
    member: int | None = None,
    sign: str | None = None,
    dry_run: bool = False,
) -> Path:
    """Cria um datain por membro (espelho + IC/SST do membro).

    Cada job aponta path_in para o seu proprio datain.<label>; assim os
    2N+1 membros do MB09 podem rodar simultaneamente sem sobrescrever o
    GANLNMC e o SSTMonthlyDirec uns dos outros. ic_naming='operacional'
    grava GANLCTR/GANL{NN}P/GANL{NN}N (com alias GANLNMC para o modelo).
    """
    safe = _safe_label(label)
    dest = Path(work_dir) / ("datain.%s" % safe)
    if dry_run:
        return dest
    dest.mkdir(parents=True, exist_ok=True)
    source_root = Path(suite) / "model" / "datain"
    if source_root.is_dir():
        for entry in sorted(source_root.iterdir()):
            if entry.name.startswith(("GANL", "SSTMonthlyDirec", "SSTWeekly")):
                continue
            target = dest / entry.name
            if not target.exists() and not target.is_symlink():
                try:
                    os.symlink(entry.resolve(), target)
                except OSError:
                    if entry.is_dir():
                        shutil.copytree(entry, target)
                    else:
                        shutil.copy2(entry, target)
    link_inputs(
        suite,
        ic_path,
        sst_path,
        labeli or label,
        resolution=resolution,
        mode=mode,
        datain=dest,
        ic_naming=ic_naming,
        control=control,
        member=member,
        sign=sign,
    )
    return dest


# ---------------------------------------------------------------------------
# PBS: submissao e espera
# ---------------------------------------------------------------------------
def _resolve_tool(preferred: str, fallback: str) -> str:
    if preferred and Path(preferred).exists():
        return preferred
    found = shutil.which(fallback)
    return found or preferred or fallback


def submit(
    pbs_path: str | Path,
    block: bool = False,
    *,
    qsub: str = DEFAULT_QSUB,
    dry_run: bool = False,
    cwd: str | Path | None = None,
) -> str:
    """Submete o script PBS e devolve o job id.

    block=True usa qsub -W block=true (o qsub so retorna quando o job termina).
    dry_run devolve a linha de comando entre parenteses.
    """
    tool = _resolve_tool(qsub, "qsub")
    cmd = [tool]
    if block:
        cmd += ["-W", "block=true"]
    cmd.append(str(pbs_path))
    printable = " ".join(_q(part) for part in cmd)
    if dry_run:
        return "(dry-run) %s" % printable
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError("qsub falhou (%d): %s" % (proc.returncode, proc.stderr.strip()))
    text = proc.stdout.strip()
    if not text:
        raise RuntimeError("qsub nao devolveu job id")
    return text.splitlines()[0].split()[0]


#: estados PBS que indicam job terminado (C = complete, F = finished)
_FINAL_STATES = frozenset({"C", "F"})
_RUNNING_STATES = frozenset({"R", "E", "S"})
_WAITING_STATES = frozenset({"Q", "W", "T", "M"})
_HELD_STATES = frozenset({"H"})


@dataclass(frozen=True)
class JobResult:
    """Resultado de wait_job."""

    job_id: str
    state: str
    exit_status: int | None = None
    disappeared: bool = False
    polls: int = 0


def parse_qstat_state(text: str, job_id: str | None = None) -> str | None:
    """Extrai o estado (R, Q, C...) de uma saida de qstat.

    Aceita o formato de uma linha (PBS Pro/Torque), qstat -f (job_state = R) e
    reconhece o job sumido (Unknown Job Id).
    """
    if not text or not text.strip():
        return None
    if re.search(r"(?i)(unknown job|no such job|invalid job|cannot find job|job not found)", text):
        return None
    match = re.search(r"(?im)^\s*job_state\s*=\s*([A-Za-z])", text)
    if match:
        return match.group(1).upper()
    states = _FINAL_STATES | _RUNNING_STATES | _WAITING_STATES | _HELD_STATES
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("Job id", "Job_ID", "---")):
            continue
        if job_id and not stripped.startswith(str(job_id)):
            continue
        tokens = stripped.split()
        # colunas do qstat de uma linha: ... [S] [Time]; procura do fim
        for token in reversed(tokens):
            if len(token) == 1 and token.upper() in states:
                return token.upper()
    return None


def _parse_exit_status(text: str) -> int | None:
    match = re.search(r"(?im)^\s*exit_status\s*=\s*(-?\d+)", text)
    return int(match.group(1)) if match else None


def wait_job(
    job_id: str,
    poll: float = 60,
    timeout: float | None = None,
    *,
    qstat: str = DEFAULT_QSTAT,
    dry_run: bool = False,
    sleep=time.sleep,
) -> JobResult:
    """Espera o job sair de Q/R.

    Considera concluido quando o estado e C/F ou quando o job desaparece do
    qstat (job purgado apos terminar). Levanta TimeoutError se timeout estourar
    e RuntimeError para job em H (held nunca roda sozinho).
    """
    if float(poll) < 0:
        raise ValueError("poll deve ser >= 0")
    tool = _resolve_tool(qstat, "qstat")
    if dry_run:
        return JobResult(job_id=job_id, state="DRY", disappeared=True)
    start = time.monotonic()
    polls = 0
    while True:
        polls += 1
        proc = subprocess.run([tool, str(job_id)], capture_output=True, text=True)
        combined = proc.stdout + "\n" + proc.stderr
        state = parse_qstat_state(combined, str(job_id))
        if proc.returncode != 0 and state is None:
            # job purgado: o qstat deixa de conhecer o id depois do termino
            return JobResult(job_id, "C", _parse_exit_status(proc.stdout), True, polls)
        if state is None:
            return JobResult(job_id, "C", None, True, polls)
        if state in _FINAL_STATES:
            return JobResult(job_id, state, _parse_exit_status(proc.stdout), False, polls)
        if state in _HELD_STATES:
            raise RuntimeError("job %s esta em H (held) e nao vai rodar" % job_id)
        if timeout is not None and (time.monotonic() - start) >= float(timeout):
            raise TimeoutError(
                "timeout de %.0fs esperando o job %s (ultimo estado %s)"
                % (float(timeout), job_id, state)
            )
        sleep(poll)


# ---------------------------------------------------------------------------
# saidas do modelo e pos-processamento
# ---------------------------------------------------------------------------
def collect_outputs(
    dataout_dir: str | Path,
    labeli: str | datetime,
    labelf: str | datetime,
    *,
    pattern: str = "GFCTNMC*",
) -> list[Path]:
    """Lista as saidas de uma execucao em <dataout>/<LABELI>-<LABELF>/."""
    label = run_id(labeli, labelf)
    base = Path(dataout_dir)
    candidates = [base / label]
    if base.name == label:
        candidates.insert(0, base)
    candidates += [base / "model" / "dataout" / label]
    run_dir = None
    for candidate in candidates:
        if candidate.is_dir():
            run_dir = candidate
            break
    if run_dir is None:
        for candidate in base.glob("**/%s" % label):
            if candidate.is_dir():
                run_dir = candidate
                break
    if run_dir is None:
        raise FileNotFoundError("saida de %s nao encontrada em %s" % (label, base))
    return sorted(
        path for path in run_dir.iterdir()
        if path.is_file() and fnmatch.fnmatch(path.name, pattern)
    )


@dataclass(frozen=True)
class PostGribResult:
    """Caminhos produzidos pelo PostGrib."""

    unf: Path
    ctl: Path
    log: Path


def render_postin(template_path: str | Path, out_path: str | Path, **values) -> Path:
    """Troca {chave}/@CHAVE@ no POSTIN do PostGrib e grava out_path."""
    text = Path(template_path).read_text(encoding="utf-8")
    pattern = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}|@([A-Za-z_][A-Za-z0-9_]*)@")

    def repl(match: re.Match) -> str:
        key = match.group(1) or match.group(2)
        return str(values[key]) if key in values else match.group(0)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(pattern.sub(repl, text), encoding="utf-8")
    return out


def run_postgrib(
    pos_dir: str | Path,
    labeli: str | datetime,
    labelf: str | datetime,
    postin_path: str | Path,
    *,
    resolution: str = DEFAULT_RESOLUTION,
    exe_name: str = "PostGrib",
    log_name: str = POSTGRIB_LOG,
    dry_run: bool = False,
) -> PostGribResult:
    """Roda <pos>/exec/PostGrib < POSTIN e devolve .unf e .ctl."""
    label_i = parse_date(labeli)
    label_f = parse_date(labelf)
    pos = Path(pos_dir)
    executable = pos / "exec" / exe_name
    unf = (
        pos / "dataout" / resolution / label_i
        / ("GPOSNMC%s%sC.unf.%s" % (label_i, label_f, resolution))
    )
    ctl = Path(str(unf) + ".ctl")
    log = pos / log_name
    if dry_run:
        return PostGribResult(unf=unf, ctl=ctl, log=log)
    if not executable.is_file():
        raise FileNotFoundError("PostGrib nao encontrado: %s" % executable)
    postin = Path(postin_path)
    if not postin.is_file():
        raise FileNotFoundError("POSTIN nao encontrado: %s" % postin)
    unf.parent.mkdir(parents=True, exist_ok=True)
    with postin.open("rb") as stdin, log.open("wb") as stdout:
        proc = subprocess.run(
            [str(executable.resolve())], stdin=stdin, stdout=stdout,
            stderr=subprocess.STDOUT, cwd=str(pos),
        )
    if proc.returncode != 0:
        raise RuntimeError("PostGrib falhou (%d); veja %s" % (proc.returncode, log))
    if not unf.exists():
        raise RuntimeError("PostGrib nao gerou %s; veja %s" % (unf, log))
    return PostGribResult(unf=unf, ctl=ctl, log=log)


def ensemble_mean_spread(
    ctl_paths: Sequence[str | Path],
    out_dir: str | Path,
    *,
    prefix: str = "ens",
    dry_run: bool = False,
) -> list[Path]:
    """Media e espalhamento (desvio-padrao) de um ensemble de arquivos GrADS.

    Fallback usado pelo scripts/mb09.sh quando bamci.score nao esta disponivel.
    Grava um <prefix>_<var>.npz por variavel com mean e spread no mesmo formato
    de GradsFile.read_volume: (nt, nz, ny, nx); variaveis 2-D declaradas como
    "pslc 0 99" saem como (nt, 1, ny, nx).
    """
    paths = [Path(p) for p in ctl_paths]
    if not paths:
        raise ValueError("nenhum .ctl de membro informado")
    if dry_run:
        return []
    import numpy as np

    from .ctl import read_ctl
    from .grads import GradsFile

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    first_ctl = read_ctl(paths[0])
    for var in first_ctl.variables:
        cubes = []
        for ctl_path in paths:
            ctl = read_ctl(ctl_path)
            if ctl.nt != first_ctl.nt:
                raise ValueError(
                    "membros com numeros de tempos diferentes: %d e %d"
                    % (first_ctl.nt, ctl.nt)
                )
            grads = GradsFile(ctl)
            cubes.append(np.stack(
                [grads.read_volume(var.name, t=t) for t in range(1, ctl.nt + 1)],
                axis=0,
            ))
        stack = np.stack(cubes, axis=0)  # (membros, nt, nz, ny, nx)
        target = out / ("%s_%s.npz" % (prefix, var.name))
        np.savez_compressed(
            target, mean=np.nanmean(stack, axis=0),
            spread=np.nanstd(stack, axis=0), members=len(paths),
        )
        written.append(target)
    return written


# ---------------------------------------------------------------------------
# CLI (usado por scripts/mb09.sh via python -m bamci.bamrun)
# ---------------------------------------------------------------------------
def _parse_set(items: Iterable[str] | None) -> dict:
    values: dict = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError("--set espera CHAVE=VALOR, recebeu %r" % (item,))
        key, value = item.split("=", 1)
        values[key.strip()] = value
    return values


def _add_dry(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--dry-run", action="store_true", help="nao executa/escreve nada")


def _add_ic_naming_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ic-naming", default="simple", choices=("simple", "operacional"),
                        help="simple: GANLNMC...; operacional: GANLCTR/GANL01P/GANL01N")
    parser.add_argument("--control", action="store_true",
                        help="IC de controle (nome GANLCTR no modo operacional)")
    parser.add_argument("--member", type=int, default=None,
                        help="indice k=1..N da perturbacao (modo operacional)")
    parser.add_argument("--sign", choices=("+", "-", "P", "N"), default=None,
                        help="sinal da perturbacao (modo operacional)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bamci.bamrun", description="driver do modelo BAM 1.2 e do ensaio MB09"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("namelist", help="renderiza o MODELIN")
    p.add_argument("--template", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--date", default=None, help="YYYYMMDDHH (idatei e idatef)")
    p.add_argument("--idatei", default=None)
    p.add_argument("--idatef", default=None)
    p.add_argument("--dhfct", type=int, default=None)
    p.add_argument("--trunc", type=int, default=None)
    p.add_argument("--vert", type=int, default=None)
    p.add_argument("--dt", type=int, default=None)
    p.add_argument("--path-in", default=None)
    p.add_argument("--dir-output", default=None)
    p.add_argument("--set", dest="sets", action="append", default=None,
                   help="CHAVE=VALOR adicional (repetivel)")
    _add_dry(p)

    p = sub.add_parser("pbs", help="renderiza/grava o script PBS do modelo")
    p.add_argument("--out", default=None, help="sem --out imprime na tela")
    p.add_argument("--job-name", required=True)
    p.add_argument("--workdir", default=".")
    p.add_argument("--exec-dir", required=True)
    p.add_argument("--modelin", default=MODELIN_NAME)
    p.add_argument("--aprun-args", default=None)
    p.add_argument("--queue", default=DEFAULT_QUEUE)
    p.add_argument("--account", default=DEFAULT_ACCOUNT)
    p.add_argument("--walltime", default=DEFAULT_WALLTIME)
    p.add_argument("--nodes", type=int, default=None)
    p.add_argument("--ppn", type=int, default=None)
    p.add_argument("--omp-threads", type=int, default=DEFAULT_OMP_THREADS)
    p.add_argument("--env-source", default=DEFAULT_ENV_SOURCE)
    p.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    p.add_argument("--ranks-per-node", type=int, default=DEFAULT_RANKS_PER_NODE)
    p.add_argument("--pbs-style", default="mpp", choices=("mpp", "nodes", "select"),
                   help="sintaxe do recurso PBS (default mpp para o ensaio)")
    p.add_argument("--mppwidth", type=int, default=None,
                   help="PBS: #PBS -l mppwidth=N (ensaio operacional)")
    p.add_argument("--mppnppn", type=int, default=None,
                   help="PBS: #PBS -l mppnppn=N")
    p.add_argument("--mppdepth", type=int, default=None,
                   help="PBS: #PBS -l mppdepth=N")
    _add_dry(p)

    p = sub.add_parser("prepare", help="cria exec_NMC.<label> e copia o executavel")
    p.add_argument("--work", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--executable", default=DEFAULT_MODEL_EXE)
    p.add_argument("--modelin-name", default=MODELIN_NAME)
    _add_dry(p)

    p = sub.add_parser("link", help="cria GANLNMC... e SSTMonthlyDirec...")
    p.add_argument("--suite", default=DEFAULT_SUITE)
    p.add_argument("--ic", required=True)
    p.add_argument("--sst", required=True)
    p.add_argument("--labeli", required=True)
    p.add_argument("--resolution", default=None)
    p.add_argument("--mode", default="link", choices=("link", "copy"))
    p.add_argument("--datain", default=None)
    p.add_argument("--nlat", type=int, default=None)
    _add_ic_naming_args(p)
    _add_dry(p)

    p = sub.add_parser("stage", help="cria um datain por membro do ensemble")
    p.add_argument("--suite", default=DEFAULT_SUITE)
    p.add_argument("--work", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--ic", required=True)
    p.add_argument("--sst", required=True)
    p.add_argument("--labeli", default=None)
    p.add_argument("--resolution", default=None)
    p.add_argument("--mode", default="link", choices=("link", "copy"))
    _add_ic_naming_args(p)
    _add_dry(p)

    p = sub.add_parser("submit", help="submete o script PBS")
    p.add_argument("pbs_path")
    p.add_argument("--block", action="store_true")
    p.add_argument("--qsub", default=DEFAULT_QSUB)
    _add_dry(p)

    p = sub.add_parser("wait", help="espera um job PBS terminar")
    p.add_argument("job_id")
    p.add_argument("--poll", type=float, default=60.0)
    p.add_argument("--timeout", type=float, default=None)
    p.add_argument("--qstat", default=DEFAULT_QSTAT)
    _add_dry(p)

    p = sub.add_parser("collect", help="lista as saidas de uma execucao")
    p.add_argument("--dataout", required=True)
    p.add_argument("--labeli", required=True)
    p.add_argument("--labelf", required=True)
    p.add_argument("--pattern", default="GFCTNMC*")

    p = sub.add_parser("postin", help="renderiza o POSTIN do PostGrib")
    p.add_argument("--template", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--set", dest="sets", action="append", default=None)

    p = sub.add_parser("postgrib", help="roda o PostGrib e devolve .unf/.ctl")
    p.add_argument("--pos-dir", default=None,
                   help="<suite>/pos (default: DEFAULT_SUITE/pos)")
    p.add_argument("--labeli", required=True)
    p.add_argument("--labelf", required=True)
    p.add_argument("--postin", required=True)
    p.add_argument("--resolution", default=DEFAULT_RESOLUTION)
    _add_dry(p)

    p = sub.add_parser("ensemble", help="media/espalhamento (fallback sem score)")
    p.add_argument("--out", required=True)
    p.add_argument("--prefix", default="ens")
    p.add_argument("ctl", nargs="+")
    _add_dry(p)

    p = sub.add_parser("shift", help="soma horas a uma data YYYYMMDDHH")
    p.add_argument("--date", required=True)
    p.add_argument("--hours", type=int, required=True)

    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "namelist":
            overrides = _parse_set(args.sets)
            for key in ("dhfct", "trunc", "vert", "dt"):
                value = getattr(args, key)
                if value is not None:
                    overrides[key] = value
            if args.idatei is not None:
                overrides["idatei"] = args.idatei
            if args.idatef is not None:
                overrides["idatef"] = args.idatef
            if args.path_in is not None:
                overrides["path_in"] = args.path_in
            if args.dir_output is not None:
                overrides["dirfNameOutput"] = args.dir_output
            if args.date is not None:
                overrides.setdefault("idatei", args.date)
                overrides.setdefault("idatef", args.date)
            if args.dry_run:
                print("(dry-run) render_modelin %s -> %s (%s)"
                      % (args.template, args.out, sorted(overrides)))
            else:
                print(render_modelin(args.template, args.out, **overrides))
        elif args.command == "pbs":
            if args.dry_run:
                print("(dry-run) render_pbs %s" % args.job_name)
            else:
                text = render_pbs(
                    job_name=args.job_name,
                    workdir=args.workdir,
                    exec_dir=args.exec_dir,
                    modelin=args.modelin,
                    aprun_args=args.aprun_args,
                    queue=args.queue,
                    account=args.account,
                    walltime=args.walltime,
                    nodes=args.nodes,
                    ppn=args.ppn,
                    omp_threads=args.omp_threads,
                    env_source=args.env_source,
                    ranks=args.ranks,
                    ranks_per_node=args.ranks_per_node,
                    pbs_style=args.pbs_style,
                    mppwidth=args.mppwidth,
                    mppnppn=args.mppnppn,
                    mppdepth=args.mppdepth,
                )
                if args.out:
                    print(write_pbs(args.out, text))
                else:
                    print(text)
        elif args.command == "prepare":
            result = prepare_run_dir(
                args.work, args.label, args.executable,
                modelin_name=args.modelin_name, dry_run=args.dry_run,
            )
            print(result.run_dir)
        elif args.command == "link":
            result = link_inputs(
                args.suite, args.ic, args.sst, args.labeli,
                resolution=args.resolution, mode=args.mode,
                datain=args.datain, nlat=args.nlat,
                ic_naming=args.ic_naming, control=args.control,
                member=args.member, sign=args.sign, dry_run=args.dry_run,
            )
            print(result.ic)
            print(result.sst)
        elif args.command == "stage":
            print(stage_datain(
                args.work, args.label, args.suite, args.ic, args.sst,
                labeli=args.labeli, resolution=args.resolution,
                mode=args.mode, ic_naming=args.ic_naming,
                control=args.control, member=args.member, sign=args.sign,
                dry_run=args.dry_run,
            ))
        elif args.command == "submit":
            print(submit(args.pbs_path, block=args.block,
                         qsub=args.qsub, dry_run=args.dry_run))
        elif args.command == "wait":
            result = wait_job(
                args.job_id, poll=args.poll, timeout=args.timeout,
                qstat=args.qstat, dry_run=args.dry_run,
            )
            print("%s %s%s" % (
                result.job_id, result.state,
                " (desapareceu)" if result.disappeared else "",
            ))
        elif args.command == "collect":
            for path in collect_outputs(
                args.dataout, args.labeli, args.labelf, pattern=args.pattern
            ):
                print(path)
        elif args.command == "postin":
            print(render_postin(args.template, args.out, **_parse_set(args.sets)))
        elif args.command == "postgrib":
            pos_dir = args.pos_dir or str(Path(DEFAULT_SUITE) / "pos")
            result = run_postgrib(
                pos_dir, args.labeli, args.labelf, args.postin,
                resolution=args.resolution, dry_run=args.dry_run,
            )
            print(result.unf)
            print(result.ctl)
        elif args.command == "ensemble":
            for path in ensemble_mean_spread(
                args.ctl, args.out, prefix=args.prefix, dry_run=args.dry_run
            ):
                print(path)
        elif args.command == "shift":
            print(shift_date(args.date, args.hours))
        else:  # pragma: no cover - argparse garante o subcomando
            raise ValueError("subcomando desconhecido: %s" % args.command)
    except (ValueError, KeyError, FileNotFoundError, RuntimeError, TimeoutError) as exc:
        print("bamrun: erro: %s" % exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
