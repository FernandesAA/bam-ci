#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# mb09.sh - ensaio MB09 de ponta a ponta no BAM 1.2 subsazonal (JACI)
#
# Passos:
#   0 ambiente + WORK
#   1 condicao inicial de controle (bamci run)
#   2 SST semanal + SSTMonthlyDirec
#   3 N ICs aleatorias (bamci.mb09 random, N impar)
#   4 submete controle + N (DHFCT=3, 36 h) e espera
#   5 pos-processa as N+1 saidas (PostGrib)
#   6 constroi EOFs e escreve os 2N+1 ICs finais (bamci.mb09 build/apply)
#   7 submete as 2N+1 previsoes ate LEAD_DAYS e espera
#   8 pos-processa as 2N+1 e monta media/espalhamento (+ score, se existir)
#   9 resumo
#
# Uso:
#   DRY=1 scripts/mb09.sh                 # so imprime os comandos
#   DRY=0 scripts/mb09.sh --steps 0,1,2   # roda um trecho
#   DRY=0 scripts/mb09.sh --steps 4,5 --skip-ic
#
# Nada e apagado com rm -rf: o lixo vai para WORK/trash.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- configuracao (sobrescreva por ambiente) -------------------------------
SUITE_ROOT="${SUITE_ROOT:-/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper}"
WORK="${WORK:-$PWD/mb09_work}"
DATE="${DATE:-$(date -u +%Y%m%d%H)}"
TRUNC="${TRUNC:-126}"
VERT="${VERT:-42}"
PRESET="${PRESET:-etsa}"
LEAD_DAYS="${LEAD_DAYS:-30}"
# defaults do ensaio de perturbacao (operacional oensMB09_bam_xc50)
QUEUE="${QUEUE:-pesq}"
ACCOUNT="${ACCOUNT:-CPTEC}"
WALLTIME="${WALLTIME:-02:00:00}"
MPPWIDTH="${MPPWIDTH:-48}"
MPPNPPN="${MPPNPPN:-4}"
MPPDEPTH="${MPPDEPTH:-6}"
OMP_THREADS="${OMP_THREADS:-6}"
RANKS="${RANKS:-48}"
RANKS_PER_NODE="${RANKS_PER_NODE:-4}"
# previsoes longas: outra fila/walltime (sobrescreva por ambiente)
QUEUE_FCST="${QUEUE_FCST:-workq}"
WALLTIME_FCST="${WALLTIME_FCST:-23:00:00}"
# mpp: oensMB09 | nodes: bam_SeasonalClima | select: PBS Pro
PBS_STYLE="${PBS_STYLE:-mpp}"
NODES="${NODES:-$(( (RANKS + RANKS_PER_NODE - 1) / RANKS_PER_NODE ))}"
PPN="${PPN:-$RANKS_PER_NODE}"
DRY="${DRY:-1}"
STEPS="${STEPS:-0,1,2,3,4,5,6,7,8,9}"

CASE="${CASE:-ERA5}"
GRIB="${GRIB:-}"
DELTA_SIGMA="${DELTA_SIGMA:-$SUITE_ROOT/pre/datain/DeltaSigma.L$(printf '%03d' "$VERT")}"
SST_SOURCE="${SST_SOURCE:-oisst}"
SST_INPUT="${SST_INPUT:-}"
SST_DAYS="${SST_DAYS:-7}"
NLAT="${NLAT:-192}"
SEED="${SEED:-0}"
DT="${DT:-450}"
DHFCT_MB09="${DHFCT_MB09:-3}"
DHFCT_FCST="${DHFCT_FCST:-24}"
# N_PERTURBATIONS tem de ser IMPAR (membros = 2*N+1); N_RANDOM fica como alias antigo
N_PERTURBATIONS="${N_PERTURBATIONS:-${N_RANDOM:-7}}"
IC_NAMING="${IC_NAMING:-simple}"
SKIP_IC="${SKIP_IC:-0}"
PS_UNITS="${PS_UNITS:-hpa}"
MB09_VARS="${MB09_VARS:-}"

MODEL_EXE="${MODEL_EXE:-/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper/model/exec/ParModel_MPI}"
ENV_SOURCE="${ENV_SOURCE:-/p/projetos/ioper/fix/sazonal/load_modulei_intel_2022_oneapi_mpi.bash}"
QSUB="${QSUB:-/opt/pbs/default/bin/qsub}"
QSTAT="${QSTAT:-/opt/pbs/bin/qstat}"
# resolvidos em resolve_templates (abaixo), aceitando o caminho real da suite
MODELIN_TEMPLATE="${MODELIN_TEMPLATE:-}"
POSTIN_TEMPLATE="${POSTIN_TEMPLATE:-}"
PYTHON="${PYTHON:-python3}"

# --- argumentos de linha de comando ----------------------------------------
usage() {
  cat <<'USAGE'
uso: mb09.sh [--steps LISTA] [--date AAAAMMDDHH] [--work DIR] [--dry|--no-dry]
             [--skip-ic] [--ic-naming {simple,operacional}]

  --steps LISTA  passos separados por virgula (ex.: 0,1,2 ou 4,5)
  --date DATA    data de analise AAAAMMDDHH (default: agora em UTC)
  --work DIR     diretorio de trabalho (default: $PWD/mb09_work)
  --dry          apenas imprime os comandos (default)
  --no-dry       executa de verdade
  --skip-ic      reaproveita a IC de controle se ela ja existir
  --ic-naming    simple (GANLNMC...) ou operacional (GANLCTR/GANL01P/GANL01N)

Variaveis principais: SUITE_ROOT, WORK, DATE, TRUNC, VERT, PRESET, LEAD_DAYS,
N_PERTURBATIONS (impar, default 7), IC_NAMING, QUEUE, ACCOUNT, WALLTIME,
PBS_STYLE (mpp|nodes|select), MPPWIDTH, MPPNPPN, MPPDEPTH, OMP_THREADS, RANKS,
RANKS_PER_NODE, QUEUE_FCST, WALLTIME_FCST, GRIB, DELTA_SIGMA, SST_INPUT,
SST_SOURCE, MODELIN_TEMPLATE, POSTIN_TEMPLATE, MODEL_EXE, ENV_SOURCE, PYTHON,
DRY, STEPS.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --steps) STEPS="$2"; shift 2 ;;
    --steps=*) STEPS="${1#*=}"; shift ;;
    --date) DATE="$2"; shift 2 ;;
    --date=*) DATE="${1#*=}"; shift ;;
    --work) WORK="$2"; shift 2 ;;
    --work=*) WORK="${1#*=}"; shift ;;
    --dry) DRY=1; shift ;;
    --no-dry) DRY=0; shift ;;
    --skip-ic) SKIP_IC=1; shift ;;
    --ic-naming) IC_NAMING="$2"; shift 2 ;;
    --ic-naming=*) IC_NAMING="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "mb09.sh: argumento desconhecido: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# --- validacao de parametros -----------------------------------------------
if ! [[ "$N_PERTURBATIONS" =~ ^[0-9]+$ ]] || (( N_PERTURBATIONS < 1 )); then
  echo "ERRO: N_PERTURBATIONS deve ser um inteiro positivo (recebido '$N_PERTURBATIONS')" >&2
  exit 2
fi
if (( N_PERTURBATIONS % 2 == 0 )); then
  echo "ERRO: N_PERTURBATIONS deve ser IMPAR (membros = 2*N+1); recebido $N_PERTURBATIONS" >&2
  exit 2
fi
case "$IC_NAMING" in
  simple|operacional) ;;
  *)
    echo "ERRO: IC_NAMING deve ser 'simple' ou 'operacional' (recebido '$IC_NAMING')" >&2
    exit 2
    ;;
esac
case "$PBS_STYLE" in
  mpp|nodes|select) ;;
  *)
    echo "ERRO: PBS_STYLE deve ser 'mpp', 'nodes' ou 'select' (recebido '$PBS_STYLE')" >&2
    exit 2
    ;;
esac
N_MEMBERS=$((2 * N_PERTURBATIONS + 1))

# --- derivados --------------------------------------------------------------
RESOLUTION="TQ$(printf '%04d' "$TRUNC")L$(printf '%03d' "$VERT")"
DATE8="${DATE:0:8}"
IC_CONTROL="$SUITE_ROOT/model/datain/GANL${CASE}${DATE}S.unf.${RESOLUTION}"
SST_WEEKLY="$SUITE_ROOT/model/datain/SSTWeekly${DATE8}.G$(printf '%05d' "$NLAT")"
MB09_DIR="$WORK/mb09"
RANDOM_DIR="$MB09_DIR/random"
FINAL_DIR="$MB09_DIR/ic"
RANDOM_JSON="$RANDOM_DIR/random.json"
FINAL_JSON="$FINAL_DIR/members.json"
CTLS=()
BUILD_ARGS=()
JOBS=()
RANDOM_ICS=()
FINAL_ICS=()

# --- helpers ----------------------------------------------------------------
date_shift() {
  local d="$1" h="$2"
  # GNU date: "AAAA-MM-DD HH:00 <n> hours" (sem o '+' inicial)
  date -u -d "${d:0:4}-${d:4:2}-${d:6:2} ${d:8:2}:00 ${h} hours" +%Y%m%d%H
}

label_after() {
  date_shift "$DATE" "$1"
}

want_step() {
  [[ ",$STEPS," == *",$1,"* ]]
}

log_step() {
  printf '\n===== PASSO %s: %s =====\n' "$1" "$2"
}

run() {
  if [[ "$DRY" == "1" ]]; then
    printf 'DRY: ' >&2
    printf '%q ' "$@" >&2
    printf '\n' >&2
    return 0
  fi
  printf '+ ' >&2
  printf '%q ' "$@" >&2
  printf '\n' >&2
  "$@"
}

bamci_cmd() {
  if [[ -n "${BAMCI:-}" ]]; then
    "$BAMCI" "$@"
  elif command -v bamci >/dev/null 2>&1; then
    bamci "$@"
  else
    "$PYTHON" -m bamci "$@"
  fi
}

py_bamrun() {
  "$PYTHON" -m bamci.bamrun "$@"
}

check_file() {
  if [[ -e "$1" ]]; then
    return 0
  fi
  if [[ "$DRY" == "1" ]]; then
    echo "AVISO (dry): nao encontrado $1" >&2
  else
    echo "ERRO: nao encontrado $1" >&2
    exit 1
  fi
}

check_cmd() {
  if command -v "$1" >/dev/null 2>&1; then
    return 0
  fi
  if [[ "$DRY" == "1" ]]; then
    echo "AVISO (dry): comando '$1' fora do PATH" >&2
  else
    echo "ERRO: comando '$1' nao encontrado" >&2
    exit 1
  fi
}

trash() {
  local target="$1"
  [[ -e "$target" ]] || return 0
  run mkdir -p "$WORK/trash"
  run mv "$target" "$WORK/trash/$(basename "$target").$(date -u +%Y%m%d%H%M%S)"
}

launch_run() {
  local tag="$1" ic="$2" labelf="$3" dhfct="$4"
  local control="${5:-0}" member="${6:-}" sign="${7:-}"
  local run_dir="$WORK/runs/exec_NMC.$tag"
  local datain_dir="$WORK/datain/$tag"
  local out_dir="$WORK/out/$tag"
  local pbs_file="$WORK/pbs/$tag.pbs"
  local queue="${RUN_QUEUE:-$QUEUE}"
  local walltime="${RUN_WALLTIME:-$WALLTIME}"
  local stage_args=(--suite "$SUITE_ROOT" --work "$WORK/datain" \
      --label "$tag" --ic "$ic" --sst "$SST_WEEKLY" --labeli "$DATE" \
      --resolution "$RESOLUTION" --mode link)
  local job=""

  if [[ "$IC_NAMING" == "operacional" ]]; then
    stage_args+=(--ic-naming operacional)
    if [[ "$control" == "1" ]]; then
      stage_args+=(--control)
    else
      stage_args+=(--member "$member" --sign "$sign")
    fi
  fi
  run py_bamrun stage "${stage_args[@]}"
  run py_bamrun prepare --work "$WORK/runs" --label "$tag" \
      --executable "$MODEL_EXE"
  run py_bamrun namelist --template "$MODELIN_TEMPLATE" \
      --out "$run_dir/MODELIN" --date "$DATE" --idatef "$labelf" \
      --dhfct "$dhfct" --trunc "$TRUNC" --vert "$VERT" --dt "$DT" \
      --path-in "$datain_dir" --dir-output "$out_dir"
  run py_bamrun pbs --out "$pbs_file" --job-name "bam-$tag" \
      --workdir "$WORK" --exec-dir "$run_dir" --queue "$queue" \
      --account "$ACCOUNT" --walltime "$walltime" --nodes "$NODES" \
      --ppn "$PPN" --omp-threads "$OMP_THREADS" --env-source "$ENV_SOURCE" \
      --ranks "$RANKS" --ranks-per-node "$RANKS_PER_NODE" --pbs-style "$PBS_STYLE" \
      --mppwidth "$MPPWIDTH" --mppnppn "$MPPNPPN" --mppdepth "$MPPDEPTH"
  job="$(run py_bamrun submit --qsub "$QSUB" "$pbs_file")"
  if [[ -n "$job" ]]; then
    echo "  job $tag -> $job"
    JOBS+=("$job")
  fi
}

wait_jobs() {
  local job
  if [[ "${#JOBS[@]}" -eq 0 ]]; then
    echo "  (nenhum job submetido para esperar)"
    return 0
  fi
  for job in "${JOBS[@]}"; do
    run py_bamrun wait "$job" --poll 60 --timeout 172800 --qstat "$QSTAT"
  done
  JOBS=()
}

setup_post_dir() {
  local pp_dir="$1" entry name
  run mkdir -p "$pp_dir/dataout"
  if [[ -d "$SUITE_ROOT/pos" ]]; then
    for entry in "$SUITE_ROOT/pos"/*; do
      name="$(basename "$entry")"
      [[ "$name" == "dataout" ]] && continue
      if [[ ! -e "$pp_dir/$name" ]]; then
        run ln -s "$entry" "$pp_dir/$name"
      fi
    done
  fi
}

postprocess_run() {
  local tag="$1" labelf="$2"
  local pp_dir="$WORK/post/$tag"
  local postin="$WORK/postin/$tag.POSTIN"
  local ctl
  ctl="$(expected_ctl "$tag" "$labelf")"

  run mkdir -p "$WORK/postin"
  setup_post_dir "$pp_dir"
  run py_bamrun collect --dataout "$WORK/out/$tag" \
      --labeli "$DATE" --labelf "$labelf"
  run py_bamrun postin --template "$POSTIN_TEMPLATE" --out "$postin" \
      --set "LABELI=$DATE" --set "LABELF=$labelf" \
      --set "RESOLUTION=$RESOLUTION" --set "MEMBER=$tag" \
      --set "DATAOUT=$WORK/out/$tag" --set "POSOUT=$pp_dir/dataout"
  run py_bamrun postgrib --pos-dir "$pp_dir" --labeli "$DATE" \
      --labelf "$labelf" --postin "$postin" --resolution "$RESOLUTION"
  CTLS+=("$ctl")
}

load_random_ics() {
  local k kk
  RANDOM_ICS=()
  if [[ -f "$RANDOM_JSON" ]]; then
    mapfile -t RANDOM_ICS < <("$PYTHON" -c '
import json, os, sys
d = json.load(open(sys.argv[1]))
base = os.path.dirname(os.path.abspath(sys.argv[1]))
for item in d["members"]:
    print(item if os.path.isabs(item) else os.path.join(base, item))
' "$RANDOM_JSON")
    if [[ "${#RANDOM_ICS[@]}" -eq 0 ]]; then
      echo "ERRO: random.json sem membros: $RANDOM_JSON" >&2
      exit 1
    fi
  else
    for k in $(seq 1 "$N_PERTURBATIONS"); do
      kk="$(printf '%02d' "$k")"
      RANDOM_ICS+=("$RANDOM_DIR/mb09rand_${kk}.unf.${RESOLUTION}")
    done
  fi
}

expected_ctl() {
  local tag="$1" labelf="$2"
  echo "$WORK/post/$tag/dataout/$RESOLUTION/$DATE/GPOSNMC${DATE}${labelf}C.unf.${RESOLUTION}.ctl"
}

load_final_ics() {
  local k idx
  FINAL_ICS=()
  if [[ -f "$FINAL_JSON" ]]; then
    mapfile -t FINAL_ICS < <("$PYTHON" -c '
import json, os, sys
d = json.load(open(sys.argv[1]))
base = os.path.dirname(os.path.abspath(sys.argv[1]))
def full(p):
    return p if os.path.isabs(p) else os.path.join(base, p)
print(full(d["control"]))
for item in d["members"]:
    print(full(item["path"]))
' "$FINAL_JSON")
    if [[ "${#FINAL_ICS[@]}" -eq 0 ]]; then
      echo "ERRO: members.json sem membros: $FINAL_JSON" >&2
      exit 1
    fi
  else
    for k in $(seq 0 $((2 * N_PERTURBATIONS))); do
      idx="$(printf '%03d' "$k")"
      FINAL_ICS+=("$FINAL_DIR/GANLMB09_${idx}.unf.${RESOLUTION}")
    done
  fi
}

resolve_file() {
  local fallback="$1" candidate
  shift
  for candidate in "$@"; do
    if [[ -f "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  echo "$fallback"
}

if [[ -z "$MODELIN_TEMPLATE" ]]; then
  MODELIN_TEMPLATE="$(resolve_file "$SUITE_ROOT/model/exec/MODELIN" \
      "$SUITE_ROOT/model/exec/MODELIN" "$SUITE_ROOT/model/MODELIN" \
      "$SUITE_ROOT/pre/MODELIN" "$SUITE_ROOT/pre/exec/MODELIN")"
fi
if [[ -z "$POSTIN_TEMPLATE" ]]; then
  POSTIN_TEMPLATE="$(resolve_file "$SUITE_ROOT/pos/POSTIN-GRIB" \
      "$SUITE_ROOT/pos/POSTIN-GRIB" "$SUITE_ROOT/pos/exec/POSTIN-GRIB" \
      "$SUITE_ROOT/pos/POSTIN")"
fi

# ===========================================================================
# PASSO 0 - ambiente e WORK
# ===========================================================================
if want_step 0; then
  log_step 0 "ambiente e criacao de WORK"
  check_cmd "$PYTHON"
  check_cmd aprun
  check_file "$MODEL_EXE"
  check_file "$ENV_SOURCE"
  check_file "$MODELIN_TEMPLATE"
  check_file "$QSUB"
  check_file "$QSTAT"
  if want_step 1 && [[ "$SKIP_IC" != "1" && -z "$GRIB" ]]; then
    if [[ "$DRY" == "1" ]]; then
      echo "AVISO (dry): GRIB vazio; o passo 1 nao tem analise de entrada" >&2
    else
      echo "ERRO: defina GRIB=<arquivo> ou use --skip-ic" >&2
      exit 1
    fi
  fi
  if want_step 2 && [[ -z "$SST_INPUT" ]]; then
    if [[ "$DRY" == "1" ]]; then
      echo "AVISO (dry): SST_INPUT vazio; o passo 2 nao tem SST de entrada" >&2
    else
      echo "ERRO: defina SST_INPUT=<arquivo ou glob>" >&2
      exit 1
    fi
  fi
  run mkdir -p "$WORK" "$WORK/runs" "$WORK/datain" "$WORK/out" \
      "$WORK/pbs" "$WORK/post" "$WORK/postin" "$WORK/trash" \
      "$MB09_DIR"
fi

# ===========================================================================
# PASSO 1 - condicao inicial de controle
# ===========================================================================
if want_step 1; then
  log_step 1 "condicao inicial de controle (bamci run)"
  if [[ "$SKIP_IC" == "1" || -f "$IC_CONTROL" ]]; then
    echo "  reaproveitando IC de controle: $IC_CONTROL"
  else
    run bamci_cmd run --case "$CASE" --date "$DATE" -r "$RESOLUTION" \
        --grib "$GRIB" --delta-sigma "$DELTA_SIGMA" \
        --out-dir "$SUITE_ROOT/model/datain"
  fi
fi

# ===========================================================================
# PASSO 2 - SST semanal + SSTMonthlyDirec
# ===========================================================================
if want_step 2; then
  log_step 2 "SST semanal (bamci sst) e nome lido pelo modelo"
  run bamci_cmd sst --source "$SST_SOURCE" --input "$SST_INPUT" \
      --date "$DATE" --days "$SST_DAYS" -r "$RESOLUTION" --out "$SST_WEEKLY"
  run py_bamrun link --suite "$SUITE_ROOT" --ic "$IC_CONTROL" \
      --sst "$SST_WEEKLY" --labeli "$DATE" --resolution "$RESOLUTION" \
      --mode link --ic-naming "$IC_NAMING" --control
fi

# ===========================================================================
# PASSO 3 - 7 ICs aleatorias
# ===========================================================================
if want_step 3; then
  log_step 3 "$N_PERTURBATIONS ICs aleatorias (bamci.mb09 random)"
  run mkdir -p "$RANDOM_DIR" "$MB09_DIR/forecast"
  run "$PYTHON" -m bamci.mb09 random --ic "$IC_CONTROL" \
      --out-dir "$RANDOM_DIR" --n-modes "$N_PERTURBATIONS" --seed "$SEED" \
      --preset "$PRESET" --resolution "$RESOLUTION" --prefix mb09rand
fi

# ===========================================================================
# PASSO 4 - controle + 7 previsoes curtas (36 h, DHFCT=3)
# ===========================================================================
if want_step 4; then
  log_step 4 "submete controle + $N_PERTURBATIONS (36 h, DHFCT=$DHFCT_MB09)"
  LABEL36="$(label_after 36)"
  RUN_QUEUE="$QUEUE"
  RUN_WALLTIME="$WALLTIME"
  load_random_ics
  launch_run ctrl "$IC_CONTROL" "$LABEL36" "$DHFCT_MB09" 1
  for k in $(seq 1 "$N_PERTURBATIONS"); do
    kk="$(printf '%02d' "$k")"
    launch_run "rand${kk}" "${RANDOM_ICS[$((k - 1))]}" "$LABEL36" "$DHFCT_MB09" 0 "$k" "+"
  done
  wait_jobs
fi

# ===========================================================================
# PASSO 5 - pos-processa as 8 saidas
# ===========================================================================
if want_step 5; then
  log_step 5 "pos-processamento (PostGrib) das $((N_PERTURBATIONS + 1)) saidas"
  LABEL36="$(label_after 36)"
  CTLS=()
  load_random_ics
  postprocess_run ctrl "$LABEL36"
  for k in $(seq 1 "$N_PERTURBATIONS"); do
    kk="$(printf '%02d' "$k")"
    postprocess_run "rand${kk}" "$LABEL36"
  done
fi

# ===========================================================================
# PASSO 6 - EOFs e 15 ICs finais
# ===========================================================================
if want_step 6; then
  log_step 6 "EOFs (bamci.mb09 build) e $N_MEMBERS ICs finais (bamci.mb09 apply)"
  LABEL36="$(label_after 36)"
  if [[ "${#CTLS[@]}" -lt $((N_PERTURBATIONS + 1)) ]]; then
    # permite rodar --steps 6 sozinho, reaproveitando o pos do passo 5
    CTLS=("$(expected_ctl ctrl "$LABEL36")")
    for k in $(seq 1 "$N_PERTURBATIONS"); do
      kk="$(printf '%02d' "$k")"
      CTLS+=("$(expected_ctl "rand${kk}" "$LABEL36")")
    done
  fi
  if [[ "${#CTLS[@]}" -lt $((N_PERTURBATIONS + 1)) ]]; then
    echo "ERRO: esperava $((N_PERTURBATIONS + 1)) .ctl do passo 5, tenho ${#CTLS[@]}" >&2
    exit 1
  fi
  MEMBER_CTLS=("${CTLS[@]:1}")
  BUILD_ARGS=(--work "$MB09_DIR" --preset "$PRESET" --out "$MB09_DIR/eof.npz" \
      --control "${CTLS[0]}" --members "${MEMBER_CTLS[@]}" --ps-units "$PS_UNITS" \
      --analysis-time "$DATE" --n-modes "$N_PERTURBATIONS")
  if [[ -n "$MB09_VARS" ]]; then
    BUILD_ARGS+=(--vars "$MB09_VARS")
  fi
  run "$PYTHON" -m bamci.mb09 build "${BUILD_ARGS[@]}"
  run "$PYTHON" -m bamci.mb09 apply --eof "$MB09_DIR/eof.npz" \
      --ic "$IC_CONTROL" --out-dir "$FINAL_DIR" --preset "$PRESET" \
      --resolution "$RESOLUTION" --prefix GANLMB09 --n-modes "$N_PERTURBATIONS"
fi

# ===========================================================================
# PASSO 7 - 15 previsoes ate LEAD_DAYS
# ===========================================================================
if want_step 7; then
  log_step 7 "$N_MEMBERS previsoes ate $LEAD_DAYS dias (DHFCT=$DHFCT_FCST)"
  LABEL_FCST="$(label_after "$((LEAD_DAYS * 24))")"
  RUN_QUEUE="$QUEUE_FCST"
  RUN_WALLTIME="$WALLTIME_FCST"
  load_final_ics
  if [[ "${#FINAL_ICS[@]}" -ne "$N_MEMBERS" ]]; then
    echo "ERRO: esperava $N_MEMBERS ICs finais, tenho ${#FINAL_ICS[@]}" >&2
    exit 1
  fi
  for k in $(seq 0 $((N_MEMBERS - 1))); do
    idx="$(printf '%03d' "$k")"
    if (( k == 0 )); then
      launch_run "ens${idx}" "${FINAL_ICS[$k]}" "$LABEL_FCST" "$DHFCT_FCST" 1
    else
      member_idx=$(( (k + 1) / 2 ))
      if (( k % 2 == 1 )); then sign="+"; else sign="-"; fi
      launch_run "ens${idx}" "${FINAL_ICS[$k]}" "$LABEL_FCST" "$DHFCT_FCST" 0 "$member_idx" "$sign"
    fi
  done
  wait_jobs
fi

# ===========================================================================
# PASSO 8 - pos-processa as 15, media/espalhamento e metricas
# ===========================================================================
if want_step 8; then
  log_step 8 "pos-processamento das $N_MEMBERS e media/espalhamento"
  LABEL_FCST="$(label_after "$((LEAD_DAYS * 24))")"
  CTLS=()
  for k in $(seq 0 $((N_MEMBERS - 1))); do
    idx="$(printf '%03d' "$k")"
    postprocess_run "ens${idx}" "$LABEL_FCST"
  done
  if "$PYTHON" -c 'import bamci.score' >/dev/null 2>&1; then
    run "$PYTHON" -m bamci.score --out "$WORK/score" "${CTLS[@]}"
  else
    run py_bamrun ensemble --out "$WORK/mean_spread" --prefix ens "${CTLS[@]}"
  fi
fi

# ===========================================================================
# PASSO 9 - resumo
# ===========================================================================
if want_step 9; then
  log_step 9 "resumo"
  echo "  SUITE_ROOT : $SUITE_ROOT"
  echo "  WORK       : $WORK"
  echo "  DATE       : $DATE"
  echo "  RESOLUTION : $RESOLUTION"
  echo "  N_PERTURB  : $N_PERTURBATIONS (membros = $N_MEMBERS)"
  echo "  IC_NAMING  : $IC_NAMING"
  echo "  pbs style  : $PBS_STYLE"
  echo "  fila curta : $QUEUE / $WALLTIME (mppwidth=$MPPWIDTH mppnppn=$MPPNPPN mppdepth=$MPPDEPTH)"
  echo "  fila longa : $QUEUE_FCST / $WALLTIME_FCST"
  echo "  IC controle: $IC_CONTROL"
  echo "  SST        : $SST_WEEKLY"
  echo "  random     : $RANDOM_JSON"
  echo "  EOFs       : $MB09_DIR/eof.npz"
  echo "  ICs finais : $FINAL_DIR (members.json: $FINAL_JSON)"
  echo "  saidas     : $WORK/out"
  echo "  pos        : $WORK/post"
  echo "  score      : $WORK/score (ou $WORK/mean_spread)"
fi
