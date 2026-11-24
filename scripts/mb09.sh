#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# mb09.sh - ensaio MB09 de ponta a ponta no BAM 1.2 subsazonal (JACI)
#
# Passos:
#   0 ambiente + WORK
#   1 condicao inicial de controle (bamci run)
#   2 SST semanal + SSTMonthlyDirec
#   3 7 ICs aleatorias (bamci.mb09 random)
#   4 submete controle + 7 (DHFCT=3, 36 h) e espera
#   5 pos-processa as 8 saidas (PostGrib)
#   6 constroi EOFs e escreve os 15 ICs finais (bamci.mb09 build/apply)
#   7 submete as 15 previsoes ate LEAD_DAYS e espera
#   8 pos-processa as 15 e monta media/espalhamento (+ score, se existir)
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
QUEUE="${QUEUE:-workq}"
ACCOUNT="${ACCOUNT:-CPTEC}"
WALLTIME="${WALLTIME:-23:00:00}"
NODES="${NODES:-3}"
PPN="${PPN:-40}"
OMP_THREADS="${OMP_THREADS:-10}"
RANKS="${RANKS:-12}"
RANKS_PER_NODE="${RANKS_PER_NODE:-4}"
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
N_RANDOM="${N_RANDOM:-7}"
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
             [--skip-ic]

  --steps LISTA  passos separados por virgula (ex.: 0,1,2 ou 4,5)
  --date DATA    data de analise AAAAMMDDHH (default: agora em UTC)
  --work DIR     diretorio de trabalho (default: $PWD/mb09_work)
  --dry          apenas imprime os comandos (default)
  --no-dry       executa de verdade
  --skip-ic      reaproveita a IC de controle se ela ja existir

Variaveis principais: SUITE_ROOT, WORK, DATE, TRUNC, VERT, PRESET, LEAD_DAYS,
QUEUE, ACCOUNT, WALLTIME, NODES, PPN, OMP_THREADS, GRIB, DELTA_SIGMA,
SST_INPUT, SST_SOURCE, MODELIN_TEMPLATE, POSTIN_TEMPLATE, MODEL_EXE,
ENV_SOURCE, PYTHON, DRY, STEPS.
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
    -h|--help) usage; exit 0 ;;
    *) echo "mb09.sh: argumento desconhecido: $1" >&2; usage >&2; exit 2 ;;
  esac
done

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
  local run_dir="$WORK/runs/exec_NMC.$tag"
  local datain_dir="$WORK/datain/$tag"
  local out_dir="$WORK/out/$tag"
  local pbs_file="$WORK/pbs/$tag.pbs"
  local job=""

  run py_bamrun stage --suite "$SUITE_ROOT" --work "$WORK/datain" \
      --label "$tag" --ic "$ic" --sst "$SST_WEEKLY" --labeli "$DATE" \
      --resolution "$RESOLUTION" --mode link
  run py_bamrun prepare --work "$WORK/runs" --label "$tag" \
      --executable "$MODEL_EXE"
  run py_bamrun namelist --template "$MODELIN_TEMPLATE" \
      --out "$run_dir/MODELIN" --date "$DATE" --idatef "$labelf" \
      --dhfct "$dhfct" --trunc "$TRUNC" --vert "$VERT" --dt "$DT" \
      --path-in "$datain_dir" --dir-output "$out_dir"
  run py_bamrun pbs --out "$pbs_file" --job-name "bam-$tag" \
      --workdir "$WORK" --exec-dir "$run_dir" --queue "$QUEUE" \
      --account "$ACCOUNT" --walltime "$WALLTIME" --nodes "$NODES" \
      --ppn "$PPN" --omp-threads "$OMP_THREADS" --env-source "$ENV_SOURCE" \
      --ranks "$RANKS" --ranks-per-node "$RANKS_PER_NODE"
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
    for k in $(seq 1 "$N_RANDOM"); do
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
    for k in $(seq 0 14); do
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
      --mode link
fi

# ===========================================================================
# PASSO 3 - 7 ICs aleatorias
# ===========================================================================
if want_step 3; then
  log_step 3 "7 ICs aleatorias (bamci.mb09 random)"
  run mkdir -p "$RANDOM_DIR" "$MB09_DIR/forecast"
  run "$PYTHON" -m bamci.mb09 random --ic "$IC_CONTROL" \
      --out-dir "$RANDOM_DIR" --n "$N_RANDOM" --seed "$SEED" \
      --preset "$PRESET" --resolution "$RESOLUTION" --prefix mb09rand
fi

# ===========================================================================
# PASSO 4 - controle + 7 previsoes curtas (36 h, DHFCT=3)
# ===========================================================================
if want_step 4; then
  log_step 4 "submete controle + $N_RANDOM (36 h, DHFCT=$DHFCT_MB09)"
  LABEL36="$(label_after 36)"
  load_random_ics
  launch_run ctrl "$IC_CONTROL" "$LABEL36" "$DHFCT_MB09"
  for k in $(seq 1 "$N_RANDOM"); do
    kk="$(printf '%02d' "$k")"
    launch_run "rand${kk}" "${RANDOM_ICS[$((k - 1))]}" "$LABEL36" "$DHFCT_MB09"
  done
  wait_jobs
fi

# ===========================================================================
# PASSO 5 - pos-processa as 8 saidas
# ===========================================================================
if want_step 5; then
  log_step 5 "pos-processamento (PostGrib) das 8 saidas"
  LABEL36="$(label_after 36)"
  CTLS=()
  load_random_ics
  postprocess_run ctrl "$LABEL36"
  for k in $(seq 1 "$N_RANDOM"); do
    kk="$(printf '%02d' "$k")"
    postprocess_run "rand${kk}" "$LABEL36"
  done
fi

# ===========================================================================
# PASSO 6 - EOFs e 15 ICs finais
# ===========================================================================
if want_step 6; then
  log_step 6 "EOFs (bamci.mb09 build) e 15 ICs finais (bamci.mb09 apply)"
  LABEL36="$(label_after 36)"
  if [[ "${#CTLS[@]}" -lt 8 ]]; then
    # permite rodar --steps 6 sozinho, reaproveitando o pos do passo 5
    CTLS=("$(expected_ctl ctrl "$LABEL36")")
    for k in $(seq 1 "$N_RANDOM"); do
      kk="$(printf '%02d' "$k")"
      CTLS+=("$(expected_ctl "rand${kk}" "$LABEL36")")
    done
  fi
  if [[ "${#CTLS[@]}" -lt 8 ]]; then
    echo "ERRO: esperava 8 .ctl do passo 5, tenho ${#CTLS[@]}" >&2
    exit 1
  fi
  MEMBER_CTLS=("${CTLS[@]:1}")
  BUILD_ARGS=(--work "$MB09_DIR" --preset "$PRESET" --out "$MB09_DIR/eof.npz" \
      --control "${CTLS[0]}" --members "${MEMBER_CTLS[@]}" --ps-units "$PS_UNITS" \
      --analysis-time "$DATE")
  if [[ -n "$MB09_VARS" ]]; then
    BUILD_ARGS+=(--vars "$MB09_VARS")
  fi
  run "$PYTHON" -m bamci.mb09 build "${BUILD_ARGS[@]}"
  run "$PYTHON" -m bamci.mb09 apply --eof "$MB09_DIR/eof.npz" \
      --ic "$IC_CONTROL" --out-dir "$FINAL_DIR" --preset "$PRESET" \
      --resolution "$RESOLUTION" --prefix GANLMB09
fi

# ===========================================================================
# PASSO 7 - 15 previsoes ate LEAD_DAYS
# ===========================================================================
if want_step 7; then
  log_step 7 "15 previsoes ate $LEAD_DAYS dias (DHFCT=$DHFCT_FCST)"
  LABEL_FCST="$(label_after "$((LEAD_DAYS * 24))")"
  load_final_ics
  if [[ "${#FINAL_ICS[@]}" -ne 15 ]]; then
    echo "ERRO: esperava 15 ICs finais, tenho ${#FINAL_ICS[@]}" >&2
    exit 1
  fi
  for k in $(seq 0 14); do
    idx="$(printf '%03d' "$k")"
    launch_run "ens${idx}" "${FINAL_ICS[$k]}" "$LABEL_FCST" "$DHFCT_FCST"
  done
  wait_jobs
fi

# ===========================================================================
# PASSO 8 - pos-processa as 15, media/espalhamento e metricas
# ===========================================================================
if want_step 8; then
  log_step 8 "pos-processamento das 15 e media/espalhamento"
  LABEL_FCST="$(label_after "$((LEAD_DAYS * 24))")"
  CTLS=()
  for k in $(seq 0 14); do
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
  echo "  IC controle: $IC_CONTROL"
  echo "  SST        : $SST_WEEKLY"
  echo "  random     : $RANDOM_JSON"
  echo "  EOFs       : $MB09_DIR/eof.npz"
  echo "  ICs finais : $FINAL_DIR (members.json: $FINAL_JSON)"
  echo "  saidas     : $WORK/out"
  echo "  pos        : $WORK/post"
  echo "  score      : $WORK/score (ou $WORK/mean_spread)"
fi
