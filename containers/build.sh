#!/bin/bash
# Constroi a imagem apptainer do bam-ci.
#   ./containers/build.sh            -> bam-ci.sif no diretorio atual
#   APPTAINER="apptainer --fakeroot" ./containers/build.sh
set -euo pipefail
APPTAINER="${APPTAINER:-apptainer}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
$APPTAINER build --force bam-ci.sif containers/bam-ci.def
echo "imagem: $ROOT/bam-ci.sif"
"$APPTAINER" run bam-ci.sif --help | head -5 || true
