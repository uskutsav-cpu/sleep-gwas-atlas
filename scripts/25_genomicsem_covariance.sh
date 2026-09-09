#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${R_BIN:-$ROOT_DIR/.r-env/bin/Rscript}"

test -x "$R_BIN" || {
  echo "ERROR: pinned R runtime missing at $R_BIN; run scripts/25_setup_genomicsem.sh" >&2
  exit 1
}

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

cd "$ROOT_DIR"
exec "$R_BIN" scripts/25_genomicsem_covariance.R "$@"
