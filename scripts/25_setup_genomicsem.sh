#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BIN="${CONDA_BIN:-conda}"
R_ENV="${R_ENV:-$ROOT_DIR/.r-env}"
GENOMICSEM_COMMIT="6b65ca5db39fdade08b0d811477be1cdd57b5039"
SIMSALAPAR_URL="https://cran.r-project.org/src/contrib/simsalapar_1.0-13.tar.gz"

"$CONDA_BIN" env create --prefix "$R_ENV" --file "$ROOT_DIR/environment/genomicsem.yml"
"$R_ENV/bin/Rscript" -e \
  "install.packages('$SIMSALAPAR_URL', repos=NULL, type='source')"

SOURCE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sleep-atlas-genomicsem.XXXXXX")"
trap 'rm -rf "$SOURCE_DIR"' EXIT
git clone --filter=blob:none https://github.com/GenomicSEM/GenomicSEM.git "$SOURCE_DIR"
git -C "$SOURCE_DIR" checkout --detach "$GENOMICSEM_COMMIT"
test "$(git -C "$SOURCE_DIR" rev-parse HEAD)" = "$GENOMICSEM_COMMIT"
"$R_ENV/bin/R" CMD INSTALL "$SOURCE_DIR"

"$R_ENV/bin/Rscript" -e \
  'stopifnot(getRversion() == "4.3.3", as.character(packageVersion("GenomicSEM")) == "0.0.5", as.character(packageVersion("simsalapar")) == "1.0-13")'
