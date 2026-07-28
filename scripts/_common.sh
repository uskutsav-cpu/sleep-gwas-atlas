#!/usr/bin/env bash
# Shared helpers for the executable Phase 0/1 shell entry points.
# This file is sourced; do not invoke it directly.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}
LDSC_PYTHON=${LDSC_PYTHON:-$PYTHON_BIN}
LDSC_DIR=${LDSC_DIR:-ldsc}
CONFIG=${CONFIG:-config/traits.tsv}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_file() {
  [ -f "$1" ] || die "required file not found: $1"
}

trait_field() {
  local trait=$1
  local column=$2
  awk -F'\t' -v trait_id="$trait" -v wanted="$column" '
    NR == 1 {
      for (i = 1; i <= NF; ++i) if ($i == wanted) column = i
      if (!column) exit 2
      next
    }
    $1 == trait_id { print $column; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$CONFIG"
}

require_ldsc() {
  require_file "$LDSC_DIR/ldsc.py"
  require_file "$LDSC_DIR/munge_sumstats.py"
  require_file "ref/w_hm3.snplist"
  require_file "ref/eur_w_ld_chr/1.l2.ldscore.gz"
}
