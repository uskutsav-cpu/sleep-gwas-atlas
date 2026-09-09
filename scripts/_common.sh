#!/usr/bin/env bash
# Shared helpers for the executable Phase 0/1 shell entry points.
# This file is sourced; do not invoke it directly.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}
LDSC_PYTHON=${LDSC_PYTHON:-.ldsc-env/bin/python}
LDSC_DIR=${LDSC_DIR:-ldsc}
# Production entry points always use the locked atlas-v1.0 manifest. Keeping
# this path non-overridable prevents an environment variable from silently
# changing the analysed phenotype set.
CONFIG=config/analysis_panel.tsv
PANEL_LOCK=config/analysis_panel.lock.json
VARIANT_MAPPINGS=config/variant_mapping_plans.tsv
LIFTOVER_PLANS=config/liftover_plans.tsv
HM3_PREFILTER_PLANS=config/hm3_prefilter_plans.tsv
PUBLIC_SOURCES=config/public_gwas_sources.tsv

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
      for (i = 1; i <= NF; ++i) if ($i == "trait_id") trait_column = i
      if (!column || !trait_column) exit 2
      next
    }
    $trait_column == trait_id { print $column; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$CONFIG"
}

variant_mapping_field() {
  local trait=$1
  local column=$2
  [ -f "$VARIANT_MAPPINGS" ] || return 1
  awk -F'\t' -v trait_id="$trait" -v wanted="$column" '
    NR == 1 {
      for (i = 1; i <= NF; ++i) if ($i == wanted) column = i
      for (i = 1; i <= NF; ++i) if ($i == "trait_id") trait_column = i
      if (!column || !trait_column) exit 2
      next
    }
    $trait_column == trait_id { print $column; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$VARIANT_MAPPINGS"
}

liftover_field() {
  local trait=$1
  local column=$2
  [ -f "$LIFTOVER_PLANS" ] || return 1
  awk -F'\t' -v trait_id="$trait" -v wanted="$column" '
    NR == 1 {
      for (i = 1; i <= NF; ++i) if ($i == wanted) column = i
      for (i = 1; i <= NF; ++i) if ($i == "trait_id") trait_column = i
      if (!column || !trait_column) exit 2
      next
    }
    $trait_column == trait_id { print $column; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$LIFTOVER_PLANS"
}

hm3_prefilter_field() {
  local trait=$1
  local column=$2
  [ -f "$HM3_PREFILTER_PLANS" ] || return 1
  awk -F'\t' -v trait_id="$trait" -v wanted="$column" '
    NR == 1 {
      for (i = 1; i <= NF; ++i) if ($i == wanted) column = i
      for (i = 1; i <= NF; ++i) if ($i == "trait_id") trait_column = i
      if (!column || !trait_column) exit 2
      next
    }
    $trait_column == trait_id { print $column; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$HM3_PREFILTER_PLANS"
}

public_source_field() {
  local source=$1
  local column=$2
  [ -f "$PUBLIC_SOURCES" ] || return 1
  awk -F'\t' -v source_id="$source" -v wanted="$column" '
    NR == 1 {
      for (i = 1; i <= NF; ++i) if ($i == wanted) column = i
      for (i = 1; i <= NF; ++i) if ($i == "source_id") source_column = i
      if (!column || !source_column) exit 2
      next
    }
    $source_column == source_id { print $column; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$PUBLIC_SOURCES"
}

require_ldsc() {
  require_file "$LDSC_DIR/ldsc.py"
  require_file "$LDSC_DIR/munge_sumstats.py"
  require_file "ref/w_hm3.snplist"
  require_file "ref/eur_w_ld_chr/1.l2.ldscore.gz"
}

validate_panel() {
  "$PYTHON_BIN" scripts/00_validate_panel.py --manifest "$CONFIG" --lock "$PANEL_LOCK"
}
