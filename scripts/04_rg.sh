#!/usr/bin/env bash
# Run sleep x disease LDSC genetic correlations after the h2 gate.
#
# Usage:
#   bash scripts/04_rg.sh --h2 results/tables/h2_summary.tsv
#   bash scripts/04_rg.sh --pair insomnia mdd  # controlled one-pair rerun
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/_common.sh

require_file "$CONFIG"
require_ldsc
LOGDIR=${LDSC_LOGDIR:-results/logs}
RG_OUT=${RG_OUT:-results/tables/rg_matrix.tsv}
INCLUSION_OUT=${INCLUSION_OUT:-results/tables/phase1_inclusion.tsv}
READINESS_OUT=${READINESS_OUT:-results/tables/trait_readiness.tsv}
mkdir -p "$LOGDIR" "$(dirname "$RG_OUT")"
validate_panel
"$PYTHON_BIN" scripts/00_validate_panel.py --manifest "$CONFIG" --lock "$PANEL_LOCK" \
  --write-provenance results/tables/analysis_panel_provenance.tsv

if [ "${1:-}" = "--pair" ]; then
  [ "$#" -eq 3 ] || die "usage: bash scripts/04_rg.sh --pair SLEEP_TRAIT DISEASE_TRAIT"
  sleep_traits=$2
  disease_traits=$3
  [ "$(trait_field "$sleep_traits" domain)" = "sleep" ] || die \
    "$sleep_traits is not a sleep/circadian trait in the locked panel"
  [ "$(trait_field "$disease_traits" domain)" != "sleep" ] || die \
    "$disease_traits is not a non-sleep trait in the locked panel"
  LOGDIR="$LOGDIR/pair_${sleep_traits}__${disease_traits}"
  mkdir -p "$LOGDIR"
  if [ "$RG_OUT" = "results/tables/rg_matrix.tsv" ]; then
    RG_OUT="results/tables/rg_pair_${sleep_traits}__${disease_traits}.tsv"
  fi
  pair_mode=1
else
  [ "${1:-}" = "--h2" ] && [ "$#" -eq 2 ] || die "usage: bash scripts/04_rg.sh --h2 results/tables/h2_summary.tsv"
  h2_table=$2
  require_file "$h2_table"
  "$PYTHON_BIN" scripts/10_phase0_audit.py --config "$CONFIG" --lock "$PANEL_LOCK" \
    --h2 "$h2_table" --out "$READINESS_OUT" --strict
  "$PYTHON_BIN" scripts/09_select_phase1_traits.py --config "$CONFIG" --lock "$PANEL_LOCK" \
    --readiness "$READINESS_OUT" --h2 "$h2_table" --out "$INCLUSION_OUT"
  sleep_traits=$(awk -F'\t' '
    NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
    $h["domain"] == "sleep" && $h["include_phase1"] == "True" {print $h["trait_id"]}
  ' "$INCLUSION_OUT" | tr '\n' ' ')
  disease_traits=$(awk -F'\t' '
    NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
    $h["domain"] != "sleep" && $h["include_phase1"] == "True" {print $h["trait_id"]}
  ' "$INCLUSION_OUT" | tr '\n' ' ')
  pair_mode=0
fi

[ -n "$sleep_traits" ] || die "no readiness-eligible sleep traits passed h2 QC"
[ -n "$disease_traits" ] || die "no readiness-eligible disease traits passed h2 QC"
for trait in $sleep_traits $disease_traits; do
  require_file "data/munged/$trait.sumstats.gz"
done

list=$(printf '%s\n' $disease_traits | sed 's|^|data/munged/|; s|$|.sumstats.gz|' | paste -sd, -)
for sleep_trait in $sleep_traits; do
  echo "==> rg $sleep_trait vs ${disease_traits}"
  if [ "$pair_mode" -eq 1 ]; then
    rg_out="$LOGDIR/rg_${sleep_trait}__${disease_traits}"
  else
    rg_out="$LOGDIR/rg_$sleep_trait"
  fi
  "$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" \
    --rg "data/munged/$sleep_trait.sumstats.gz,$list" \
    --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
    --out "$rg_out"
done

if [ "$pair_mode" -eq 1 ]; then
  "$PYTHON_BIN" scripts/05_collate.py --mode rg --config "$CONFIG" \
    --logdir "$LOGDIR" --out "$RG_OUT"
  echo "Controlled pair result: $RG_OUT"
else
  "$PYTHON_BIN" scripts/05_collate.py --mode rg --config "$CONFIG" \
    --inclusion "$INCLUSION_OUT" --logdir "$LOGDIR" --out "$RG_OUT"
  "$PYTHON_BIN" scripts/06_heatmap.py --rg "$RG_OUT" --config "$CONFIG" \
    --out results/figures/fig2_rg_heatmap.png
fi
