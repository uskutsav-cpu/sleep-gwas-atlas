#!/usr/bin/env bash
# Complete the locked sleep x disease family without promoting h2-QC failures.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/_common.sh

H2=${1:-results/tables/h2_summary.tsv}
PRIMARY=${PRIMARY_RG:-results/tables/rg_primary_phase1.tsv}
SENSITIVITY=${SENSITIVITY_RG:-results/tables/rg_qc_failed_sensitivity.tsv}
FINAL=${FINAL_RG:-results/tables/rg_matrix.tsv}
LOGDIR=${SENSITIVITY_LOGDIR:-results/logs/rg_qc_failed_sensitivity}
READINESS=results/tables/trait_readiness.tsv
INCLUSION=results/tables/phase1_inclusion.tsv

require_file "$H2"
require_file results/tables/rg_matrix.tsv
validate_panel
require_ldsc
prepare_analysis_workspace
"$PYTHON_BIN" scripts/10_phase0_audit.py --config "$CONFIG" --lock "$PANEL_LOCK" \
  --h2 "$H2" --out "$READINESS" --strict
"$PYTHON_BIN" scripts/09_select_phase1_traits.py --config "$CONFIG" --lock "$PANEL_LOCK" \
  --readiness "$READINESS" --h2 "$H2" --out "$INCLUSION"

sleep_traits=$(awk -F'\t' '
  NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
  $h["domain"] == "sleep" {print $h["trait_id"]}
' "$INCLUSION" | tr '\n' ' ')
failed_diseases=$(awk -F'\t' '
  NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
  $h["domain"] != "sleep" && $h["include_phase1"] == "False" && \
    $h["ldsc_ready"] == "True" {print $h["trait_id"]}
' "$INCLUSION" | tr '\n' ' ')
[ -n "$sleep_traits" ] || die "locked panel contains no sleep traits"
[ -n "$failed_diseases" ] || die "no LDSC-ready non-sleep h2-QC failures to analyze"
for trait in $sleep_traits $failed_diseases; do
  require_file "$MUNGED_DIR/$trait.sumstats.gz"
done

mkdir -p "$LOGDIR" "$(dirname "$SENSITIVITY")"
cp results/tables/rg_matrix.tsv "$PRIMARY"
list=$(printf '%s\n' $failed_diseases | sed "s|^|$MUNGED_DIR/|; s|$|.sumstats.gz|" | paste -sd, -)
for sleep_trait in $sleep_traits; do
  echo "==> sensitivity rg $sleep_trait vs $failed_diseases"
  "$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" \
    --rg "$MUNGED_DIR/$sleep_trait.sumstats.gz,$list" \
    --ref-ld-chr "$REF_DIR/eur_w_ld_chr/" --w-ld-chr "$REF_DIR/eur_w_ld_chr/" \
    --out "$LOGDIR/rg_$sleep_trait"
done

"$PYTHON_BIN" scripts/05_collate.py --mode rg --config "$CONFIG" \
  --logdir "$LOGDIR" --out "$SENSITIVITY"
"$PYTHON_BIN" scripts/23_merge_rg_families.py --config "$CONFIG" --h2 "$H2" \
  --primary "$PRIMARY" --sensitivity "$SENSITIVITY" --out "$FINAL"
