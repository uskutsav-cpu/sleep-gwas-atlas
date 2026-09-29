#!/usr/bin/env bash
# Run LDSC h2 and collate the QC gate for selected munged traits.
#
# By default binary traits run on the liability scale only when the config
# contains a cited population prevalence. Use --observed-scale only for the
# interim h2 power gate / Phase 1 rg work; it is not a substitute for a final
# liability-scale h2 table.
#
# Usage:
#   bash scripts/03_h2_qc.sh insomnia bipolar
#   bash scripts/03_h2_qc.sh --observed-scale insomnia bipolar
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/_common.sh

observed_scale=0
if [ "${1:-}" = "--observed-scale" ]; then
  observed_scale=1
  shift
fi
[ "$#" -gt 0 ] || die "usage: bash scripts/03_h2_qc.sh [--observed-scale] TRAIT [TRAIT ...]"
require_file "$CONFIG"
validate_panel
require_ldsc
prepare_analysis_workspace

LOGDIR=${LDSC_LOGDIR:-results/logs}
H2_OUT=${H2_OUT:-results/tables/h2_summary.tsv}
mkdir -p "$LOGDIR" "$(dirname "$H2_OUT")"

for trait in "$@"; do
  trait_type=$(trait_field "$trait" type) || die "trait '$trait' is not in $CONFIG"
  source_status=$(trait_field "$trait" source_status) || die "trait '$trait' has no source status in $CONFIG"
  [ "$source_status" = "SOURCE_VERIFIED" ] || die \
    "$trait has source_status=$source_status; complete source verification before h2"
  require_file "$MUNGED_DIR/$trait.sumstats.gz"
  echo "==> h2 $trait ($trait_type)"

  command=("$LDSC_PYTHON" "$LDSC_DIR/ldsc.py"
    --h2 "$MUNGED_DIR/$trait.sumstats.gz"
    --ref-ld-chr "$REF_DIR/eur_w_ld_chr/" --w-ld-chr "$REF_DIR/eur_w_ld_chr/"
    --out "$LOGDIR/h2_$trait")

  if [ "$trait_type" = "binary" ] && [ "$observed_scale" -eq 0 ]; then
    population_prevalence=$(trait_field "$trait" pop_prev)
    ncase=$(trait_field "$trait" ncase)
    ncontrol=$(trait_field "$trait" ncontrol)
    prevalence_citation=$(trait_field "$trait" pop_prev_citation 2>/dev/null || true)
    [ -n "$prevalence_citation" ] && [ "$prevalence_citation" != "UNRESOLVED" ] || die \
      "$trait has no cited pop_prev_citation. Refusing an unsourced liability-scale h2; rerun with --observed-scale only for interim Phase 1 QC."
    sample_prevalence=$("$PYTHON_BIN" -c "print(round(float('$ncase') / (float('$ncase') + float('$ncontrol')), 6))")
    command+=(--samp-prev "$sample_prevalence" --pop-prev "$population_prevalence")
  elif [ "$trait_type" = "binary" ]; then
    echo "  NOTE: observed-scale h2 selected; do not report it as liability-scale h2." >&2
  fi
  "${command[@]}"
done

# The canonical table is cumulative.  Restricting this collation to "$@"
# would silently replace earlier valid rows whenever h2 is run incrementally.
# The loop above already fails if any newly requested LDSC run fails, so now
# rebuild the table from every locked-panel h2 log present in LOGDIR.
"$PYTHON_BIN" scripts/05_collate.py --mode h2 --config "$CONFIG" \
  --logdir "$LOGDIR" --out "$H2_OUT"
