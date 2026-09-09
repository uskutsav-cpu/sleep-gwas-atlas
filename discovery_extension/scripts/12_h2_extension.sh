#!/usr/bin/env bash
# Run observed-scale LDSC h2 for the locked extension and collate the primary gate.
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/_common.sh

[ "$#" -gt 0 ] || die "usage: bash discovery_extension/scripts/12_h2_extension.sh --all|TRAIT [TRAIT ...]"
python3 discovery_extension/scripts/05_validate_extension_panel.py
require_ldsc

all_mode=0
if [ "$1" = "--all" ]; then
  [ "$#" -eq 1 ] || die "--all cannot be combined with trait IDs"
  all_mode=1
  traits=$(awk -F'\t' 'NR > 1 {print $1}' discovery_extension/config/candidate_traits.tsv)
else
  traits="$*"
fi

logdir=discovery_extension/logs/h2
mkdir -p "$logdir" discovery_extension/results/ldsc
for trait in $traits; do
  awk -F'\t' -v wanted="$trait" 'NR > 1 && $1 == wanted {found=1} END {exit !found}' \
    discovery_extension/config/candidate_traits.tsv || die "$trait is not in the locked extension panel"
  sumstats="discovery_extension/data/munged/$trait.sumstats.gz"
  require_file "$sumstats"
  echo "==> extension observed-scale h2 $trait"
  "$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" \
    --h2 "$sumstats" \
    --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
    --out "$logdir/h2_$trait"
done

if [ "$all_mode" -eq 1 ]; then
  python3 discovery_extension/scripts/13_collate_extension_ldsc.py \
    --mode h2 --logdir "$logdir" \
    --out discovery_extension/results/ldsc/extension_trait_readiness.tsv \
    --h2-failed-out discovery_extension/results/ldsc/extension_trait_qc_failed.tsv
else
  echo "Controlled partial h2 run complete; the locked 100-trait h2 table was not regenerated."
fi
