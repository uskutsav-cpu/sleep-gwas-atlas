#!/usr/bin/env bash
# Run the isolated 12 x extension-primary-pass LDSC rg family.
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/_common.sh

[ "$#" -eq 2 ] && [ "$1" = "--h2" ] || die \
  "usage: bash discovery_extension/scripts/14_rg_extension.sh --h2 EXTENSION_H2_TABLE"
h2_table=$2
require_file "$h2_table"
python3 discovery_extension/scripts/00_verify_core_checkpoint.py
python3 discovery_extension/scripts/05_validate_extension_panel.py
require_ldsc

extension_traits=$(awk -F'\t' '
  NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
  $h["primary_rg_eligibility"] == "PRIMARY_PASS" {print $h["extension_trait_id"]}
' "$h2_table")
[ -n "$extension_traits" ] || die "no extension traits passed the rerun h2 gate"
sleep_traits=$(awk -F'\t' '
  NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
  $h["domain"] == "sleep" {print $h["trait_id"]}
' config/analysis_panel.tsv)
[ "$(printf '%s\n' $sleep_traits | wc -l | tr -d ' ')" -eq 12 ] || die "core sleep family is not 12 traits"

list=""
for trait in $extension_traits; do
  file="discovery_extension/data/munged/$trait.sumstats.gz"
  require_file "$file"
  if [ -z "$list" ]; then list=$file; else list="$list,$file"; fi
done

logdir=discovery_extension/logs/rg
mkdir -p "$logdir" discovery_extension/results/ldsc
for sleep in $sleep_traits; do
  core_file="data/munged/$sleep.sumstats.gz"
  require_file "$core_file"
  echo "==> extension rg $sleep vs rerun-h2-pass panel"
  "$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" \
    --rg "$core_file,$list" \
    --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
    --out "$logdir/rg_$sleep"
done

python3 discovery_extension/scripts/13_collate_extension_ldsc.py \
  --mode rg --logdir "$logdir" --h2 "$h2_table" \
  --out discovery_extension/results/ldsc/extension_rg_primary.tsv \
  --pair-universe-out discovery_extension/results/ldsc/extension_pair_universe.tsv
