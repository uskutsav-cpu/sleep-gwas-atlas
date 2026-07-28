#!/usr/bin/env bash
# Pairwise genetic correlation: every sleep trait x every disease trait.
# Usage: bash scripts/04_rg.sh                (uses all traits with status=PASS)
#        bash scripts/04_rg.sh insomnia mdd   (explicit pair)
set -euo pipefail
cd "$(dirname "$0")/.."
CFG=config/traits.tsv
mkdir -p results/logs
if [ "$#" -eq 2 ]; then
  SLEEP="$1"; DIS="$2"
else
  SLEEP=$(awk -F'\t' 'NR>1 && $3=="sleep"  && $12=="PASS" {print $1}' $CFG | tr '\n' ' ')
  DIS=$(  awk -F'\t' 'NR>1 && $3!="sleep"  && $12=="PASS" {print $1}' $CFG | tr '\n' ' ')
fi
echo "sleep traits : $SLEEP"
echo "disease traits: $DIS"
for S in $SLEEP; do
  # LDSC takes a comma-separated list: one sleep trait vs many diseases in
  # a single call, which is far faster than looping pair by pair.
  LIST=$(echo $DIS | tr ' ' '\n' | sed 's|^|data/munged/|; s|$|.sumstats.gz|' | paste -sd, -)
  echo "==> rg $S vs all"
  conda run -n ldsc python ldsc/ldsc.py \
      --rg data/munged/$S.sumstats.gz,$LIST \
      --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
      --out results/logs/rg_$S
done
python3 scripts/05_collate.py --mode rg --logdir results/logs --out results/tables/rg_matrix.tsv
python3 scripts/06_heatmap.py --rg results/tables/rg_matrix.tsv --config $CFG \
    --out results/figures/fig2_rg_heatmap.png
