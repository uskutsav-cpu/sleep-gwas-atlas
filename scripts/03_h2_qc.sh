#!/usr/bin/env bash
# SNP heritability + the QC gate that decides which traits survive.
# Usage: bash scripts/03_h2_qc.sh insomnia mdd chronotype ...
set -euo pipefail
cd "$(dirname "$0")/.."
CFG=config/traits.tsv
mkdir -p results/logs
for T in "$@"; do
  TYPE=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $4}' $CFG)
  PREV=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $10}' $CFG)
  NCASE=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $7}' $CFG)
  NCON=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $8}' $CFG)
  echo "==> h2 $T ($TYPE)"
  if [ "$TYPE" = "binary" ]; then
    # Liability scale: sample prevalence from the effective-N case fraction,
    # population prevalence from config. Both are YOUR assumptions - be able
    # to defend the population prevalence number in your meeting.
    SAMPPREV=$(python3 -c "print(round($NCASE/($NCASE+$NCON),6))")
    conda run -n ldsc python ldsc/ldsc.py \
        --h2 data/munged/$T.sumstats.gz \
        --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
        --samp-prev $SAMPPREV --pop-prev $PREV \
        --out results/logs/h2_$T
  else
    conda run -n ldsc python ldsc/ldsc.py \
        --h2 data/munged/$T.sumstats.gz \
        --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
        --out results/logs/h2_$T
  fi
done
python3 scripts/05_collate.py --mode h2 --logdir results/logs --out results/tables/h2_summary.tsv
