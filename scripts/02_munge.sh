#!/usr/bin/env bash
# Harmonize + munge one or more traits.  Usage: bash scripts/02_munge.sh insomnia mdd
set -euo pipefail
cd "$(dirname "$0")/.."
CFG=config/traits.tsv
for T in "$@"; do
  echo "==> $T : harmonize"
  RAW=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $6}' $CFG)
  [ -n "$RAW" ] || { echo "  !! $T not in $CFG"; exit 1; }
  [ -f "data/raw/$RAW" ] || { echo "  !! missing data/raw/$RAW"; exit 1; }
  python3 scripts/01_harmonize.py --trait "$T" --config $CFG \
      --infile "data/raw/$RAW" --outdir data/harmonized

  echo "==> $T : munge"
  # --merge-alleles restricts to HapMap3 and fixes allele orientation. [LDSC]
  conda run -n ldsc python ldsc/munge_sumstats.py \
      --sumstats data/harmonized/$T.harmonized.tsv.gz \
      --merge-alleles ref/w_hm3.snplist \
      --chunksize 500000 \
      --out data/munged/$T
done
echo "Munged files in data/munged/. Check each .log for 'Writing summary statistics'."
