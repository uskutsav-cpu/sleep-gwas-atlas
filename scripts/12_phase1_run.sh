#!/usr/bin/env bash
# Phase 1 driver: harmonized sumstats -> munged -> h2 -> QC gate.
#
# Assumes scripts/01_harmonize.py has already produced
# data/harmonized/<trait>.harmonized.tsv.gz. Run 09_fetch + 01_harmonize first;
# raw files that were evicted are re-fetched by 09_fetch from the recorded URL,
# and the SHA-256 in config/public_gwas_sources.tsv proves you got the same
# bytes back.
#
#   bash scripts/12_phase1_run.sh bmi alz mdd t2d
#
# LDSC lives behind symlinks (ldsc, .ldsc-env, ref) pointing at a shared
# install. They are gitignored; see docs for how to recreate them.
set -uo pipefail
cd "$(dirname "$0")/.."

PY=.ldsc-env/bin/python
CFG=config/traits.tsv
mkdir -p data/munged results/logs results/tables

for T in "$@"; do
  H="data/harmonized/$T.harmonized.tsv.gz"
  if [ ! -f "$H" ]; then
    echo "!! $T: no harmonized file ($H) -- run 01_harmonize.py first"
    continue
  fi

  if [ ! -f "data/munged/$T.sumstats.gz" ]; then
    echo "==> $T : munge"
    $PY ldsc/munge_sumstats.py --sumstats "$H" \
        --merge-alleles ref/w_hm3.snplist --chunksize 500000 \
        --out "data/munged/$T" > "results/logs/munge_$T.log" 2>&1 \
      || { echo "   !! munge failed, see results/logs/munge_$T.log"; continue; }
  fi

  echo "==> $T : h2"
  # Liability scale needs BOTH prevalences. Every pop_prev in this registry is
  # UNKNOWN by design (rule 1: uncited beats invented), so h2 is computed on
  # the OBSERVED scale and the liability conversion is deferred until the
  # mentor supplies sourced prevalences. Z = h2/SE is invariant to that choice,
  # so the QC gate below is unaffected.
  TYPE=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $4}' $CFG)
  PREV=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $10}' $CFG)
  if [ "$TYPE" = "binary" ] && [ "$PREV" != "UNKNOWN" ] && [ "$PREV" != "NA" ]; then
    NCASE=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $7}' $CFG)
    NCON=$(awk -F'\t' -v t="$T" 'NR>1 && $1==t {print $8}' $CFG)
    SAMPPREV=$($PY -c "print(round($NCASE/($NCASE+$NCON),6))")
    $PY ldsc/ldsc.py --h2 "data/munged/$T.sumstats.gz" \
        --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
        --samp-prev "$SAMPPREV" --pop-prev "$PREV" \
        --out "results/logs/h2_$T" > /dev/null 2>&1
  else
    $PY ldsc/ldsc.py --h2 "data/munged/$T.sumstats.gz" \
        --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
        --out "results/logs/h2_$T" > /dev/null 2>&1
  fi
  grep -E 'Total (Observed|Liability) scale h2|^Intercept' "results/logs/h2_$T.log" \
    | sed "s/^/   $T  /"
done

echo
echo "==> collate + QC gate"
$PY scripts/05_collate.py --mode h2 --logdir results/logs \
    --out results/tables/h2_summary.tsv
