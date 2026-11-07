#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
MUNGED='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged'
FRAILTY='/Volumes/Extreme SSD/sleep-frailty-full-tp-v3-20261009/munged_full_native'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT="$PWD/outputs/shared-genetics-atlas/global_audit/brain6_frailty_covariance"
export MPLCONFIGDIR="$PWD/work/ldsc_mpl" TMPDIR='/Volumes/Brain6Work/tmp'
for spec in 'insomnia,frailty' 'insomnia,fried_frailty_score' 'adhd,frailty'; do
  a="${spec%,*}"; b="${spec#*,}"
  "$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$MUNGED/$a.sumstats.gz,$FRAILTY/$b.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$OUT/${a}_${b}"
done
