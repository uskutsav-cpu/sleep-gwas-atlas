#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
DATA='/Volumes/Brain6Work/brain6-frailty-paired'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT="$PWD/outputs/shared-genetics-atlas/global_audit/brain6_frailty_covariance/three_frailty_common_support/insomnia_n_sensitivity"
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
for TRAIT in FI FFS HFRS; do
 "$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$DATA/$TRAIT.sumstats.gz,$OUT/insomnia.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --print-delete-vals --out "$OUT/$TRAIT" > "$OUT/${TRAIT}_console.log" 2>&1
done
