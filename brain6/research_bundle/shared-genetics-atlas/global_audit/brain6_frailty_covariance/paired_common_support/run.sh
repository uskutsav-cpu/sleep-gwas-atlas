#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
DATA='/Volumes/Brain6Work/brain6-frailty-paired'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT="$PWD/outputs/shared-genetics-atlas/global_audit/brain6_frailty_covariance/paired_common_support"
export MPLCONFIGDIR="$PWD/work/ldsc_mpl" TMPDIR='/Volumes/Brain6Work/tmp'
for trait in FI FFS; do
 "$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$DATA/$trait.sumstats.gz,$DATA/insomnia.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --print-delete-vals --out "$OUT/$trait"
done
