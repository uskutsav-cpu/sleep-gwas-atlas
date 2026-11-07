#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
DATA='/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT="$PWD/outputs/shared-genetics-atlas/global_audit/insomnia_frailty_same_support/intercept_sensitivity"
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
for int in 0.165 0.193 0; do
 "$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$DATA/HFRS.sumstats.gz,$DATA/S.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --intercept-gencov "0,$int" --out "$OUT/intercept_$int"
done
