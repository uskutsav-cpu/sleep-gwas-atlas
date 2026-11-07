#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
DATA='/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT='outputs/shared-genetics-atlas/global_audit/insomnia_frailty_same_support'
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
for trait in FI FFS HFRS; do
 "$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$DATA/$trait.sumstats.gz,$DATA/S.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$PWD/$OUT/${trait}_insomnia"
done
