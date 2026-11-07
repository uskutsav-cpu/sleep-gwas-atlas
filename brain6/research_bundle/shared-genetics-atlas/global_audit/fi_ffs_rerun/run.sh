#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
FRAILTY='/Volumes/Extreme SSD/sleep-frailty-full-tp-v3-20261009/munged_full_native'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT='outputs/shared-genetics-atlas/global_audit/fi_ffs_rerun/fi_ffs'
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
"$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$FRAILTY/frailty.sumstats.gz,$FRAILTY/fried_frailty_score.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$PWD/$OUT"
