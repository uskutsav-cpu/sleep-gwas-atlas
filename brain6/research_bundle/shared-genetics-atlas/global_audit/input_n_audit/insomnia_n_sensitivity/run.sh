#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ADHD='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/adhd.sumstats.gz'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT="$PWD/outputs/shared-genetics-atlas/global_audit/input_n_audit/insomnia_n_sensitivity"
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
for SCENARIO in constantN rawN; do
 "$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$OUT/matched_${SCENARIO}.sumstats.gz,$ADHD" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$OUT/rg_${SCENARIO}" > "$OUT/${SCENARIO}_console.log" 2>&1
done
