#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
reference='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
data='/Volumes/Extreme SSD/brain6-physical-factor-20261101'
out="$(cd "$(dirname "$0")" && pwd)"
export MPLCONFIGDIR="$out" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for a in walking grip activity; do
 "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$data/$a.common.sumstats.gz,$data/questionnaire.common.sumstats.gz" --ref-ld-chr "$reference" --w-ld-chr "$reference" --print-delete-vals --out "$out/${a}_questionnaire"
done
