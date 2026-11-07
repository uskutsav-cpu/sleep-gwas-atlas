#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
falls='/Volumes/Brain6NHANES/neale-falls-20261028/falls.ldsc.bgz'
sleep='/Volumes/Brain6Work/physical-function-source-sensitivity-20261028'
out="$PWD/outputs/shared-genetics-atlas/component_exclusion/falls_sleep_ldsc/paired_source_sensitivity"
export TMPDIR='/Volumes/Brain6NHANES/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for trait in questionnaire registry; do
 "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$falls,$sleep/$trait.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --print-delete-vals --out "$out/$trait"
done
