#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
reference='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
data='/Volumes/Brain6Work/ukb-hfrs-comparison-20261028'
source='/Volumes/Extreme SSD/brain6-ukb-hfrs-without-dementia-20261028.sumstats.gz'
out="$PWD/outputs/shared-genetics-atlas/global_audit/ukb_hfrs_public/dementia_sensitivity"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for sleep in questionnaire registry; do
  "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$source,$data/$sleep.sumstats.gz" --ref-ld-chr "$reference" --w-ld-chr "$reference" --print-delete-vals --out "$out/without_${sleep}"
done
