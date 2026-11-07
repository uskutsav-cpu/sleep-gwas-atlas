#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
reference='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
data='/Volumes/Brain6Work/physical-function-source-sensitivity-20261028'
out="$PWD/outputs/shared-genetics-atlas/component_exclusion/physical_function_ldsc/insomnia_source_sensitivity"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for physical in walking grip; do
  for sleep in questionnaire registry; do
    "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$data/$physical.sumstats.gz,$data/$sleep.sumstats.gz" --ref-ld-chr "$reference" --w-ld-chr "$reference" --print-delete-vals --out "$out/${physical}_${sleep}"
  done
done
