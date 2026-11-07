#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
reference='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
data='/Volumes/Brain6Work/ukb-hfrs-comparison-20261028'
out="$PWD/outputs/shared-genetics-atlas/global_audit/ukb_hfrs_public"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for hfrs in ukb finngen; do
  if [ "$hfrs" = ukb ]; then h='/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz'; else h="$data/finngen_hfrs.sumstats.gz"; fi
  for sleep in questionnaire registry; do
    "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$h,$data/$sleep.sumstats.gz" --ref-ld-chr "$reference" --w-ld-chr "$reference" --print-delete-vals --out "$out/${hfrs}_${sleep}"
  done
done
