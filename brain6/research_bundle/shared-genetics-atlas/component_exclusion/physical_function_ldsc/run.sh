#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
reference='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
insomnia='/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz'
out="$PWD/outputs/shared-genetics-atlas/component_exclusion/physical_function_ldsc"
export TMPDIR='/Volumes/Brain6Work/tmp'
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for trait in walking grip; do
  case "$trait" in
    walking) source='/Volumes/Brain6Work/neale_924_walking.ldsc.bgz';;
    grip) source='/Volumes/Brain6Work/neale_46irnt_grip.ldsc.bgz';;
  esac
  "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$source,$insomnia" --ref-ld-chr "$reference" --w-ld-chr "$reference" --out "$out/$trait"
done
