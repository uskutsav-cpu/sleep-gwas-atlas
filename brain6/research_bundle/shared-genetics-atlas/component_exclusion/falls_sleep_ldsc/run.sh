#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
falls='/Volumes/Brain6NHANES/neale-falls-20261028/falls.ldsc.bgz'
out="$PWD/outputs/shared-genetics-atlas/component_exclusion/falls_sleep_ldsc"
export TMPDIR='/Volumes/Brain6NHANES/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for trait in questionnaire registry; do
 case "$trait" in
 questionnaire) sleep='/Volumes/Brain6Work/brain6-frailty-paired/insomnia.sumstats.gz';;
 registry) sleep='/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz';;
 esac
 "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$falls,$sleep" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/$trait"
done
