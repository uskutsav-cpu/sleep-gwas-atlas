#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
inputs='/Volumes/Brain6Work/physical-function-source-sensitivity-20261028'
out="$PWD/outputs/shared-genetics-atlas/component_exclusion/sarcopenia_public_candidate/ldsc_feasibility"
export TMPDIR='/Volumes/Brain6NHANES/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for t in questionnaire registry walking grip; do
 export STREAM_TAG="$t"
 "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$out/stream_source.py"),"$inputs/$t.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/run2_$t"
 test -s "$out/stream_receipt_$t.json"
done
