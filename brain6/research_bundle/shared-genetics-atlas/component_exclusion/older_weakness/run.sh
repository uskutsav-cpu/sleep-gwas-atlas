#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
inputs='/Volumes/Brain6Work/physical-function-source-sensitivity-20261028'
out="$PWD/outputs/shared-genetics-atlas/component_exclusion/older_weakness"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for t in questionnaire registry; do
  export STREAM_TAG="retry_$t"
  "$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$out/stream.py"),"$inputs/$t.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/retry_$t"
  test -s "$out/receipt_retry_$t.json"
done
