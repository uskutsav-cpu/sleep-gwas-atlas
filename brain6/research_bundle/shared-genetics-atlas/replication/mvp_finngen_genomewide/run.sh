#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
adhd='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz'
out="$PWD/outputs/shared-genetics-atlas/replication/mvp_finngen_genomewide"
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
export STREAM_RECEIPT="$out/stream_receipt-$tag.json"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$out/stream.py"),"$adhd" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/mvp_adhd-$tag"
test -s "$STREAM_RECEIPT"
