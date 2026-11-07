#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
base="$PWD/outputs/shared-genetics-atlas"
out="$base/component_exclusion/mvp_older_weakness"
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
export STREAM_TAG="$tag" STREAM_RECEIPT="$out/mvp_receipt-$tag.json"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$out/weakness_stream.py"),<(/opt/miniconda3/bin/python3 "$base/replication/mvp_finngen_genomewide/stream.py") --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/ldsc-$tag"
test -s "$STREAM_RECEIPT"
test -s "$out/receipt_$tag.json"
/opt/miniconda3/bin/python3 - "$STREAM_RECEIPT" "$base/replication/mvp_finngen_genomewide/stream_receipt.json" "$out/receipt_$tag.json" "$base/component_exclusion/older_weakness/receipt_retry_registry.json" <<'PY'
import json,sys
assert json.load(open(sys.argv[1]))==json.load(open(sys.argv[2])), 'MVP receipt mismatch'
a=json.load(open(sys.argv[3])); b=json.load(open(sys.argv[4])); assert a==b, 'Weakness receipt mismatch'
print('Both complete-source receipts match previous audits.')
PY
