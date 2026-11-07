#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
walking='/Volumes/Brain6Work/neale_924_walking.ldsc.bgz'
grip='/Volumes/Brain6Work/neale_46irnt_grip.ldsc.bgz'
base="$PWD/outputs/shared-genetics-atlas"
out="$base/component_exclusion/mvp_physical_components"
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
/opt/miniconda3/bin/python3 - "$base/component_exclusion/physical_function_ldsc/input_receipt.json" "$walking" "$grip" <<'PY'
import hashlib,json,sys
a=json.load(open(sys.argv[1]))
for k,p in zip(('walking','grip'),sys.argv[2:]):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 assert h.hexdigest()==a[k]['sha256'], k+' hash mismatch'
print('Both physical-input hashes verified.')
PY
export STREAM_RECEIPT="$out/stream-$tag.json"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$base/replication/mvp_finngen_genomewide/stream.py"),"$walking","$grip" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/ldsc-$tag"
test -s "$STREAM_RECEIPT"
/opt/miniconda3/bin/python3 - "$STREAM_RECEIPT" "$base/replication/mvp_finngen_genomewide/stream_receipt.json" <<'PY'
import json,sys
assert json.load(open(sys.argv[1]))==json.load(open(sys.argv[2])), 'MVP receipt mismatch'
print('MVP full-stream receipt verified.')
PY
