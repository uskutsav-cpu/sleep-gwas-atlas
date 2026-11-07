#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
inp='/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity'
base="$PWD/outputs/shared-genetics-atlas/replication"
out="$base/mvp_frailty_genomewide"
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
/usr/bin/shasum -a 256 "$inp/FI.sumstats.gz" "$inp/FFS.sumstats.gz" "$inp/HFRS.sumstats.gz" > "$out/input-hashes-$tag.txt"
/opt/miniconda3/bin/python3 - "$inp/support.json" "$out/input-hashes-$tag.txt" <<'PY'
import json,sys
expected=json.load(open(sys.argv[1]))['output_sha256']
rows=open(sys.argv[2]).read().splitlines()
for trait,row in zip(('FI','FFS','HFRS'),rows):
    assert row.split()[0]==expected[trait], f'{trait} hash mismatch'
PY
export STREAM_RECEIPT="$out/stream-$tag.json"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$base/mvp_finngen_genomewide/stream.py"),"$inp/FI.sumstats.gz","$inp/FFS.sumstats.gz","$inp/HFRS.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/ldsc-$tag"
test -s "$STREAM_RECEIPT"
