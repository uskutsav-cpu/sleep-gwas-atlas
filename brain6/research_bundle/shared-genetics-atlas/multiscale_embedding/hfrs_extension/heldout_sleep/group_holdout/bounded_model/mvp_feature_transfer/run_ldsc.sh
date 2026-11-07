#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
m='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
base="$PWD/outputs/shared-genetics-atlas"
out="$base/multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer"
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
/usr/bin/shasum -a 256 "$m/adhd.sumstats.gz" "$m/mdd.sumstats.gz" "$m/scz.sumstats.gz" "$m/bipolar.sumstats.gz" "$m/alz.sumstats.gz" "$m/parkinson.sumstats.gz" > "$out/input-hashes-$tag.txt"
export STREAM_RECEIPT="$out/stream-$tag.json"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg <(/opt/miniconda3/bin/python3 "$base/replication/mvp_finngen_genomewide/stream.py"),"$m/adhd.sumstats.gz","$m/mdd.sumstats.gz","$m/scz.sumstats.gz","$m/bipolar.sumstats.gz","$m/alz.sumstats.gz","$m/parkinson.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/ldsc-$tag"
test -s "$STREAM_RECEIPT"
/opt/miniconda3/bin/python3 - "$STREAM_RECEIPT" "$base/replication/mvp_finngen_genomewide/stream_receipt.json" <<'PY'
import json,sys
assert json.load(open(sys.argv[1]))==json.load(open(sys.argv[2])), 'MVP stream mismatch'
print('MVP full-stream receipt verified')
PY
