#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
base="$PWD/outputs/shared-genetics-atlas"
inp="$base/multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run"
out="$base/global_audit/mvp_adhd_frailty_partial"
/opt/miniconda3/bin/python3 - "$inp" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]);r=json.load(open(p/'prepared_receipt.json'))
for n in ('adhd','FI','FFS','HFRS'):
 h=hashlib.sha256()
 with open(p/(n+'.sumstats.gz'),'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 assert h.hexdigest()==r['prepared_sha256'][n]
print('Four input hashes verified.')
PY
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$inp/adhd.sumstats.gz,$inp/FI.sumstats.gz,$inp/FFS.sumstats.gz,$inp/HFRS.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --print-delete-vals --out "$out/ldsc-$tag"
