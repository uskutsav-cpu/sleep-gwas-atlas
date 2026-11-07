#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
base="$PWD/outputs/shared-genetics-atlas"
inp="$base/multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run"
out="$base/global_audit/mvp_adhd_mdd_frailty_partial"
/opt/miniconda3/bin/python3 - "$inp" "$ref" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]);ref=Path(sys.argv[2]);r=json.load(open(p/'prepared_receipt.json'))
for n in ('mdd','adhd','FI','FFS','HFRS'):
 h=hashlib.sha256()
 with open(p/(n+'.sumstats.gz'),'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 assert h.hexdigest()==r['prepared_sha256'][n],n
for n,h in r['ldscore_sha256'].items():assert hashlib.sha256((ref/(n+'.l2.ldscore.gz')).read_bytes()).hexdigest()==h,n
print('Five matched input files and 22 LD-score hashes verified.')
PY
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$inp/mdd.sumstats.gz,$inp/adhd.sumstats.gz,$inp/FI.sumstats.gz,$inp/FFS.sumstats.gz,$inp/HFRS.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --print-delete-vals --out "$out/ldsc-$tag"
