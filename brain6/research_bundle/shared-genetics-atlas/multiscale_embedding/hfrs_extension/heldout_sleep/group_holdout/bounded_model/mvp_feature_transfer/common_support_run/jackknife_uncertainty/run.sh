#!/bin/bash
set -euo pipefail
root='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
ref='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
out="$PWD/outputs/shared-genetics-atlas/multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run"
/opt/miniconda3/bin/python3 - "$out" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]);r=json.load(open(p/'prepared_receipt.json'))
for n,h in r['prepared_sha256'].items():
 x=hashlib.sha256()
 with open(p/(n+'.sumstats.gz'),'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):x.update(b)
 assert x.hexdigest()==h,n+' hash mismatch'
print('Prepared file hashes verified.')
PY
tag="$(date -u +%Y%m%dT%H%M%S)-$$"
jack="$out/jackknife_uncertainty"
export TMPDIR='/Volumes/Brain6Work/tmp' MPLCONFIGDIR="$PWD/work/ldsc_mpl" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$root/.ldsc-env/bin/python3.9" "$root/ldsc/ldsc.py" --rg "$out/MVP.sumstats.gz,$out/adhd.sumstats.gz,$out/mdd.sumstats.gz,$out/scz.sumstats.gz,$out/bipolar.sumstats.gz,$out/alz.sumstats.gz,$out/parkinson.sumstats.gz,$out/FI.sumstats.gz,$out/FFS.sumstats.gz,$out/HFRS.sumstats.gz" --ref-ld-chr "$ref" --w-ld-chr "$ref" --print-delete-vals --out "$jack/ldsc-$tag"
