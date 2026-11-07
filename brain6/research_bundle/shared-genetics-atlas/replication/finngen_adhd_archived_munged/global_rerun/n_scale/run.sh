#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
INSOMNIA='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/insomnia.sumstats.gz'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT="$PWD/outputs/shared-genetics-atlas/replication/finngen_adhd_archived_munged/global_rerun/n_scale"
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
"$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --h2 "$OUT/synthetic_totalN.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$OUT/adhd_h2_totalN" > "$OUT/h2_console.log" 2>&1
"$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$INSOMNIA,$OUT/synthetic_totalN.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$OUT/rg_totalN" > "$OUT/rg_console.log" 2>&1
