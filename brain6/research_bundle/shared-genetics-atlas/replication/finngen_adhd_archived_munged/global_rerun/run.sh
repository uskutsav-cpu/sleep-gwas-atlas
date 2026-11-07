#!/bin/bash
set -euo pipefail
ROOT='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas'
MUNGED='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged'
REF='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr/'
OUT='outputs/shared-genetics-atlas/replication/finngen_adhd_archived_munged/global_rerun/insomnia_adhd'
export MPLCONFIGDIR="$PWD/work/ldsc_mpl"
"$ROOT/.ldsc-env/bin/python3.9" "$ROOT/ldsc/ldsc.py" --rg "$MUNGED/insomnia.sumstats.gz,$MUNGED/track_b_finngen_r13_F5_ADHD.sumstats.gz" --ref-ld-chr "$REF" --w-ld-chr "$REF" --out "$PWD/$OUT"
