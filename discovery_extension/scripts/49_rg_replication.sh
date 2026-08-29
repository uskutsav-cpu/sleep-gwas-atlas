#!/usr/bin/env bash
# Run only h2-eligible pairs from the fixed 217-pair independent-replication family.
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/_common.sh

python3 discovery_extension/scripts/00_verify_core_checkpoint.py
python3 discovery_extension/scripts/22_lock_replication_manifest.py --validate-only
require_ldsc
python3 discovery_extension/scripts/49_prepare_replication_rg_jobs.py

logdir=discovery_extension/logs/replication/rg
mkdir -p "$logdir"
while IFS=$'\t' read -r sleep source_ids pair_ids source_count; do
  [ "$sleep" != "sleep_trait" ] || continue
  core="data/munged/$sleep.sumstats.gz"
  require_file "$core"
  list=$core
  old_ifs=$IFS
  IFS=','
  for source_id in $source_ids; do
    sumstats="discovery_extension/data/replication/munged/$source_id.sumstats.gz"
    require_file "$sumstats"
    list="$list,$sumstats"
  done
  IFS=$old_ifs
  echo "==> replication rg $sleep vs $source_count h2-pass sources"
  "$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" \
    --rg "$list" \
    --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
    --out "$logdir/rg_replication_$sleep"
done < discovery_extension/results/replication/replication_rg_jobs.tsv
python3 discovery_extension/scripts/49_collate_replication_rg.py --logdir "$logdir"
