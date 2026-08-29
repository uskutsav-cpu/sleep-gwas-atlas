#!/usr/bin/env bash
# Run the source-level h2 gate before any replication pair test.
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/_common.sh

python3 discovery_extension/scripts/00_verify_core_checkpoint.py
python3 discovery_extension/scripts/22_lock_replication_manifest.py --validate-only
require_ldsc

logdir=discovery_extension/logs/replication/h2
mkdir -p "$logdir" discovery_extension/results/replication
sources=$(awk -F'\t' '
  NR == 1 {for (i=1; i<=NF; i++) h[$i]=i; next}
  $h["source_curation_status"] == "COMPLETE_BEFORE_RESULTS" && !seen[$h["replication_source_id"]]++ {print $h["replication_source_id"]}
' discovery_extension/config/replication_manifest.tsv)
[ "$(printf '%s\n' $sources | wc -l | tr -d ' ')" -eq 13 ] || die "replication manifest is not the exact 13-source family"
for source_id in $sources; do
  sumstats="discovery_extension/data/replication/munged/$source_id.sumstats.gz"
  require_file "$sumstats"
  echo "==> replication h2 $source_id"
  "$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" \
    --h2 "$sumstats" \
    --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ \
    --out "$logdir/h2_$source_id"
done
python3 discovery_extension/scripts/48_collate_replication_h2.py --logdir "$logdir"
