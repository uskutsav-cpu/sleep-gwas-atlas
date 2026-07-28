#!/usr/bin/env bash
# Harmonize and munge one or more curated, EUR hg19 traits.
# Usage: bash scripts/02_munge.sh insomnia mdd
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/_common.sh

[ "$#" -gt 0 ] || die "usage: bash scripts/02_munge.sh TRAIT [TRAIT ...]"
require_file "$CONFIG"
require_ldsc
mkdir -p data/harmonized data/munged

for trait in "$@"; do
  raw_file=$(trait_field "$trait" raw_file) || die "trait '$trait' is not in $CONFIG"
  source_build=$(trait_field "$trait" build) || die "trait '$trait' has no build in $CONFIG"
  curation_status=$(trait_field "$trait" status) || die "trait '$trait' has no curation status in $CONFIG"
  [ "$curation_status" = "CURATED" ] || die "$trait has status=$curation_status; complete and audit Phase 0 curation before munging"
  require_file "data/raw/$raw_file"

  echo "==> $trait: Phase 0 harmonization"
  "$PYTHON_BIN" scripts/01_harmonize.py \
    --trait "$trait" --config "$CONFIG" --source-build "$source_build" \
    --infile "data/raw/$raw_file" --outdir data/harmonized

  echo "==> $trait: HapMap3 munging"
  "$LDSC_PYTHON" "$LDSC_DIR/munge_sumstats.py" \
    --sumstats "data/harmonized/$trait.harmonized.tsv.gz" \
    --merge-alleles ref/w_hm3.snplist --chunksize 500000 \
    --out "data/munged/$trait"
  require_file "data/munged/$trait.sumstats.gz"
  require_file "data/munged/$trait.log"
  if grep -q 'WARNING' "data/munged/$trait.log"; then
    echo "  WARNING: LDSC reported warnings for $trait; inspect data/munged/$trait.log before h2." >&2
  fi
done

echo "Munging complete. Review each data/harmonized/*.qc.txt and data/munged/*.log before h2."
