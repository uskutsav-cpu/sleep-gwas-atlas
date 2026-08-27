#!/usr/bin/env bash
# Harmonize and munge one or more source-verified, EUR hg19 traits.
# Usage: bash scripts/02_munge.sh insomnia bipolar
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/_common.sh

[ "$#" -gt 0 ] || die "usage: bash scripts/02_munge.sh TRAIT [TRAIT ...]"
require_file "$CONFIG"
validate_panel
require_ldsc
mkdir -p data/harmonized data/munged

for trait in "$@"; do
  raw_file=$(trait_field "$trait" raw_file) || die "trait '$trait' is not in $CONFIG"
  source_build=$(trait_field "$trait" build) || die "trait '$trait' has no build in $CONFIG"
  source_status=$(trait_field "$trait" source_status) || die "trait '$trait' has no source status in $CONFIG"
  ancestry=$(trait_field "$trait" ancestry) || die "trait '$trait' has no ancestry in $CONFIG"
  [ "$source_status" = "SOURCE_VERIFIED" ] || die \
    "$trait has source_status=$source_status; verify and register the selected GWAS before munging"
  [ "$ancestry" = "EUR" ] || die \
    "$trait has ancestry=$ancestry; this EUR workflow refuses a different or unresolved ancestry"
  [ "$source_build" = "hg19" ] || [ "$source_build" = "GRCh37" ] || die \
    "$trait has build=$source_build; make an explicit build decision before munging"
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
