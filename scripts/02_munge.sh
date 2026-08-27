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
  mapping_strategy=$(variant_mapping_field "$trait" strategy 2>/dev/null || true)
  mapping_path=$(variant_mapping_field "$trait" map_path 2>/dev/null || true)
  mapping_bytes=$(variant_mapping_field "$trait" map_bytes 2>/dev/null || true)
  mapping_sha256=$(variant_mapping_field "$trait" map_sha256 2>/dev/null || true)
  liftover_strategy=$(liftover_field "$trait" strategy 2>/dev/null || true)
  liftover_path=$(liftover_field "$trait" chain_path 2>/dev/null || true)
  liftover_bytes=$(liftover_field "$trait" chain_bytes 2>/dev/null || true)
  liftover_sha256=$(liftover_field "$trait" chain_sha256 2>/dev/null || true)
  [ "$source_status" = "SOURCE_VERIFIED" ] || die \
    "$trait has source_status=$source_status; verify and register the selected GWAS before munging"
  [ "$ancestry" = "EUR" ] || die \
    "$trait has ancestry=$ancestry; this EUR workflow refuses a different or unresolved ancestry"
  if [ "$liftover_strategy" = "UCSC_CHAIN_POINT" ]; then
    [ "$source_build" = "hg38" ] || [ "$source_build" = "GRCh38" ] || die \
      "$trait has build=$source_build; its registered liftover plan requires hg38/GRCh38"
  elif [ "$mapping_strategy" = "BY_RSID_ALLELES" ]; then
    [ "$source_build" = "UNRESOLVED" ] || [ "$source_build" = "hg19" ] || [ "$source_build" = "GRCh37" ] || die \
      "$trait has build=$source_build; an rsID map cannot excuse a conflicting source build"
  else
    [ "$source_build" = "hg19" ] || [ "$source_build" = "GRCh37" ] || die \
      "$trait has build=$source_build; make an explicit build decision before munging"
  fi
  require_file "data/raw/$raw_file"

  source_build_args=()
  if [ "$source_build" = "hg19" ] || [ "$source_build" = "GRCh37" ] || \
     [ "$source_build" = "hg38" ] || [ "$source_build" = "GRCh38" ]; then
    source_build_args=(--source-build "$source_build")
  fi
  variant_map_args=()
  if [ -n "$mapping_strategy" ] || [ -n "$mapping_path" ]; then
    [ -n "$mapping_strategy" ] && [ -n "$mapping_path" ] && \
      [ -n "$mapping_bytes" ] && [ -n "$mapping_sha256" ] || die \
      "$trait has an incomplete variant-mapping plan"
    require_file "$mapping_path"
    require_file "$mapping_path.provenance.json"
    variant_map_args=(
      --variant-map "$mapping_path"
      --variant-map-strategy "$mapping_strategy"
      --expected-variant-map-bytes "$mapping_bytes"
      --expected-variant-map-sha256 "$mapping_sha256"
    )
  fi
  liftover_args=()
  if [ -n "$liftover_strategy" ] || [ -n "$liftover_path" ]; then
    [ "$liftover_strategy" = "UCSC_CHAIN_POINT" ] && [ -n "$liftover_path" ] && \
      [ -n "$liftover_bytes" ] && [ -n "$liftover_sha256" ] || die \
      "$trait has an incomplete or unsupported liftover plan"
    require_file "$liftover_path"
    liftover_args=(
      --liftover-chain "$liftover_path"
      --expected-liftover-chain-bytes "$liftover_bytes"
      --expected-liftover-chain-sha256 "$liftover_sha256"
    )
  fi

  echo "==> $trait: Phase 0 harmonization"
  "$PYTHON_BIN" scripts/01_harmonize.py \
    --trait "$trait" --config "$CONFIG" \
    ${source_build_args[@]+"${source_build_args[@]}"} \
    ${variant_map_args[@]+"${variant_map_args[@]}"} \
    ${liftover_args[@]+"${liftover_args[@]}"} \
    --infile "data/raw/$raw_file" --outdir data/harmonized

  echo "==> $trait: HapMap3 munging"
  ldsc_ignore_args=()
  if grep -Fq 'FRQ column absent - source-level MAF QC must be documented' \
      "data/harmonized/$trait.qc.txt"; then
    # The harmonized contract carries an explicit FRQ=NA placeholder when the
    # source has no allele-frequency field.  LDSC otherwise recognizes FRQ as
    # a required numeric column and drops every row before HapMap3 matching.
    # Ignore it only when the harmonization ledger proves the source lacked it;
    # a present-but-invalid frequency column must still fail visibly.
    ldsc_ignore_args=(--ignore FRQ)
    echo "  source has no FRQ; ignoring the all-missing FRQ placeholder in LDSC"
  fi
  "$LDSC_PYTHON" "$LDSC_DIR/munge_sumstats.py" \
    --sumstats "data/harmonized/$trait.harmonized.tsv.gz" \
    --merge-alleles ref/w_hm3.snplist --chunksize 500000 \
    ${ldsc_ignore_args[@]+"${ldsc_ignore_args[@]}"} \
    --out "data/munged/$trait"
  require_file "data/munged/$trait.sumstats.gz"
  require_file "data/munged/$trait.log"
  if grep -q 'WARNING' "data/munged/$trait.log"; then
    echo "  WARNING: LDSC reported warnings for $trait; inspect data/munged/$trait.log before h2." >&2
  fi
done

echo "Munging complete. Review each data/harmonized/*.qc.txt and data/munged/*.log before h2."
