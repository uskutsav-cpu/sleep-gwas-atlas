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
prepare_analysis_workspace

for trait in "$@"; do
  raw_file=$(trait_field "$trait" raw_file) || die "trait '$trait' is not in $CONFIG"
  source_id=$(trait_field "$trait" source_id) || die "trait '$trait' has no source_id in $CONFIG"
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

  harmonize_infile="data/raw/$raw_file"
  prefilter_args=()
  prefilter_strategy=$(hm3_prefilter_field "$trait" strategy 2>/dev/null || true)
  if [ -n "$prefilter_strategy" ]; then
    [ "$prefilter_strategy" = "HAPMAP3_RSID_ALLOWLIST" ] || die \
      "$trait has unsupported prefilter strategy=$prefilter_strategy"
    planned_source=$(hm3_prefilter_field "$trait" source_id)
    snp_column=$(hm3_prefilter_field "$trait" snp_column)
    allowlist_path=$(hm3_prefilter_field "$trait" allowlist_path)
    [ "$planned_source" = "$source_id" ] || die \
      "$trait prefilter source=$planned_source does not match selected source=$source_id"
    require_file "$allowlist_path"
    source_bytes=$(public_source_field "$source_id" archive_bytes)
    source_sha256=$(public_source_field "$source_id" archive_sha256)
    [ -n "$source_bytes" ] && [ -n "$source_sha256" ] || die \
      "$trait prefilter source lacks registered bytes/SHA-256"
    mkdir -p "$HARMONIZED_DIR/.prefilter"
    harmonize_infile="$HARMONIZED_DIR/.prefilter/$trait.hm3.tsv.gz"
    prefilter_provenance="$HARMONIZED_DIR/$trait.prefilter.provenance.json"
    echo "==> $trait: streaming registered HapMap3 prefilter"
    "$PYTHON_BIN" scripts/21_prefilter_hm3.py \
      --input "data/raw/$raw_file" --output "$harmonize_infile" \
      --provenance "$prefilter_provenance" --snp-column "$snp_column" \
      --allowlist "$allowlist_path" --expected-input-bytes "$source_bytes" \
      --expected-input-sha256 "$source_sha256"
    prefilter_args=(--prefilter-provenance "$prefilter_provenance")
  fi

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
    ${prefilter_args[@]+"${prefilter_args[@]}"} \
    --infile "$harmonize_infile" --outdir "$HARMONIZED_DIR"

  echo "==> $trait: HapMap3 munging"
  ldsc_ignore_args=()
  if grep -Fq 'FRQ column absent - source-level MAF QC must be documented' \
      "$HARMONIZED_DIR/$trait.qc.txt"; then
    # The harmonized contract carries an explicit FRQ=NA placeholder when the
    # source has no allele-frequency field.  LDSC otherwise recognizes FRQ as
    # a required numeric column and drops every row before HapMap3 matching.
    # Ignore it only when the harmonization ledger proves the source lacked it;
    # a present-but-invalid frequency column must still fail visibly.
    ldsc_ignore_args=(--ignore FRQ)
    echo "  source has no FRQ; ignoring the all-missing FRQ placeholder in LDSC"
  fi
  "$LDSC_PYTHON" "$LDSC_DIR/munge_sumstats.py" \
    --sumstats "$HARMONIZED_DIR/$trait.harmonized.tsv.gz" \
    --merge-alleles "$REF_DIR/w_hm3.snplist" --chunksize 500000 \
    ${ldsc_ignore_args[@]+"${ldsc_ignore_args[@]}"} \
    --out "$MUNGED_DIR/$trait"
  require_file "$MUNGED_DIR/$trait.sumstats.gz"
  require_file "$MUNGED_DIR/$trait.log"
  if grep -q 'WARNING' "$MUNGED_DIR/$trait.log"; then
    echo "  WARNING: LDSC reported warnings for $trait; inspect $MUNGED_DIR/$trait.log before h2." >&2
  fi
done

echo "Munging complete. Review each $HARMONIZED_DIR/*.qc.txt and $MUNGED_DIR/*.log before h2."
