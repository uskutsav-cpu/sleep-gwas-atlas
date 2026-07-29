#!/usr/bin/env bash
# Download and materialize only source-registered, public GWAS releases.
#
# Usage from repo root:
#   bash scripts/11_materialize_public_gwas.sh --download dashti_2019_sleep_duration
#   bash scripts/11_materialize_public_gwas.sh --materialize dashti_2019_sleep_duration
#   bash scripts/11_materialize_public_gwas.sh --verify
#
# This intentionally does not infer a source, accept a form, perform a
# liftover, or set any trait to CURATED. The Phase 0 audit remains the gate.
set -euo pipefail
cd "$(dirname "$0")/.."

SOURCES=${SOURCES:-config/public_gwas_sources.tsv}
RAW_DIR=${RAW_DIR:-data/raw}
ARCHIVE_DIR=${ARCHIVE_DIR:-$RAW_DIR/.archives}
PYTHON_BIN=${PYTHON_BIN:-python3}
MODE=${1:---verify}
SOURCE_ID=${2:-}

usage() {
  cat <<'EOF'
Usage:
  bash scripts/11_materialize_public_gwas.sh --verify
  bash scripts/11_materialize_public_gwas.sh --download SOURCE_ID
  bash scripts/11_materialize_public_gwas.sh --materialize SOURCE_ID

`--download` only permits rows marked PUBLIC. `--materialize` streams the
registered archive member to gzip and creates hard-linked alias files when the
registry intentionally names the same exact phenotype twice. A source marked
``GZIP_WRAPPED_ZIP_COLUMNS`` is a reviewed multi-phenotype nested archive; it
is handled only by its named, source-specific materializer. A source marked
``DIRECT_TSV`` is a single, uncompressed tabular release and is losslessly
gzip-wrapped after its registered integrity check. A source marked
``DIRECT_GZIP`` is already a single gzipped tabular release and is hard-linked
into the raw directory after its registered integrity check. Large public
Google Drive releases with a reviewed source-specific materializer are fetched
by validated byte ranges rather than a blind resume.
EOF
}

case "$MODE" in
  --verify) [ -z "$SOURCE_ID" ] || { usage >&2; exit 2; } ;;
  --download|--materialize) [ -n "$SOURCE_ID" ] || { usage >&2; exit 2; } ;;
  *) usage >&2; exit 2 ;;
esac

[ -f "$SOURCES" ] || { echo "ERROR: source registry not found: $SOURCES" >&2; exit 1; }
mkdir -p "$RAW_DIR" "$ARCHIVE_DIR"

read_source() {
  awk -F '\t' -v target="$SOURCE_ID" '
    NR == 1 { next }
    $1 == target { print; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$SOURCES"
}

if [ "$MODE" = --verify ]; then
  awk -F '\t' 'NR == 1 { print "source_id\tarchive_present\tarchive_sha256\traw_files_present"; next }
    {
      cmd = "test -s data/raw/.archives/" $6
      archive = (system(cmd) == 0 ? "yes" : "no")
      sha = "UNAVAILABLE"
      if (archive == "yes") {
        cmd = "shasum -a 256 data/raw/.archives/" $6 " | awk '\''{print $1}'\''"
        cmd | getline sha
        close(cmd)
      }
      n = split($10, raw, ","); all = "yes"
      for (i = 1; i <= n; i++) if (system("test -s data/raw/" raw[i]) != 0) all = "no"
      print $1 "\t" archive "\t" sha "\t" all
    }' "$SOURCES"
  exit 0
fi

SOURCE_ROW=$(read_source) || { echo "ERROR: unknown source_id: $SOURCE_ID" >&2; exit 1; }
IFS=$'\t' read -r source_id trait_ids source_page download_url access archive_name archive_bytes archive_sha256 archive_member raw_files pmid ancestry build_status acquisition_status notes <<< "$SOURCE_ROW"
[ "$access" = PUBLIC ] || { echo "ERROR: $source_id is not approved as PUBLIC" >&2; exit 1; }
ARCHIVE="$ARCHIVE_DIR/$archive_name"

require_archive_integrity() {
  local actual_bytes actual_sha256
  actual_bytes=$(wc -c < "$ARCHIVE" | tr -d ' ')
  if [ "$actual_bytes" != "$archive_bytes" ]; then
    echo "ERROR: archive byte count mismatch for $source_id: expected $archive_bytes, got $actual_bytes" >&2
    exit 1
  fi
  actual_sha256=$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')
  if [ "$actual_sha256" != "$archive_sha256" ]; then
    echo "ERROR: archive SHA-256 mismatch for $source_id" >&2
    exit 1
  fi
}

if [ "$MODE" = --download ]; then
  if [ "$archive_member" = "BCAC_2020_META_RSID" ]; then
    "$PYTHON_BIN" scripts/14_ranged_download.py \
      --url "$download_url" --out "$ARCHIVE" \
      --expected-bytes "$archive_bytes" --expected-sha256 "$archive_sha256" \
      --max-chunks 4
  else
    curl --fail --location --retry 3 --continue-at - --output "$ARCHIVE" "$download_url"
  fi
  # A byte-range resume can leave a superficially plausible archive after a
  # proxy interruption. Validate the registered container before it can be
  # streamed into a raw input; retain a failed archive for inspection/re-download.
  if [ "$archive_member" = "GZIP_WRAPPED_ZIP_COLUMNS" ]; then
    "$PYTHON_BIN" scripts/12_materialize_accelerometer_sleep.py --source "$ARCHIVE" --verify-only
  elif [ "$archive_member" = "NEALE_GRIP_RSID_JOIN" ]; then
    NEALE_VARIANTS="$ARCHIVE_DIR/neale_ukbb_round2_variants.tsv.bgz"
    curl --fail --location --retry 3 --continue-at - --output "$NEALE_VARIANTS" \
      "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
    "$PYTHON_BIN" scripts/13_materialize_neale_grip.py \
      --results "$ARCHIVE" --variants "$NEALE_VARIANTS" --verify-only
  elif [ "$archive_member" = "BCAC_2020_META_RSID" ]; then
    "$PYTHON_BIN" scripts/15_materialize_bcac_breast.py --source "$ARCHIVE" --verify-only
  elif [ "$archive_member" = "PHELAN_2017_OVARIAN_OVERALL_RSID" ]; then
    "$PYTHON_BIN" scripts/16_materialize_phelan_ovarian.py \
      --source "$ARCHIVE" --expected-sha256 "$archive_sha256" --verify-only
  elif [ "$archive_member" = "DIRECT_TSV" ]; then
    test -s "$ARCHIVE" || { echo "ERROR: downloaded direct TSV is empty" >&2; exit 1; }
  elif [ "$archive_member" = "DIRECT_GZIP" ]; then
    gzip -t "$ARCHIVE"
  else
    unzip -t "$ARCHIVE" > /dev/null
  fi
  require_archive_integrity
  echo "Downloaded $source_id -> $ARCHIVE"
  exit 0
fi

test -s "$ARCHIVE" || { echo "ERROR: archive missing; run --download $source_id first" >&2; exit 1; }
require_archive_integrity
LOCK_DIR="$RAW_DIR/.${source_id}.materialize.lock"
if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  echo "ERROR: $source_id is already being materialized; wait for it to finish" >&2
  exit 1
fi
TEMP=""
cleanup_materialization() {
  [ -z "$TEMP" ] || rm -f "$TEMP"
  rmdir "$LOCK_DIR" 2>/dev/null || true
}
trap cleanup_materialization EXIT
IFS=',' read -r -a outputs <<< "$raw_files"
if [ "$archive_member" = "BCAC_2020_META_RSID" ]; then
  [ "${#outputs[@]}" -eq 1 ] || {
    echo "ERROR: BCAC_2020_META_RSID must register exactly one raw output" >&2
    exit 1
  }
  PRIMARY="$RAW_DIR/${outputs[0]}"
  "$PYTHON_BIN" scripts/15_materialize_bcac_breast.py --source "$ARCHIVE" --verify-only
  if [ ! -s "$PRIMARY" ]; then
    "$PYTHON_BIN" scripts/15_materialize_bcac_breast.py --source "$ARCHIVE" --out "$PRIMARY"
  fi
  gzip -t "$PRIMARY"
  echo "Materialized $source_id"
  printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
  printf '  raw sha256: '; shasum -a 256 "$PRIMARY" | awk '{print $1}'
  printf '  files: %s\n' "$raw_files"
  exit 0
fi
if [ "$archive_member" = "PHELAN_2017_OVARIAN_OVERALL_RSID" ]; then
  [ "${#outputs[@]}" -eq 1 ] || {
    echo "ERROR: PHELAN_2017_OVARIAN_OVERALL_RSID must register exactly one raw output" >&2
    exit 1
  }
  PRIMARY="$RAW_DIR/${outputs[0]}"
  "$PYTHON_BIN" scripts/16_materialize_phelan_ovarian.py \
    --source "$ARCHIVE" --expected-sha256 "$archive_sha256" --verify-only
  if [ ! -s "$PRIMARY" ]; then
    "$PYTHON_BIN" scripts/16_materialize_phelan_ovarian.py \
      --source "$ARCHIVE" --expected-sha256 "$archive_sha256" --out "$PRIMARY"
  fi
  gzip -t "$PRIMARY"
  echo "Materialized $source_id"
  printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
  printf '  raw sha256: '; shasum -a 256 "$PRIMARY" | awk '{print $1}'
  printf '  files: %s\n' "$raw_files"
  exit 0
fi
if [ "$archive_member" = "NEALE_GRIP_RSID_JOIN" ]; then
  [ "${#outputs[@]}" -eq 1 ] || {
    echo "ERROR: NEALE_GRIP_RSID_JOIN must register exactly one raw output" >&2
    exit 1
  }
  PRIMARY="$RAW_DIR/${outputs[0]}"
  "$PYTHON_BIN" scripts/13_materialize_neale_grip.py \
    --results "$ARCHIVE" --variants "$ARCHIVE_DIR/neale_ukbb_round2_variants.tsv.bgz" --verify-only
  if [ ! -s "$PRIMARY" ]; then
    "$PYTHON_BIN" scripts/13_materialize_neale_grip.py \
      --results "$ARCHIVE" --variants "$ARCHIVE_DIR/neale_ukbb_round2_variants.tsv.bgz" \
      --out "$PRIMARY"
  fi
  gzip -t "$PRIMARY"
  echo "Materialized $source_id"
  printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
  printf '  raw sha256: '; shasum -a 256 "$PRIMARY" | awk '{print $1}'
  printf '  files: %s\n' "$raw_files"
  exit 0
fi
if [ "$archive_member" = "GZIP_WRAPPED_ZIP_COLUMNS" ]; then
  existing_outputs=0
  for output in "${outputs[@]}"; do
    if [ -e "$RAW_DIR/$output" ]; then
      test -s "$RAW_DIR/$output" || {
        echo "ERROR: existing materialized output is empty: $RAW_DIR/$output" >&2
        exit 1
      }
      existing_outputs=$((existing_outputs + 1))
    fi
  done
  if [ "$existing_outputs" -eq 0 ]; then
    "$PYTHON_BIN" scripts/12_materialize_accelerometer_sleep.py --source "$ARCHIVE" --out-dir "$RAW_DIR"
  elif [ "$existing_outputs" -ne "${#outputs[@]}" ]; then
    echo "ERROR: only some outputs already exist for $source_id; inspect before re-materializing" >&2
    exit 1
  fi
  for output in "${outputs[@]}"; do
    test -s "$RAW_DIR/$output" || {
      echo "ERROR: source-specific materializer did not produce $output" >&2
      exit 1
    }
  done
  echo "Materialized $source_id"
  printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
  printf '  files: %s\n' "$raw_files"
  exit 0
fi
if [ "$archive_member" = "DIRECT_TSV" ]; then
  [ "${#outputs[@]}" -eq 1 ] || {
    echo "ERROR: DIRECT_TSV must register exactly one raw output" >&2
    exit 1
  }
  PRIMARY="$RAW_DIR/${outputs[0]}"
  if [ ! -s "$PRIMARY" ]; then
    TEMP="$RAW_DIR/.${outputs[0]}.partial.$$"
    gzip -c "$ARCHIVE" > "$TEMP"
    gzip -t "$TEMP"
    mv "$TEMP" "$PRIMARY"
    TEMP=""
  fi
  echo "Materialized $source_id"
  printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
  printf '  raw sha256: '; shasum -a 256 "$PRIMARY" | awk '{print $1}'
  printf '  files: %s\n' "$raw_files"
  exit 0
fi
if [ "$archive_member" = "DIRECT_GZIP" ]; then
  [ "${#outputs[@]}" -eq 1 ] || {
    echo "ERROR: DIRECT_GZIP must register exactly one raw output" >&2
    exit 1
  }
  PRIMARY="$RAW_DIR/${outputs[0]}"
  if [ ! -e "$PRIMARY" ]; then
    ln "$ARCHIVE" "$PRIMARY"
  elif ! cmp -s "$ARCHIVE" "$PRIMARY"; then
    echo "ERROR: existing raw file differs from registered direct gzip: $PRIMARY" >&2
    exit 1
  fi
  gzip -t "$PRIMARY"
  echo "Materialized $source_id"
  printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
  printf '  raw sha256: '; shasum -a 256 "$PRIMARY" | awk '{print $1}'
  printf '  files: %s\n' "$raw_files"
  exit 0
fi
PRIMARY="$RAW_DIR/${outputs[0]}"
if [ ! -s "$PRIMARY" ]; then
  TEMP="$RAW_DIR/.${outputs[0]}.partial.$$"
  unzip -p "$ARCHIVE" "$archive_member" | gzip -c > "$TEMP"
  gzip -t "$TEMP"
  mv "$TEMP" "$PRIMARY"
  TEMP=""
fi

for output in "${outputs[@]:1}"; do
  TARGET="$RAW_DIR/$output"
  if [ ! -e "$TARGET" ]; then
    ln "$PRIMARY" "$TARGET"
  elif ! cmp -s "$PRIMARY" "$TARGET"; then
    echo "ERROR: existing alias differs from $PRIMARY: $TARGET" >&2
    exit 1
  fi
done

echo "Materialized $source_id"
printf '  archive sha256: '; shasum -a 256 "$ARCHIVE" | awk '{print $1}'
printf '  raw sha256: '; shasum -a 256 "$PRIMARY" | awk '{print $1}'
printf '  files: %s\n' "$raw_files"
