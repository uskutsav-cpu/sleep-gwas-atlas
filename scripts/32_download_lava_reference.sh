#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

SOURCE_TABLE=config/lava_reference_sources.tsv
TARGET=ref/lava/ukb_v1.1
MIN_FREE_BYTES=37580963840
DOWNLOAD=false

usage() {
  echo "Usage: bash scripts/32_download_lava_reference.sh [--target DIR] [--download]"
  echo "Without --download this only reports the exact 13.14 GiB archive plan and disk preflight."
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --target) TARGET=$2; shift 2 ;;
    --download) DOWNLOAD=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ ! -f "$SOURCE_TABLE" ]; then
  echo "ERROR: missing source registry: $SOURCE_TABLE" >&2
  exit 1
fi
TOTAL_BYTES=$(awk -F '\t' 'NR>1 {sum += $4} END {printf "%.0f", sum}' "$SOURCE_TABLE")
echo "Registered compressed archives: $TOTAL_BYTES bytes (13.14 GiB)"
echo "Official uncompressed reference: 15 GiB"
mkdir -p "$TARGET"
FREE_KIB=$(df -Pk "$TARGET" | awk 'NR==2 {print $4}')
FREE_BYTES=$((FREE_KIB * 1024))
PRESENT_KIB=$(find "$TARGET" -maxdepth 2 -type f \( -name 'lava-ukb-v1.1_chr*.zip' -o -name 'lava-ukb-v1.1_chr*.zip.part' -o -name 'lava-ukb-v1.1_chr*.info' -o -name 'lava-ukb-v1.1_chr*.bcor' \) -exec du -k {} + 2>/dev/null | awk '{sum += $1} END {printf "%.0f", sum+0}')
EFFECTIVE_BYTES=$((FREE_BYTES + PRESENT_KIB * 1024))
echo "Free space at target: $FREE_BYTES bytes"
echo "State-aware free-plus-verified-target-progress: $EFFECTIVE_BYTES bytes"
echo "Preflight minimum (archives + extraction + margin): $MIN_FREE_BYTES bytes (35 GiB)"

if [ "$DOWNLOAD" != true ]; then
  echo "No download requested; reference storage was not changed."
  exit 0
fi
if [ "$TARGET" != "ref/lava/ukb_v1.1" ]; then
  echo "ERROR: production acquisition must use the policy-locked target ref/lava/ukb_v1.1" >&2
  exit 1
fi
if [ "$EFFECTIVE_BYTES" -lt "$MIN_FREE_BYTES" ]; then
  echo "ERROR: insufficient free space for the checksum-ledgered download and extraction" >&2
  exit 1
fi
if [ -f "$TARGET/reference.provenance.json" ]; then
  python3 scripts/lava_contract.py --verify-reference
  echo "The immutable LAVA reference is already complete; no files were replaced."
  exit 0
fi
if [ -f "$TARGET/download_manifest.tsv" ] || [ -f "$TARGET/extracted_manifest.tsv" ]; then
  echo "ERROR: partial unsealed LAVA manifest exists; inspect it before any replacement" >&2
  exit 1
fi
if ! command -v unzip >/dev/null 2>&1; then
  echo "ERROR: unzip is required" >&2
  exit 1
fi

ARCHIVE_DIR="$TARGET/archives"
mkdir -p "$ARCHIVE_DIR"
MANIFEST_TMP="$TARGET/download_manifest.tsv.tmp"
printf 'archive_id\tchromosomes\turl\tarchive_bytes\tarchive_filename\tsha256\n' > "$MANIFEST_TMP"
tail -n +2 "$SOURCE_TABLE" | while IFS=$'\t' read -r archive_id chromosomes url expected_bytes filename; do
  destination="$ARCHIVE_DIR/$filename"
  partial="$destination.part"
  if [ -f "$destination" ] && [ "$(wc -c < "$destination" | tr -d ' ')" != "$expected_bytes" ]; then
    echo "ERROR: completed archive has the wrong size; refusing replacement: $destination" >&2
    exit 1
  fi
  if [ ! -f "$destination" ]; then
    if [ -f "$partial" ] && [ "$(wc -c < "$partial" | tr -d ' ')" -gt "$expected_bytes" ]; then
      echo "ERROR: partial archive exceeds the registered byte count: $partial" >&2
      exit 1
    fi
    curl -fL --retry 8 --retry-all-errors --continue-at - "$url" -o "$partial"
    actual_bytes=$(wc -c < "$partial" | tr -d ' ')
    if [ "$actual_bytes" != "$expected_bytes" ]; then
      echo "ERROR: byte-count mismatch for $archive_id" >&2
      exit 1
    fi
    mv "$partial" "$destination"
  fi
  archive_sha=$(shasum -a 256 "$destination" | awk '{print $1}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$archive_id" "$chromosomes" "$url" "$expected_bytes" "$filename" "$archive_sha" >> "$MANIFEST_TMP"
  needs_extraction=false
  case "$archive_id" in
    lava_ukb_v1.1_01_02) first=1; last=2 ;;
    lava_ukb_v1.1_03_04) first=3; last=4 ;;
    lava_ukb_v1.1_05_06) first=5; last=6 ;;
    lava_ukb_v1.1_07_09) first=7; last=9 ;;
    lava_ukb_v1.1_10_12) first=10; last=12 ;;
    lava_ukb_v1.1_13_16) first=13; last=16 ;;
    lava_ukb_v1.1_17_23) first=17; last=22 ;;
    *) echo "ERROR: unrecognized registered archive ID: $archive_id" >&2; exit 1 ;;
  esac
  for chromosome in $(seq "$first" "$last"); do
    if [ ! -s "$TARGET/lava-ukb-v1.1_chr${chromosome}.info" ] || [ ! -s "$TARGET/lava-ukb-v1.1_chr${chromosome}.bcor" ]; then
      needs_extraction=true
    fi
  done
  if [ "$needs_extraction" = true ]; then
    unzip -oq "$destination" -d "$TARGET"
  else
    echo "Archive $archive_id is already extracted; existing chromosome payloads were preserved."
  fi
done

EXTRACTED_TMP="$TARGET/extracted_manifest.tsv.tmp"
printf 'chromosome\tfile_type\tpath\tbytes\tsha256\n' > "$EXTRACTED_TMP"
for chromosome in $(seq 1 22); do
  if [ ! -f "$TARGET/lava-ukb-v1.1_chr${chromosome}.info" ] || \
     [ ! -f "$TARGET/lava-ukb-v1.1_chr${chromosome}.bcor" ]; then
    echo "ERROR: exact LAVA .info/.bcor pair is missing for chromosome $chromosome" >&2
    exit 1
  fi
  for suffix in info bcor; do
    extracted="$TARGET/lava-ukb-v1.1_chr${chromosome}.${suffix}"
    extracted_bytes=$(wc -c < "$extracted" | tr -d ' ')
    extracted_sha=$(shasum -a 256 "$extracted" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$chromosome" "$suffix" "lava-ukb-v1.1_chr${chromosome}.${suffix}" "$extracted_bytes" "$extracted_sha" >> "$EXTRACTED_TMP"
  done
done
mv "$MANIFEST_TMP" "$TARGET/download_manifest.tsv"
mv "$EXTRACTED_TMP" "$TARGET/extracted_manifest.tsv"
python3 scripts/lava_contract.py --seal-reference
echo "Downloaded and extracted the complete registered LAVA UKB v1.1 reference."
