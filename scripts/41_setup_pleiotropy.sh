#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

PLACO_URL=https://raw.githubusercontent.com/RayDebashree/PLACO/3ba3cae1d323ad117fb4540e620bcefa79f70663/PLACO_v0.2.0.R
PLACO_PATH=.r-env/share/placo/PLACO_v0.2.0.R
PLACO_SHA=fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124
PLEIOFDR_URL=https://github.com/precimed/pleiofdr.git
PLEIOFDR_COMMIT=0da963cac22fe8de9030166d7aea974bcbb0a367
PLEIOFDR_DIR=work/pleiofdr
REFERENCE_URL=https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
REFERENCE_PATH=ref/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
REFERENCE_BYTES=2383912974
TEMPLATE_URL=https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/9545380.ref
TEMPLATE_PATH=ref/pleiofdr/9545380.ref
TEMPLATE_BYTES=274423819
TEMPLATE_SHA=06268420a0ec04e4529e832e1d4f4a53231b078cc5a3741a3eb215a5a4e1a9d5
DOWNLOAD_SOFTWARE=false
DOWNLOAD_REFERENCE=false
DOWNLOAD_TEMPLATE=false

for argument in "$@"; do
  case "$argument" in
    --download-software) DOWNLOAD_SOFTWARE=true ;;
    --download-reference) DOWNLOAD_REFERENCE=true ;;
    --download-template) DOWNLOAD_TEMPLATE=true ;;
    *) echo "Usage: bash scripts/41_setup_pleiotropy.sh [--download-software] [--download-template] [--download-reference]" >&2; exit 2 ;;
  esac
done

echo "PLACO+ source: 7,515 bytes; SHA256 $PLACO_SHA"
echo "pleioFDR commit: $PLEIOFDR_COMMIT"
echo "pleioFDR reference: $REFERENCE_BYTES bytes (2.22 GiB download)"
echo "pleioFDR variant template: $TEMPLATE_BYTES bytes (0.26 GiB download); SHA256 $TEMPLATE_SHA"

if [ -s "$PLACO_PATH" ]; then
  observed=$(shasum -a 256 "$PLACO_PATH" | awk '{print $1}')
  [ "$observed" = "$PLACO_SHA" ] || { echo "ERROR: installed PLACO+ source checksum differs" >&2; exit 1; }
elif [ "$DOWNLOAD_SOFTWARE" = true ]; then
  mkdir -p "$(dirname "$PLACO_PATH")"
  temporary="${PLACO_PATH}.partial"
  trap 'rm -f "$temporary"' EXIT
  curl -fL --retry 4 --retry-delay 2 "$PLACO_URL" -o "$temporary"
  observed=$(shasum -a 256 "$temporary" | awk '{print $1}')
  [ "$observed" = "$PLACO_SHA" ] || { echo "ERROR: PLACO+ source checksum mismatch" >&2; exit 1; }
  install -m 0444 "$temporary" "$PLACO_PATH"
  rm -f "$temporary"
  trap - EXIT
else
  echo "BLOCKED: PLACO+ source is absent; use --download-software for the small pinned source."
fi

if [ -d "$PLEIOFDR_DIR/.git" ]; then
  observed=$(git -C "$PLEIOFDR_DIR" rev-parse HEAD)
  [ "$observed" = "$PLEIOFDR_COMMIT" ] || { echo "ERROR: pleioFDR checkout differs from the pin" >&2; exit 1; }
elif [ "$DOWNLOAD_SOFTWARE" = true ]; then
  mkdir -p "$(dirname "$PLEIOFDR_DIR")"
  git clone --filter=blob:none "$PLEIOFDR_URL" "$PLEIOFDR_DIR"
  git -C "$PLEIOFDR_DIR" checkout --detach "$PLEIOFDR_COMMIT"
else
  echo "BLOCKED: pleioFDR code is absent; use --download-software for the small pinned checkout."
fi

if [ -s "$REFERENCE_PATH" ]; then
  observed_bytes=$(stat -f %z "$REFERENCE_PATH" 2>/dev/null || stat -c %s "$REFERENCE_PATH")
  [ "$observed_bytes" = "$REFERENCE_BYTES" ] || { echo "ERROR: pleioFDR reference byte count differs" >&2; exit 1; }
  shasum -a 256 "$REFERENCE_PATH"
elif [ "$DOWNLOAD_REFERENCE" = true ]; then
  free_bytes=$(df -Pk "$ROOT" | awk 'NR==2 {printf "%.0f\n", $4 * 1024}')
  minimum_bytes=21474836480
  [ "$free_bytes" -ge "$minimum_bytes" ] || {
    echo "ERROR: reference download requires at least 20 GiB free; observed $free_bytes bytes" >&2
    exit 1
  }
  mkdir -p "$(dirname "$REFERENCE_PATH")"
  curl -fL --retry 4 --retry-delay 5 -C - "$REFERENCE_URL" -o "$REFERENCE_PATH"
  observed_bytes=$(stat -f %z "$REFERENCE_PATH" 2>/dev/null || stat -c %s "$REFERENCE_PATH")
  [ "$observed_bytes" = "$REFERENCE_BYTES" ] || { echo "ERROR: pleioFDR reference byte count differs" >&2; exit 1; }
  shasum -a 256 "$REFERENCE_PATH"
else
  echo "No reference download requested. Re-run with --download-reference only after reviewing the 2.22 GiB plan."
fi

if [ -s "$TEMPLATE_PATH" ]; then
  observed_bytes=$(stat -f %z "$TEMPLATE_PATH" 2>/dev/null || stat -c %s "$TEMPLATE_PATH")
  [ "$observed_bytes" = "$TEMPLATE_BYTES" ] || { echo "ERROR: pleioFDR variant-template byte count differs" >&2; exit 1; }
  observed=$(shasum -a 256 "$TEMPLATE_PATH" | awk '{print $1}')
  [ "$observed" = "$TEMPLATE_SHA" ] || { echo "ERROR: pleioFDR variant-template checksum differs" >&2; exit 1; }
elif [ "$DOWNLOAD_TEMPLATE" = true ]; then
  mkdir -p "$(dirname "$TEMPLATE_PATH")"
  temporary="${TEMPLATE_PATH}.partial"
  trap 'rm -f "$temporary"' EXIT
  curl -fL --retry 4 --retry-delay 5 "$TEMPLATE_URL" -o "$temporary"
  observed_bytes=$(stat -f %z "$temporary" 2>/dev/null || stat -c %s "$temporary")
  [ "$observed_bytes" = "$TEMPLATE_BYTES" ] || { echo "ERROR: pleioFDR variant-template byte count differs" >&2; exit 1; }
  observed=$(shasum -a 256 "$temporary" | awk '{print $1}')
  [ "$observed" = "$TEMPLATE_SHA" ] || { echo "ERROR: pleioFDR variant-template checksum differs" >&2; exit 1; }
  install -m 0444 "$temporary" "$TEMPLATE_PATH"
  rm -f "$temporary"
  trap - EXIT
else
  echo "No variant-template download requested. Re-run with --download-template after reviewing the 0.26 GiB plan."
fi

if [ -f ref/pleiofdr/runtime.provenance.json ]; then
  python3 scripts/pleiotropy_contract.py --verify-runtime
elif [ -s "$PLACO_PATH" ] && [ -d "$PLEIOFDR_DIR/.git" ] && \
     [ -s "$REFERENCE_PATH" ] && [ -s "$TEMPLATE_PATH" ]; then
  python3 scripts/pleiotropy_contract.py --seal-runtime
else
  echo "Pleiotropy runtime remains unsealed until all four pinned components are present."
fi
