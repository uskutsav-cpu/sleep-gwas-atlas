#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$PWD}"
STORAGE_PATH="${2:-$ROOT/data/raw}"
if [[ ! -e "$STORAGE_PATH" ]]; then STORAGE_PATH="$ROOT"; fi
available_kib="$(df -Pk "$STORAGE_PATH" | awk 'NR == 2 {print $4}')"
if [[ ! "$available_kib" =~ ^[0-9]+$ ]]; then
  echo "ERROR: could not determine available storage for $STORAGE_PATH" >&2
  exit 2
fi

required_kib=$((20 * 1024 * 1024))
available_gib=$((available_kib / 1024 / 1024))
if (( available_kib < required_kib )); then
  echo "ERROR: collection requires at least 20 GiB free; found ${available_gib} GiB at $STORAGE_PATH" >&2
  exit 75
fi

echo "Storage gate passed: ${available_gib} GiB free at $STORAGE_PATH"
