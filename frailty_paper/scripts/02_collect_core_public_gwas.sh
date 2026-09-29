#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
cd "$ROOT"
bash frailty_paper/scripts/09_require_storage.sh "$ROOT"
[[ -f scripts/11_materialize_public_gwas.sh ]] || { echo "Missing existing repo downloader" >&2; exit 1; }

# These are already checksum-pinned in your repository. Download only if the raw file is absent.
SOURCES=(
  atkins_2021_frailty_index
  timmers_2019_parental_lifespan
  deelen_2019_longevity_90th
  neale_2018_left_grip_strength
  bellenguez_2022_alzheimer_stage1
  nalls_2019_parkinson_public_proxy
)

for src in "${SOURCES[@]}"; do
  row=$(awk -F '\t' -v s="$src" 'NR>1 && $1==s {print; exit}' config/public_gwas_sources.tsv)
  [[ -n "$row" ]] || { echo "WARNING: source not registered: $src" >&2; continue; }
  raw=$(printf '%s\n' "$row" | awk -F '\t' '{print $10}')
  all_present=1
  IFS=',' read -r -a files <<< "$raw"
  for f in "${files[@]}"; do
    [[ -s "data/raw/$f" ]] || all_present=0
  done
  if [[ "$all_present" -eq 1 ]]; then
    echo "SKIP $src: raw material already present ($raw)"
    continue
  fi
  echo "DOWNLOAD $src"
  bash scripts/11_materialize_public_gwas.sh --download "$src"
  echo "MATERIALIZE $src"
  bash scripts/11_materialize_public_gwas.sh --materialize "$src"
done

VERIFY_OUT="frailty_paper/manifests/core_public_gwas_verify.tsv"
if [[ -e "$VERIFY_OUT" ]]; then
  VERIFY_OUT="frailty_paper/manifests/core_public_gwas_verify.refreshed.tsv"
  [[ ! -e "$VERIFY_OUT" ]] || { echo "Refusing to overwrite existing audit output: $VERIFY_OUT" >&2; exit 1; }
fi
bash scripts/11_materialize_public_gwas.sh --verify > "$VERIFY_OUT"
echo "Core registered GWAS audit written to $VERIFY_OUT"
