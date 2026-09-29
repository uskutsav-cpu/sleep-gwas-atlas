#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
cd "$ROOT"
PAPER="$ROOT/frailty_paper"
# shellcheck disable=SC1091
source "$PAPER/.venv/bin/activate"
bash "$PAPER/scripts/09_require_storage.sh" "$ROOT"

# Fried Frailty Score (FFS), Ye et al. / GCST90295968, 386,565 UKB participants.
# The public Figshare route is currently challenge-blocked; do not retry it here.
FFS_DIR="$PAPER/data/gwas/fried_frailty_score"
mkdir -p "$FFS_DIR"
FFS_FILE="$FFS_DIR/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv"
if [[ ! -s "$FFS_FILE" ]]; then
  echo "BLOCKED: Fried Frailty Score remains unavailable through its challenge-blocked Figshare route; skipping without treating the zero-byte placeholder as data." >&2
else
  echo "NOTICE: Fried Frailty Score file exists locally but requires checksum/provenance verification before use: $FFS_FILE" >&2
fi

# Five physical-frailty components: weight loss, exhaustion, low activity,
# slow walking speed, low grip strength.
python "$PAPER/scripts/03_download_zenodo_record.py" \
  --record 14011550 \
  --outdir "$PAPER/data/gwas/physical_frailty_components"

# Healthspan GWAS resource used in aging-genetics studies.
python "$PAPER/scripts/03_download_zenodo_record.py" \
  --record 1302861 \
  --outdir "$PAPER/data/gwas/healthspan"

echo "Downloadable physical-component and healthspan resources collected; Fried Frailty Score remains blocked pending an accessible, verifiable source."
