#!/usr/bin/env bash
# Exercise the complete Phase 0/1 code path with deliberately fake data.
# All artifacts remain in results/_smoketest/ and carry a synthetic watermark.
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON_BIN=${PYTHON_BIN:-python3}
SMOKE_DIR=${SMOKE_DIR:-results/_smoketest}
PANEL=config/analysis_panel.tsv
export MPLCONFIGDIR=${MPLCONFIGDIR:-"${TMPDIR:-/tmp}/sleep-gwas-atlas-matplotlib"}
mkdir -p "$MPLCONFIGDIR"
case "$SMOKE_DIR" in
  results/_smoketest|results/_smoketest/*) ;;
  *) echo "ERROR: SMOKE_DIR must remain under results/_smoketest/" >&2; exit 1 ;;
esac

"$PYTHON_BIN" scripts/00_validate_panel.py --manifest "$PANEL" \
  --write-provenance "$SMOKE_DIR/analysis_panel_provenance.tsv"

"$PYTHON_BIN" scripts/make_test_data.py --logs --raw --out data/_test_raw --logdir "$SMOKE_DIR"
"$PYTHON_BIN" scripts/01_harmonize.py --trait insomnia --infile data/_test_raw/insomnia.txt.gz --outdir data/_test_harmonized
"$PYTHON_BIN" scripts/01_harmonize.py --trait bipolar --infile data/_test_raw/bipolar.txt.gz --outdir data/_test_harmonized
"$PYTHON_BIN" scripts/05_collate.py --mode h2 --logdir "$SMOKE_DIR" --out "$SMOKE_DIR/h2_summary.tsv"
"$PYTHON_BIN" scripts/make_smoke_readiness.py --h2 "$SMOKE_DIR/h2_summary.tsv" \
  --out "$SMOKE_DIR/trait_readiness.tsv"
"$PYTHON_BIN" scripts/09_select_phase1_traits.py --readiness "$SMOKE_DIR/trait_readiness.tsv" \
  --h2 "$SMOKE_DIR/h2_summary.tsv" --out "$SMOKE_DIR/phase1_inclusion.tsv"
"$PYTHON_BIN" scripts/05_collate.py --mode rg --logdir "$SMOKE_DIR" --out "$SMOKE_DIR/rg_matrix.tsv"
"$PYTHON_BIN" scripts/06_heatmap.py --rg "$SMOKE_DIR/rg_matrix.tsv" --config "$PANEL" \
  --out "$SMOKE_DIR/fig2_smoketest.png" --provenance-label "SYNTHETIC - NOT REAL DATA"
"$PYTHON_BIN" scripts/07_export_metadata.py --h2 "$SMOKE_DIR/h2_summary.tsv" --out "$SMOKE_DIR/gwas_metadata_table.tsv"
"$PYTHON_BIN" scripts/08_phase1_report.py --h2 "$SMOKE_DIR/h2_summary.tsv" \
  --rg "$SMOKE_DIR/rg_matrix.tsv" \
  --out "$SMOKE_DIR/phase1_summary_report.md"
"$PYTHON_BIN" scripts/10_phase0_audit.py --out "$SMOKE_DIR/source_readiness.tsv"
touch "$SMOKE_DIR/SMOKE_TEST_OK"
echo "Smoke test completed. All outputs are synthetic and confined to $SMOKE_DIR/."
