#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 baseline_analysis.py
python3 prepare_data.py
python3 independent_item_audit.py
Rscript survey_analysis.R
python3 assemble_supplement.py
python3 genetic_reference_check.py
python3 figures.py
python3 verify_final.py
