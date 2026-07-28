#!/usr/bin/env bash
# One-time setup. Run from the repo root: bash scripts/00_setup.sh
#
# This installs a maintained Python 3 LDSC implementation plus the EUR
# reference files. It needs a network connection and several GB of free disk;
# inspect the source URLs below before starting a real analysis.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}
VENV_DIR=${VENV_DIR:-.venv}
LDSC_DIR=${LDSC_DIR:-ldsc}

echo "==> [1/4] Python environment"
if [ ! -x "$VENV_DIR/bin/python" ]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi
"$VENV_DIR/bin/python" -m pip install --upgrade pip
"$VENV_DIR/bin/python" -m pip install -r requirements-pipeline.txt

echo "==> [2/4] LDSC (maintained Python 3 implementation)"
if [ ! -d "$LDSC_DIR/.git" ]; then
  git clone --branch ldsc39 --depth 1 https://github.com/CBIIT/ldsc.git "$LDSC_DIR"
else
  echo "    $LDSC_DIR already exists, leaving its checked-out revision unchanged"
fi
"$VENV_DIR/bin/python" -m pip install -r "$LDSC_DIR/requirements.txt"
"$VENV_DIR/bin/python" "$LDSC_DIR/ldsc.py" -h > /dev/null
"$VENV_DIR/bin/python" "$LDSC_DIR/munge_sumstats.py" -h > /dev/null
echo "    LDSC OK"

echo "==> [3/4] Reference files -> ref/"
mkdir -p ref && cd ref
# European LD scores computed on 1000G Phase 3, HapMap3 SNPs.
if [ ! -d eur_w_ld_chr ]; then
  wget -q https://zenodo.org/records/10515792/files/eur_w_ld_chr.tar.bz2 \
    || wget -q https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2
  tar -jxf eur_w_ld_chr.tar.bz2
fi
# HapMap3 SNP list used by munge_sumstats.py as the --merge-alleles target.
if [ ! -f w_hm3.snplist ]; then
  wget -q https://zenodo.org/records/7796478/files/w_hm3.snplist.bz2 \
    || wget -q https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/w_hm3.snplist.bz2
  bunzip2 -f w_hm3.snplist.bz2
fi
cd "$ROOT"
echo "    ref/ contains:"; ls ref | sed 's/^/      /'

echo "==> [4/4] Done."
cat <<'EOF'

NOTE ON DOWNLOAD MIRRORS
  The Broad hosting for LD scores has moved more than once. If both URLs
  above 404, get eur_w_ld_chr and w_hm3.snplist from the LDSC wiki:
    https://github.com/bulik/ldsc/wiki
  Do NOT substitute a different LD reference without telling your mentor —
  the LD panel must match the ancestry of your sumstats (EUR here).

NEXT
  1. Put only verified hg19, EUR raw sumstats in data/raw/, named as in
     config/traits.tsv. Do not use an hg38 file without an explicit,
     documented liftover decision.
  2. export PYTHON_BIN=.venv/bin/python LDSC_PYTHON=.venv/bin/python
     LDSC_DIR=ldsc
  3. bash scripts/02_munge.sh insomnia mdd
  4. bash scripts/03_h2_qc.sh insomnia mdd
  5. bash scripts/04_rg.sh --h2 results/tables/h2_summary.tsv
EOF
