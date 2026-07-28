#!/usr/bin/env bash
# One-time setup. Run from the repo root:  bash scripts/00_setup.sh
# Requires conda (miniforge/miniconda) and ~8 GB free disk.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)

echo "==> [1/4] LDSC (needs Python 2.7 — that's why it gets its own conda env)"
if [ ! -d ldsc ]; then
  git clone https://github.com/bulik/ldsc.git
fi
if ! conda env list | grep -q '^ldsc '; then
  conda env create --file ldsc/environment.yml   # creates env named 'ldsc'
else
  echo "    conda env 'ldsc' already exists, skipping"
fi
# Sanity check: this must print the LDSC help banner.
conda run -n ldsc python ldsc/ldsc.py -h > /dev/null && echo "    LDSC OK"

echo "==> [2/4] Reference files -> ref/"
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

echo "==> [3/4] Python deps for the harmonize/collate/plot scripts"
python3 -m pip install --quiet --break-system-packages pandas numpy scipy matplotlib

echo "==> [4/4] Done."
cat <<'EOF'

NOTE ON DOWNLOAD MIRRORS
  The Broad hosting for LD scores has moved more than once. If both URLs
  above 404, get eur_w_ld_chr and w_hm3.snplist from the LDSC wiki:
    https://github.com/bulik/ldsc/wiki
  Do NOT substitute a different LD reference without telling your mentor —
  the LD panel must match the ancestry of your sumstats (EUR here).

NEXT
  1. Put raw sumstats in data/raw/, named as in config/traits.tsv
  2. bash scripts/02_munge.sh insomnia mdd
  3. bash scripts/03_h2_qc.sh insomnia mdd
  4. bash scripts/04_rg.sh insomnia mdd
EOF
