#!/usr/bin/env bash
# One-time setup. Run from the repo root: bash scripts/00_setup.sh
#
# This installs the maintained Python 3 LDSC implementation plus the EUR
# reference files. It needs a network connection and several GB of free disk;
# inspect the source URLs below before starting a real analysis.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
CONDA_BIN=${CONDA_BIN:-conda}
LDSC_ENV_DIR=${LDSC_ENV_DIR:-.ldsc-env}
LDSC_DIR=${LDSC_DIR:-ldsc}
LDSC_REPOSITORY=https://github.com/CBIIT/ldsc.git
LDSC_COMMIT=6c673952cee74bd5c57aef1555a03b1c015399a0
REFERENCE_URL=${REFERENCE_URL:-https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1}
REFERENCE_MD5=${REFERENCE_MD5:-e2f16343c4cfaa76caa7d0c03d26b489}
REFERENCE_SHA256=${REFERENCE_SHA256:-9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069}
VARIANT_MAP_SHA256=${VARIANT_MAP_SHA256:-6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4}

echo "==> [1/5] Python 3.9 LDSC environment"
if ! command -v "$CONDA_BIN" > /dev/null 2>&1; then
  echo "ERROR: conda is required to create the reproducible Python 3.9 LDSC environment." >&2
  echo "Install Miniforge/conda, or set CONDA_BIN to its executable path." >&2
  exit 1
fi
PACKAGES=(
  python=3.9.23 numpy=1.21.5 pandas=1.3.3 scipy=1.7.3
  python-dateutil=2.8.2 pytz=2022.7.1 bitarray=2.8.3 nose=1.3.7
  pybedtools=0.10.0 flask=2.3.3 requests=2.31.0 matplotlib=3.7.5
)
if [ ! -x "$LDSC_ENV_DIR/bin/python" ]; then
  # Direct community channels avoid implicitly accepting Anaconda's commercial
  # channel Terms of Service. Versions follow CBIIT/ldsc's environment3.yml;
  # matplotlib is added for this repository's reporting scripts.
  "$CONDA_BIN" create --yes --override-channels \
    --channel conda-forge --channel bioconda --prefix "$LDSC_ENV_DIR" \
    "${PACKAGES[@]}"
else
  # Also repair a partly created environment (for example after an interrupted
  # initial setup) instead of assuming the directory alone proves completeness.
  "$CONDA_BIN" install --yes --override-channels \
    --channel conda-forge --channel bioconda --prefix "$LDSC_ENV_DIR" \
    "${PACKAGES[@]}"
fi
LDSC_PYTHON="$LDSC_ENV_DIR/bin/python"
"$LDSC_PYTHON" -c 'import numpy, pandas, scipy, matplotlib, pybedtools'

echo "==> [2/5] LDSC (maintained Python 3 implementation)"
if [ ! -d "$LDSC_DIR/.git" ]; then
  git clone "$LDSC_REPOSITORY" "$LDSC_DIR"
  git -C "$LDSC_DIR" checkout --detach "$LDSC_COMMIT"
fi
ACTUAL_LDSC_COMMIT=$(git -C "$LDSC_DIR" rev-parse HEAD)
if [ "$ACTUAL_LDSC_COMMIT" != "$LDSC_COMMIT" ]; then
  echo "ERROR: LDSC must be pinned at $LDSC_COMMIT, found $ACTUAL_LDSC_COMMIT" >&2
  echo "       Move the existing LDSC checkout aside and rerun setup." >&2
  exit 1
fi
echo "    pinned LDSC commit: $ACTUAL_LDSC_COMMIT"
# CBIIT/ldsc's ldsc39 revision opens compressed summary statistics in binary
# mode while parsing their header. Python 3 therefore raises TypeError before
# analysis. Apply this narrowly scoped, version-checked compatibility patch;
# setup deliberately aborts rather than silently applying it to another LDSC
# revision. The patch only changes the header stream to text mode.
LDSC_COMPAT_PATCH="$ROOT/patches/ldsc39-python3-compressed-header.patch"
if git -C "$LDSC_DIR" apply --unidiff-zero --reverse --check "$LDSC_COMPAT_PATCH"; then
  echo "    Python 3 compressed-header compatibility patch already applied"
elif git -C "$LDSC_DIR" apply --unidiff-zero --check "$LDSC_COMPAT_PATCH"; then
  git -C "$LDSC_DIR" apply --unidiff-zero "$LDSC_COMPAT_PATCH"
  echo "    applied Python 3 compressed-header compatibility patch"
else
  echo "ERROR: the checked-out LDSC revision does not match the reviewed compressed-header patch." >&2
  echo "       Refusing to run with an unreviewed LDSC source revision." >&2
  exit 1
fi
"$LDSC_PYTHON" "$LDSC_DIR/ldsc.py" -h > /dev/null
"$LDSC_PYTHON" "$LDSC_DIR/munge_sumstats.py" -h > /dev/null
echo "    LDSC OK"

echo "==> [3/5] Reference files -> ref/"
mkdir -p ref && cd ref
# European LD scores computed on 1000G Phase 3, HapMap3 SNPs. Zenodo record
# 8182036 explicitly documents this archive as a gzip copy of the original
# Alkes-group distribution and supplies the MD5 below.
if [ ! -d eur_w_ld_chr ]; then
  ARCHIVE=eur_w_ld_chr.tar.gz
  curl --fail --location --retry 3 --output "$ARCHIVE" "$REFERENCE_URL"
  ACTUAL_MD5=$("$LDSC_PYTHON" -c 'import hashlib, sys; h = hashlib.md5(); f = open(sys.argv[1], "rb"); [h.update(chunk) for chunk in iter(lambda: f.read(1024 * 1024), b"")]; print(h.hexdigest())' "$ARCHIVE")
  if [ "$ACTUAL_MD5" != "$REFERENCE_MD5" ]; then
    echo "ERROR: EUR LD archive MD5 mismatch: expected $REFERENCE_MD5, got $ACTUAL_MD5" >&2
    exit 1
  fi
  ACTUAL_SHA256=$("$LDSC_PYTHON" -c 'import hashlib, sys; h = hashlib.sha256(); f = open(sys.argv[1], "rb"); [h.update(chunk) for chunk in iter(lambda: f.read(1024 * 1024), b"")]; print(h.hexdigest())' "$ARCHIVE")
  if [ "$ACTUAL_SHA256" != "$REFERENCE_SHA256" ]; then
    echo "ERROR: EUR LD archive SHA-256 mismatch: expected $REFERENCE_SHA256, got $ACTUAL_SHA256" >&2
    exit 1
  fi
  tar -xzf "$ARCHIVE"
  rm -f "$ARCHIVE"
fi
# This exact HapMap3 list ships in the verified EUR archive. Link rather than
# downloading a separately versioned list so munging and LD scores agree.
for chromosome in {1..22}; do
  test -s "eur_w_ld_chr/${chromosome}.l2.ldscore.gz"
done
test -s eur_w_ld_chr/w_hm3.snplist
ln -sfn eur_w_ld_chr/w_hm3.snplist w_hm3.snplist
cd "$ROOT"
echo "    ref/ contains:"; ls ref | sed 's/^/      /'

echo "==> [4/5] Audited HapMap3 GRCh37 identity map"
"$LDSC_PYTHON" scripts/19_build_hm3_variant_map.py \
  --reference-dir ref/eur_w_ld_chr \
  --out ref/hm3_grch37_variant_map.tsv.gz
ACTUAL_MAP_SHA256=$("$LDSC_PYTHON" -c 'import hashlib, sys; h = hashlib.sha256(); f = open(sys.argv[1], "rb"); [h.update(chunk) for chunk in iter(lambda: f.read(1024 * 1024), b"")]; print(h.hexdigest())' ref/hm3_grch37_variant_map.tsv.gz)
if [ "$ACTUAL_MAP_SHA256" != "$VARIANT_MAP_SHA256" ]; then
  echo "ERROR: derived variant-map SHA-256 mismatch: expected $VARIANT_MAP_SHA256, got $ACTUAL_MAP_SHA256" >&2
  exit 1
fi

echo "==> [5/5] Done."
cat <<'EOF'

NOTE ON DOWNLOAD MIRRORS
  The script downloads Zenodo record 8182036 and verifies its published MD5.
  Its provenance is recorded in docs/reference_panel_provenance.md.
  Do NOT substitute a different LD reference without telling your mentor —
  the LD panel must match the ancestry of your sumstats (EUR here).

NEXT
  1. Put only verified hg19, EUR raw sumstats in data/raw/, named as in
     config/analysis_panel.tsv. Do not use an hg38 file without an explicit,
     documented liftover decision.
  2. export PYTHON_BIN=.ldsc-env/bin/python LDSC_PYTHON=.ldsc-env/bin/python
     LDSC_DIR=ldsc
  3. bash scripts/02_munge.sh insomnia bipolar
  4. bash scripts/03_h2_qc.sh insomnia bipolar
  5. bash scripts/04_rg.sh --h2 results/tables/h2_summary.tsv
EOF
