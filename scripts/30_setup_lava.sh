#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

COMMIT=e729a245f7b6923967a96804fbf5246eadf2d6c6
ARCHIVE_SHA256=15477c9547c3533d681cbfb6491508a29342fe99b725845054887b9d29dc2bec
LOCUS_SHA256=462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882
ARCHIVE_URL="https://github.com/josefin-werme/LAVA/archive/${COMMIT}.tar.gz"
KEEP_URL=https://cran.r-project.org/src/contrib/keep_1.0.tar.gz
KEEP_SHA256=4eaed49d4a5fe470c906435de119e9f4681a5c98302a6c762669279d211c0089
MATRIXSAMPLING_URL=https://cran.r-project.org/src/contrib/matrixsampling_2.0.0.tar.gz
MATRIXSAMPLING_SHA256=eade38bff0b162bd5dda6a8ab1cdebdac310dc54639a78240c5d1e20f3b1ab15
CPP11_URL=https://cran.r-project.org/src/contrib/Archive/cpp11/cpp11_0.5.2.tar.gz
CPP11_SHA256=0e8ac07f9d599b82e7a811f9d084e5125ae787b1ba04e5ba57f79e2642af091b
ARCHIVE=
INSTALL=true

usage() {
  echo "Usage: bash scripts/30_setup_lava.sh [--archive FILE] [--no-install]"
  echo "Installs exact LAVA ${COMMIT} and extracts its checksum-locked 2,495-locus file."
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --archive) ARCHIVE=$2; shift 2 ;;
    --no-install) INSTALL=false; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ ! -x .r-env/bin/Rscript ]; then
  echo "ERROR: .r-env/bin/Rscript is missing; create environment/genomicsem.yml first" >&2
  exit 1
fi

WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/sleep-atlas-lava-setup.XXXXXX")
cleanup() {
  rm -rf "$WORK_DIR"
}
trap cleanup EXIT

if [ -z "$ARCHIVE" ]; then
  ARCHIVE="$WORK_DIR/LAVA-${COMMIT}.tar.gz"
  curl -fL --retry 3 "$ARCHIVE_URL" -o "$ARCHIVE"
fi
if [ ! -f "$ARCHIVE" ]; then
  echo "ERROR: archive not found: $ARCHIVE" >&2
  exit 1
fi

ACTUAL_ARCHIVE_SHA256=$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')
if [ "$ACTUAL_ARCHIVE_SHA256" != "$ARCHIVE_SHA256" ]; then
  echo "ERROR: LAVA archive SHA-256 mismatch" >&2
  exit 1
fi

SOURCE_DIR="$WORK_DIR/LAVA-${COMMIT}"
tar -xzf "$ARCHIVE" -C "$WORK_DIR"
LOCUS_MEMBER="LAVA-${COMMIT}/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
mkdir -p ref/lava
LOCUS_TMP="$WORK_DIR/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
tar -xOf "$ARCHIVE" "$LOCUS_MEMBER" > "$LOCUS_TMP"
ACTUAL_LOCUS_SHA256=$(shasum -a 256 "$LOCUS_TMP" | awk '{print $1}')
if [ "$ACTUAL_LOCUS_SHA256" != "$LOCUS_SHA256" ]; then
  echo "ERROR: LAVA locus-definition SHA-256 mismatch" >&2
  exit 1
fi
if [ "$(wc -l < "$LOCUS_TMP" | tr -d ' ')" != 2496 ]; then
  echo "ERROR: expected header plus 2,495 loci" >&2
  exit 1
fi
mv "$LOCUS_TMP" ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile

if [ "$INSTALL" = true ]; then
  R_VERSION=$(.r-env/bin/Rscript -e 'cat(R.version$major, substr(R.version$minor, 1, 1), sep=".")')
  if [ "$R_VERSION" != 4.3 ]; then
    echo "ERROR: LAVA runtime must use R 4.3.x; found $R_VERSION" >&2
    exit 1
  fi
  install_cran_exact() {
    package=$1
    version=$2
    url=$3
    expected_sha=$4
    installed=$(.r-env/bin/Rscript -e "cat(if (requireNamespace('$package', quietly=TRUE)) as.character(packageVersion('$package')) else '')")
    if [ "$installed" = "$version" ]; then
      return
    fi
    dependency_archive="$WORK_DIR/${package}_${version}.tar.gz"
    curl -fL --retry 3 "$url" -o "$dependency_archive"
    dependency_sha=$(shasum -a 256 "$dependency_archive" | awk '{print $1}')
    if [ "$dependency_sha" != "$expected_sha" ]; then
      echo "ERROR: $package $version source SHA-256 mismatch" >&2
      exit 1
    fi
    .r-env/bin/R CMD INSTALL --library="$ROOT/.r-env/lib/R/library" "$dependency_archive"
  }
  install_cran_exact keep 1.0 "$KEEP_URL" "$KEEP_SHA256"
  install_cran_exact matrixsampling 2.0.0 "$MATRIXSAMPLING_URL" "$MATRIXSAMPLING_SHA256"
  install_cran_exact cpp11 0.5.2 "$CPP11_URL" "$CPP11_SHA256"
  if [ "$(uname -s)" = Darwin ] && [ "$(uname -m)" = arm64 ]; then
    R_MAKEVARS_USER="$ROOT/environment/macos-arm64.Makevars" \
      .r-env/bin/R CMD INSTALL --library="$ROOT/.r-env/lib/R/library" "$SOURCE_DIR"
  else
    .r-env/bin/R CMD INSTALL --library="$ROOT/.r-env/lib/R/library" "$SOURCE_DIR"
  fi
fi

.r-env/bin/Rscript -e '
  library(LAVA)
  stopifnot(as.character(packageVersion("LAVA")) == "0.1.5")
  stopifnot(as.character(packageVersion("matrixsampling")) == "2.0.0")
  stopifnot(as.character(packageVersion("cpp11")) == "0.5.2")
  stopifnot(as.character(packageVersion("keep")) == "1.0")
'
.r-env/bin/Rscript scripts/133_validate_track_b_lava_runtime.R
echo "LAVA 0.1.5 and the 2,495-locus GRCh37 file are ready."
echo "The 15 GiB UK Biobank LD reference is deliberately not downloaded by this setup."
