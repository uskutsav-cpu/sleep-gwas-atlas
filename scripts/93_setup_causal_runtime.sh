#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

INSTALL=true
if [ "${1:-}" = "--no-install" ]; then
  INSTALL=false
  shift
fi
if [ "$#" -ne 0 ]; then
  echo "Usage: bash scripts/93_setup_causal_runtime.sh [--no-install]" >&2
  exit 2
fi

RUNTIME_MANIFEST=config/interpretation_causal_runtime.json
RUNTIME_MANIFEST_SHA256=9a2d0106b5cd2f641b41f7d86447b8305b18abc9c01d8fc4249621860e166d94
DEPENDENCY_MANIFEST=config/interpretation_causal_packages.tsv
DEPENDENCY_MANIFEST_SHA256=6c9ae71bf61cf36de4a952c469ffe4e5dec1b56ff6ef71f80225ff650fe0e26b
INSTALLED_MANIFEST=config/interpretation_causal_installed_packages.tsv
INSTALLED_MANIFEST_SHA256=77e51c68395bf9cd0ceeb1bd04a41256b9b29734a6cda83cd5afeb82164a6e09
LIBRARY="$ROOT/.mr-env/library"
SOURCE_DIR="$ROOT/.mr-env/sources/dependencies"

hash_file() {
  shasum -a 256 "$1" | awk '{print $1}'
}

require_pin() {
  path=$1
  expected_bytes=$2
  expected_sha=$3
  label=$4
  if [ ! -f "$path" ]; then
    echo "ERROR: missing $label: $path" >&2
    exit 1
  fi
  observed_bytes=$(stat -f '%z' "$path")
  observed_sha=$(hash_file "$path")
  if [ "$observed_bytes" != "$expected_bytes" ] || [ "$observed_sha" != "$expected_sha" ]; then
    echo "ERROR: $label differs from its exact byte/SHA-256 pin" >&2
    exit 1
  fi
}

fetch_pin() {
  path=$1
  expected_bytes=$2
  expected_sha=$3
  url=$4
  label=$5
  mkdir -p "$(dirname "$path")"
  if [ ! -f "$path" ]; then
    curl -fL --retry 3 "$url" -o "$path"
  fi
  require_pin "$path" "$expected_bytes" "$expected_sha" "$label"
}

if [ "$(hash_file "$RUNTIME_MANIFEST")" != "$RUNTIME_MANIFEST_SHA256" ]; then
  echo "ERROR: causal runtime manifest differs from the setup-script pin" >&2
  exit 1
fi
if [ "$(hash_file "$DEPENDENCY_MANIFEST")" != "$DEPENDENCY_MANIFEST_SHA256" ]; then
  echo "ERROR: causal dependency source manifest differs from its pin" >&2
  exit 1
fi
if [ "$(hash_file "$INSTALLED_MANIFEST")" != "$INSTALLED_MANIFEST_SHA256" ]; then
  echo "ERROR: causal installed-package manifest differs from its pin" >&2
  exit 1
fi
if [ ! -x .r-env/bin/Rscript ]; then
  echo "ERROR: .r-env/bin/Rscript is absent; create environment/genomicsem.yml first" >&2
  exit 1
fi
R_VERSION=$(.r-env/bin/Rscript --vanilla -e 'cat(paste(R.version$major, R.version$minor, sep="."))')
R_ARCH=$(.r-env/bin/Rscript --vanilla -e 'cat(R.version$arch)')
if [ "$R_VERSION" != 4.3.3 ] || [ "$R_ARCH" != aarch64 ]; then
  echo "ERROR: causal runtime requires R 4.3.3 aarch64; found $R_VERSION $R_ARCH" >&2
  exit 1
fi

mkdir -p "$LIBRARY" "$SOURCE_DIR"
while IFS=$'\t' read -r package version archive bytes checksum url; do
  if [ "$package" = package ]; then
    continue
  fi
  fetch_pin "$ROOT/$archive" "$bytes" "$checksum" "$url" "$package $version source"
done < "$DEPENDENCY_MANIFEST"

fetch_pin \
  "$ROOT/ref/interpretation/causal/twosamplemr_0.7.9/TwoSampleMR_0.7.9.tar.gz" \
  8813664 d4ceabf89d8a21badc57b86a10fa73a6bb2d1931ba43be4f1badeed79cc6a3af \
  'https://mrcieu.r-universe.dev/src/contrib/TwoSampleMR_0.7.9.tar.gz?sha256=d4ceabf89d8a21badc57b86a10fa73a6bb2d1931ba43be4f1badeed79cc6a3af&file=/TwoSampleMR_0.7.9.tar.gz' \
  'TwoSampleMR 0.7.9'
fetch_pin \
  "$ROOT/ref/interpretation/causal/mrpresso_3e3c92d/MR-PRESSO-3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87.tar.gz" \
  453853 8bee27809fdf2ab69e549d08db0ff914d1af80f16119973c56f98c07f54f8728 \
  'https://codeload.github.com/rondolab/MR-PRESSO/tar.gz/3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87' \
  'MR-PRESSO 3e3c92d7'
fetch_pin \
  "$ROOT/ref/interpretation/causal/cause_1.2.0/cause-8d86fc8e9a4509feec84b60b03a3048fe7562c4c.tar.gz" \
  29278307 931c1f8bbdcb27557deb5e56cee1432334ad718269b430fd76e0183a646fed50 \
  'https://codeload.github.com/jean997/cause/tar.gz/8d86fc8e9a4509feec84b60b03a3048fe7562c4c' \
  'cause 1.2.0'
fetch_pin \
  "$ROOT/ref/interpretation/causal/lhcmr_aee7f348/lhcMR-aee7f348744dc7be8b4e7cd52d645ee225cd8b95.tar.gz" \
  123245 0ed84d01fca5a5119ea60a11abdb2224fbf5de5abf6a2631f9183af19a8dcec9 \
  'https://codeload.github.com/LizaDarrous/lhcMR/tar.gz/aee7f348744dc7be8b4e7cd52d645ee225cd8b95' \
  'lhcMR aee7f348'
fetch_pin \
  "$ROOT/ref/interpretation/causal/plink_1.9_b7.11/plink_mac_20250819.zip" \
  2130883 d6c129f940d7235ef24ba2dc5ec926f66690c2a439186a010343844893fda29b \
  'https://s3.amazonaws.com/plink1-assets/plink_mac_20250819.zip' \
  'PLINK 1.9 stable macOS archive'

PLINK="$ROOT/ref/interpretation/causal/plink_1.9_b7.11/plink"
if [ ! -f "$PLINK" ]; then
  unzip -p "$ROOT/ref/interpretation/causal/plink_1.9_b7.11/plink_mac_20250819.zip" plink > "$PLINK"
  chmod 755 "$PLINK"
fi
require_pin "$PLINK" 4344760 5e49297d82c680f2ccfe45c62c974773f6c6e360e5d8a453f831ad95a480b642 'PLINK 1.9 stable universal binary'
if [ "$($PLINK --version)" != 'PLINK v1.9.0-b.7.11 64-bit (19 Aug 2025)' ]; then
  echo "ERROR: PLINK executable reports an unexpected version" >&2
  exit 1
fi

if [ "$INSTALL" = false ]; then
  echo "CAUSAL_RUNTIME_SOURCES_READY install=FALSE"
  exit 0
fi

echo "CAUSAL_RUNTIME_INSTALL_START: compiling the exact dependency closure can take more than 30 minutes."
MAKEVARS=
if [ "$(uname -s)" = Darwin ] && [ "$(uname -m)" = arm64 ]; then
  gcc_root="$ROOT/.r-env/lib/gcc/arm64-apple-darwin20.0.0"
  gcc_dirs=("$gcc_root"/*)
  if [ "${#gcc_dirs[@]}" -ne 1 ] || [ ! -d "${gcc_dirs[0]}" ]; then
    echo "ERROR: expected exactly one versioned conda GCC runtime under $gcc_root" >&2
    exit 1
  fi
  MAKEVARS="$ROOT/.mr-env/Makevars"
  printf 'FLIBS = -L%s -L%s/.r-env/lib -lheapt_w -lgfortran -lquadmath -lm\n' \
    "${gcc_dirs[0]}" "$ROOT" > "$MAKEVARS"
fi

export R_LIBS_USER="$LIBRARY"
install_archive() {
  archive=$1
  if [ -n "$MAKEVARS" ]; then
    R_MAKEVARS_USER="$MAKEVARS" .r-env/bin/R CMD INSTALL --library="$LIBRARY" "$archive"
  else
    .r-env/bin/R CMD INSTALL --library="$LIBRARY" "$archive"
  fi
}

# The source list is alphabetized for auditability. Repeated passes install a
# dependency only after its own dependencies become available, without a live
# package repository or an unpinned solver.
remaining=1
pass=0
while [ "$remaining" -gt 0 ] && [ "$pass" -lt 12 ]; do
  pass=$((pass + 1))
  remaining=0
  progress=0
  while IFS=$'\t' read -r package version archive bytes checksum url; do
    if [ "$package" = package ]; then
      continue
    fi
    installed=$(.r-env/bin/Rscript --vanilla -e ".libPaths(c('$LIBRARY',.libPaths())); cat(if(requireNamespace('$package',quietly=TRUE)) as.character(packageVersion('$package')) else '')")
    if [ "$installed" = "$version" ]; then
      continue
    fi
    remaining=$((remaining + 1))
    if install_archive "$ROOT/$archive" >/dev/null 2>&1; then
      progress=$((progress + 1))
    fi
  done < "$DEPENDENCY_MANIFEST"
  if [ "$remaining" -gt 0 ] && [ "$progress" -eq 0 ]; then
    echo "ERROR: exact dependency installation reached an unresolved fixed point" >&2
    exit 1
  fi
done

install_archive "$ROOT/ref/interpretation/causal/mrpresso_3e3c92d/MR-PRESSO-3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87.tar.gz"
install_archive "$ROOT/ref/interpretation/causal/twosamplemr_0.7.9/TwoSampleMR_0.7.9.tar.gz"
install_archive "$ROOT/ref/interpretation/causal/cause_1.2.0/cause-8d86fc8e9a4509feec84b60b03a3048fe7562c4c.tar.gz"
install_archive "$ROOT/ref/interpretation/causal/lhcmr_aee7f348/lhcMR-aee7f348744dc7be8b4e7cd52d645ee225cd8b95.tar.gz"

.r-env/bin/Rscript --vanilla - "$ROOT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1L]])
.libPaths(c(file.path(root, ".mr-env/library"), .libPaths()))
expected <- read.delim(file.path(root, "config/interpretation_causal_installed_packages.tsv"), stringsAsFactors = FALSE)
observed <- installed.packages()
for (index in seq_len(nrow(expected))) {
  package <- expected$package[[index]]
  if (!package %in% rownames(observed) || observed[package, "Version"] != expected$version[[index]]) {
    stop(sprintf("installed package mismatch: %s expected %s", package, expected$version[[index]]))
  }
  expected_library <- if (expected$library[[index]] == "causal") ".mr-env/library" else ".r-env/lib/R/library"
  if (!endsWith(normalizePath(observed[package, "LibPath"]), expected_library)) {
    stop(sprintf("installed package library mismatch: %s", package))
  }
}
stopifnot(
  as.character(packageVersion("TwoSampleMR")) == "0.7.9",
  as.character(packageVersion("MRPRESSO")) == "1.0",
  as.character(packageVersion("cause")) == "1.2.0",
  as.character(packageVersion("lhcMR")) == "0.0.0.9000",
  is.function(TwoSampleMR::harmonise_data),
  is.function(MRPRESSO::mr_presso),
  is.function(cause::cause),
  is.function(lhcMR::calculate_SP)
)
RSCRIPT

echo "CAUSAL_RUNTIME_READY R=4.3.3 packages=153 plink=1.9.0-b.7.11"
