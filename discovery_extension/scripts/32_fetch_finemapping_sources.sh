#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
archive_dir="${repo_root}/.r-env/share/finemapping/source_archives"
r_binary="${repo_root}/.r-env/bin/R"
rscript_binary="${repo_root}/.r-env/bin/Rscript"
library_dir="${repo_root}/.r-env/lib/R/library"
install_packages=false
if [[ "${1:-}" == "--install" ]]; then
  install_packages=true
elif [[ $# -ne 0 ]]; then
  echo "usage: $0 [--install]" >&2
  exit 2
fi

mkdir -p "${archive_dir}"

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    echo "ERROR: neither sha256sum nor shasum is available" >&2
    return 1
  fi
}

fetch_locked() {
  local url="$1"
  local fallback_url="$2"
  local destination="$3"
  local expected="$4"
  if [[ -f "${destination}" ]] && [[ "$(sha256_file "${destination}")" == "${expected}" ]]; then
    echo "SOURCE_ARCHIVE_PRESENT $(basename "${destination}")"
    return
  fi
  local partial="${destination}.part.$$"
  trap 'rm -f "${partial}"' RETURN
  if ! curl -fL --retry 3 --output "${partial}" "${url}"; then
    if [[ "${fallback_url}" == "NA" ]]; then
      return 1
    fi
    curl -fL --retry 3 --output "${partial}" "${fallback_url}"
  fi
  local observed
  observed="$(sha256_file "${partial}")"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "ERROR: source archive checksum mismatch for ${url}: ${observed}" >&2
    return 1
  fi
  mv "${partial}" "${destination}"
  trap - RETURN
  echo "SOURCE_ARCHIVE_FETCHED $(basename "${destination}") sha256=${observed}"
}

susie_archive="${archive_dir}/susieR_0.14.2.tar.gz"
coloc_archive="${archive_dir}/coloc_v5.2.3.tar.gz"
fetch_locked \
  "https://cran.r-project.org/src/contrib/susieR_0.14.2.tar.gz" \
  "https://cran.r-project.org/src/contrib/Archive/susieR/susieR_0.14.2.tar.gz" \
  "${susie_archive}" \
  "ba02322eb1f7a7cc024c9278aa7903a34d8ad5d6f3b12c168374bc6214ed2c6e"
fetch_locked \
  "https://github.com/chr1swallace/coloc/archive/refs/tags/v5.2.3.tar.gz" \
  "NA" \
  "${coloc_archive}" \
  "9cfc21315468bbe984faf9facaaa76e0127999ff4cf798a0a41cec904dd86c60"

if [[ "${install_packages}" == true ]]; then
  if [[ ! -x "${r_binary}" || ! -x "${rscript_binary}" || ! -d "${library_dir}" ]]; then
    echo "ERROR: repository R environment is absent; bootstrap .r-env before installing fine-mapping packages" >&2
    exit 1
  fi
  "${r_binary}" CMD INSTALL --library="${library_dir}" "${susie_archive}"
  "${r_binary}" CMD INSTALL --library="${library_dir}" "${coloc_archive}"
fi

if [[ -x "${rscript_binary}" ]]; then
  "${rscript_binary}" -e \
    'suppressPackageStartupMessages(library(susieR)); suppressPackageStartupMessages(library(coloc)); stopifnot(packageVersion("susieR")=="0.14.2", packageVersion("coloc")=="5.2.3"); cat("FINEMAPPING_RUNTIME_OK susieR=0.14.2 coloc=5.2.3\n")'
else
  echo "SOURCE_ARCHIVES_OK runtime_not_checked=true"
fi
