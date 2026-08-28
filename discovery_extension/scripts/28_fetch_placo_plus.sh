#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
DESTINATION="${ROOT_DIR}/.r-env/share/placo/PLACO_v0.2.0.R"
SOURCE_URL="https://raw.githubusercontent.com/RayDebashree/PLACO/3ba3cae1d323ad117fb4540e620bcefa79f70663/PLACO_v0.2.0.R"
EXPECTED_SHA256="fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124"

mkdir -p "$(dirname "${DESTINATION}")"
TEMPORARY="${DESTINATION}.partial"
trap 'rm -f "${TEMPORARY}"' EXIT
curl -fL --retry 4 --retry-delay 2 "${SOURCE_URL}" -o "${TEMPORARY}"
OBSERVED_SHA256="$(shasum -a 256 "${TEMPORARY}" | awk '{print $1}')"
if [[ "${OBSERVED_SHA256}" != "${EXPECTED_SHA256}" ]]; then
  echo "ERROR: PLACO+ source SHA256 mismatch" >&2
  exit 1
fi
install -m 0444 "${TEMPORARY}" "${DESTINATION}"
rm -f "${TEMPORARY}"
trap - EXIT
echo "PLACO_PLUS_SOURCE_INSTALLED version=0.2.0 sha256=${OBSERVED_SHA256} path=${DESTINATION}"
