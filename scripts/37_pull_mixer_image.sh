#!/usr/bin/env bash
set -euo pipefail

IMAGE=ghcr.io/precimed/gsa-mixer@sha256:5bf54ddd6f7f81b93eeb5450b1a6dc809d1fcf926b3a815d6d514f36ce04f51c
PULL=false
REQUIRE_PRESENT=false

for argument in "$@"; do
  case "$argument" in
    --pull) PULL=true ;;
    --require-present) REQUIRE_PRESENT=true ;;
    *) echo "Usage: bash scripts/37_pull_mixer_image.sh [--pull] [--require-present]" >&2; exit 2 ;;
  esac
done

echo "Pinned MiXeR image: $IMAGE"
echo "Compressed image layers: 2,106,984,629 bytes (1.96 GiB)"
echo "Supported image platform: linux/amd64 only"
if [ "$(uname -m)" != x86_64 ]; then
  echo "ERROR: official MiXeR containers do not support this $(uname -m) host" >&2
  exit 1
fi
if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: Docker is not available" >&2
  exit 1
fi
if docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "Pinned image is already present."
  exit 0
fi
if [ "$PULL" != true ]; then
  if [ "$REQUIRE_PRESENT" = true ]; then
    echo "ERROR: pinned MiXeR image is absent and no pull was authorized" >&2
    exit 1
  fi
  echo "No pull requested. Re-run with --pull on the supported analysis host."
  exit 0
fi
docker pull "$IMAGE"
docker image inspect "$IMAGE" >/dev/null
echo "Pinned MiXeR image is ready."
