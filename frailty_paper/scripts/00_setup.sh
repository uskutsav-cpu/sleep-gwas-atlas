#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
cd "$ROOT"
if [[ ! -f config/traits.tsv || ! -f config/public_gwas_sources.tsv ]]; then
  echo "ERROR: run this from the sleep-gwas-atlas repository root." >&2
  exit 1
fi
PAPER="$ROOT/frailty_paper"
mkdir -p "$PAPER"/{review/{pubmed,manual_exports,screening,fulltext},data/{gwas,qtl,single_cell},manifests,provenance,logs,scripts,config}
VENV_DIR="${FRAILTY_VENV_DIR:-$PAPER/.venv}"
case "$VENV_DIR" in
  /*) ;;
  *) VENV_DIR="$ROOT/$VENV_DIR" ;;
esac
PYTHON_BIN="${FRAILTY_PYTHON:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: workflow Python $PYTHON_BIN was not found; set FRAILTY_PYTHON to the pinned interpreter." >&2
  exit 1
fi
PYTHON_VERSION="$($PYTHON_BIN --version 2>&1 | awk '{print $2}')"
if [[ "$PYTHON_VERSION" != "3.11.11" ]]; then
  echo "ERROR: workflow Python must be 3.11.11 (environment/tool_versions.tsv); found $PYTHON_VERSION." >&2
  echo "Install the pinned interpreter and set FRAILTY_PYTHON=/path/to/python3.11.11." >&2
  exit 1
fi
if [[ -e "$VENV_DIR" ]]; then
  if [[ ! -x "$VENV_DIR/bin/python" ]]; then
    echo "ERROR: environment target exists but is not a Python venv: $VENV_DIR" >&2
    echo "Choose an empty FRAILTY_VENV_DIR; existing contents were left unchanged." >&2
    exit 1
  fi
  EXISTING_VERSION="$("$VENV_DIR/bin/python" --version 2>&1 | awk '{print $2}')"
  if [[ "$EXISTING_VERSION" != "3.11.11" ]]; then
    echo "ERROR: refusing to layer Python 3.11.11 over existing $EXISTING_VERSION environment: $VENV_DIR" >&2
    echo "Preserve that environment and choose an empty FRAILTY_VENV_DIR for a clean rebuild." >&2
    exit 1
  fi
else
  mkdir -p "$(dirname "$VENV_DIR")"
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip wheel
python -m pip install -r "$ROOT/requirements-pipeline.txt"
python -m pip install -r "$ROOT/requirements-workflow.txt"
python -m pip install -r "$ROOT/requirements-acquisition.txt"
printf '\nEnvironment ready at %s. Activate later with:\n  source "%s/bin/activate"\n' "$VENV_DIR" "$VENV_DIR"
