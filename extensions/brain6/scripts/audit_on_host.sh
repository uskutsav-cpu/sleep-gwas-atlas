#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo 'Usage: bash audit_on_host.sh /absolute/path/to/sleep-gwas-atlas [new-output-directory]' >&2
  exit 2
fi
REPO="$(cd "$1" && pwd)"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${2:-$HOME/Downloads/brain6-host-audit-$(date +%Y%m%d-%H%M%S)}"
PY="${PYTHON:-$HERE/.venv/bin/python}"
if [ ! -x "$PY" ]; then PY="$(command -v python3)"; fi
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
"$PY" - "$REPO" <<'PYTHON'
import sys, subprocess
from pathlib import Path
from brain6.install import repository_identity
root=Path(sys.argv[1]).resolve()
g=subprocess.run(['git','-C',str(root),'rev-parse','--show-toplevel'],check=True,capture_output=True,text=True)
if Path(g.stdout.strip()).resolve()!=root:raise SystemExit('Pass the repository root, not a subdirectory')
g=subprocess.run(['git','-C',str(root),'remote','get-url','origin'],check=True,capture_output=True,text=True)
if repository_identity(g.stdout.strip())!='uskutsav-cpu/sleep-gwas-atlas':raise SystemExit('Wrong repository: no audit launched')
PYTHON
"$PY" -m brain6 host-audit24 --repo "$REPO" --out "$OUT"
printf '\nAudit saved at: %s\nNo production analyses, commits, or pushes were performed.\n' "$OUT"
