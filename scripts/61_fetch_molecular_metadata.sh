#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
POLICY="$ROOT/config/molecular_analysis_policy.json"
DEST="$ROOT/ref/molecular/eqtl_catalogue_r7"

if [ ! -f "$POLICY" ]; then
  echo "ERROR: molecular policy is missing: $POLICY" >&2
  exit 1
fi

mkdir -p "$DEST"
python3 - "$POLICY" "$DEST" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

policy_path = Path(sys.argv[1])
destination = Path(sys.argv[2])
policy = json.loads(policy_path.read_text(encoding="utf-8"))
for asset in policy["metadata_assets"]:
    final = policy_path.parent.parent / asset["path"]
    final.parent.mkdir(parents=True, exist_ok=True)
    if final.is_file():
        digest = hashlib.sha256(final.read_bytes()).hexdigest()
        if final.stat().st_size == asset["bytes"] and digest == asset["sha256"]:
            print(f"verified {asset['id']} {final}")
            continue
        raise SystemExit(f"ERROR: existing molecular metadata differs from lock: {final}")
    descriptor, temporary_name = tempfile.mkstemp(prefix=final.name + ".", dir=final.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        subprocess.run(["curl", "-fL", "--retry", "3", "--output", str(temporary), asset["url"]], check=True)
        digest = hashlib.sha256(temporary.read_bytes()).hexdigest()
        if temporary.stat().st_size != asset["bytes"] or digest != asset["sha256"]:
            raise SystemExit(f"ERROR: downloaded molecular metadata differs from lock: {asset['id']}")
        os.replace(temporary, final)
        print(f"installed {asset['id']} {final}")
    finally:
        temporary.unlink(missing_ok=True)
PY

echo "MOLECULAR_METADATA_OK"
