#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
POLICY="$ROOT/config/molecular_analysis_policy.json"
MODELS=false
RUNTIME=false
ACK_LARGE=false
for argument in "$@"; do
  case "$argument" in
    --models) MODELS=true ;;
    --runtime) RUNTIME=true ;;
    --acknowledge-large-download) ACK_LARGE=true ;;
    *) echo "ERROR: unknown argument: $argument" >&2; exit 1 ;;
  esac
done

python3 - "$ROOT" "$POLICY" "$MODELS" "$ACK_LARGE" <<'PY'
import hashlib
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

root = Path(sys.argv[1])
policy = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
fetch_models = sys.argv[3] == "true"
acknowledge_large = sys.argv[4] == "true"
twas = policy["twas"]

def digest(path, algorithm="sha256"):
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

archive = root / twas["source_archive_path"]
entrypoint = root / twas["entrypoint_path"]
archive.parent.mkdir(parents=True, exist_ok=True)
if not archive.is_file():
    descriptor, name = tempfile.mkstemp(prefix="MetaXcan_v0.8.1.", suffix=".tar.gz", dir=archive.parent)
    os.close(descriptor); temporary = Path(name)
    try:
        subprocess.run(["curl", "-fL", "--retry", "3", "--output", str(temporary), "https://github.com/hakyimlab/MetaXcan/archive/refs/tags/v0.8.1.tar.gz"], check=True)
        if temporary.stat().st_size != twas["source_archive_bytes"] or digest(temporary) != twas["source_archive_sha256"]:
            raise SystemExit("ERROR: MetaXcan source archive differs from lock")
        os.replace(temporary, archive)
    finally:
        temporary.unlink(missing_ok=True)
if archive.stat().st_size != twas["source_archive_bytes"] or digest(archive) != twas["source_archive_sha256"]:
    raise SystemExit("ERROR: MetaXcan source archive differs from lock")
install = entrypoint.parents[1]
if not entrypoint.is_file():
    if install.exists():
        raise SystemExit("ERROR: partial MetaXcan installation exists")
    staging = install.with_name("MetaXcan.building")
    if staging.exists():
        raise SystemExit("ERROR: stale MetaXcan staging directory exists")
    with tempfile.TemporaryDirectory(dir=install.parent) as temporary_directory:
        with tarfile.open(archive, "r:gz") as handle:
            members = handle.getmembers()
            if any(member.name.startswith(("/", "../")) or "/../" in member.name for member in members):
                raise SystemExit("ERROR: unsafe MetaXcan archive member")
            handle.extractall(temporary_directory, filter="data")
        extracted = Path(temporary_directory) / "MetaXcan-0.8.1"
        shutil.move(extracted, staging)
    os.replace(staging, install)
if digest(entrypoint) != twas["entrypoint_sha256"]:
    raise SystemExit("ERROR: installed MetaXcan entrypoint differs from lock")
print("verified MetaXcan v0.8.1 source and entrypoint")

if fetch_models:
    source = twas["phi_model_source"]
    warning_threshold = float(policy["runtime"]["large_download_warning_gib"]) * 1024**3
    if int(source["expected_total_bytes"]) >= warning_threshold and not acknowledge_large:
        raise SystemExit(
            "ERROR: the locked PredictDB model family is a 2.920 GiB download; "
            "warn the user first, then re-run with --acknowledge-large-download after approval"
        )
    inventory_path = root / source["inventory_path"]
    inventory_lock_path = root / source["inventory_lock_path"]
    if not inventory_path.is_file() or not inventory_lock_path.is_file():
        raise SystemExit("ERROR: lock the result-free phi-model inventory before downloading models")
    inventory_lock = json.loads(inventory_lock_path.read_text(encoding="utf-8"))
    if inventory_lock.get("inventory_sha256") != digest(inventory_path) or inventory_lock.get("policy_sha256") != digest(Path(sys.argv[2])):
        raise SystemExit("ERROR: phi-model inventory differs from its lock")
    with inventory_path.open(encoding="utf-8", newline="") as handle:
        models = list(csv.DictReader(handle, delimiter="\t"))
    if inventory_lock.get("file_ids_in_locked_order") != [row["file_id"] for row in models]:
        raise SystemExit("ERROR: phi-model inventory family/order differs from lock")
    model_dir = root / source["install_dir"]
    model_dir.mkdir(parents=True, exist_ok=True)
    download_lock_path = root / source["download_lock_path"]
    existing_lock = json.loads(download_lock_path.read_text(encoding="utf-8")) if download_lock_path.is_file() else None
    print(f"NOTICE: fetching the locked phi-enabled PredictDB release requires {source['expected_total_bytes']:,} bytes (2.920 GiB)")
    installed = {}
    for row in models:
        name, expected_bytes = row["filename"], int(row["bytes"])
        final = model_dir / name
        expected_sha = (existing_lock or {}).get("files", {}).get(name, {}).get("sha256")
        if final.is_file() and final.stat().st_size == expected_bytes and (not expected_sha or digest(final) == expected_sha):
            installed[name] = {"file_id": row["file_id"], "bytes": expected_bytes, "sha256": digest(final)}
            print(f"verified {name}"); continue
        if final.exists():
            raise SystemExit(f"ERROR: existing phi-model file differs from lock: {final}")
        descriptor, temporary_name = tempfile.mkstemp(prefix=name + ".", dir=model_dir)
        os.close(descriptor); temporary = Path(temporary_name)
        try:
            subprocess.run(["curl", "-fL", "--retry", "3", "--output", str(temporary), row["download_url"]], check=True)
            if temporary.stat().st_size != expected_bytes:
                raise SystemExit(f"ERROR: phi-model file byte size differs from inventory: {name}")
            os.replace(temporary, final)
            installed[name] = {"file_id": row["file_id"], "bytes": expected_bytes, "sha256": digest(final)}
            print(f"installed {name}")
        finally:
            temporary.unlink(missing_ok=True)
    lock = {
        "schema_version": "atlas-v1.0-twas-phi-model-downloads.1",
        "locked_utc": (existing_lock or {}).get("locked_utc", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
        "inventory_sha256": digest(inventory_path), "inventory_lock_sha256": digest(inventory_lock_path),
        "policy_sha256": digest(Path(sys.argv[2])), "file_count": len(installed),
        "total_bytes": sum(record["bytes"] for record in installed.values()),
        "files": installed, "results_accessed": False, "claim_limit": policy["claim_limit"],
    }
    serialized = json.dumps(lock, indent=2, sort_keys=True) + "\n"
    if existing_lock is not None and download_lock_path.read_text(encoding="utf-8") != serialized:
        raise SystemExit("ERROR: phi-model download lock differs from current files")
    if existing_lock is None:
        download_lock_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_lock = download_lock_path.with_suffix(".tmp")
        temporary_lock.write_text(serialized, encoding="utf-8")
        os.replace(temporary_lock, download_lock_path)
PY

if [ "$RUNTIME" = true ]; then
  if [ -e "$ROOT/.molecular-env/python" ]; then
    echo "ERROR: molecular Python environment already exists; refusing overwrite" >&2
    exit 1
  fi
  conda create -y -p "$ROOT/.molecular-env/python" \
    python=3.9 numpy=1.26 pandas=2.2 scipy=1.12 pyliftover=0.4 \
    statsmodels=0.14 h5py=3.12 pyarrow=19.0 pip=25.0
  "$ROOT/.molecular-env/python/bin/python" -m pip install \
    cyvcf2==0.30.28 bgen-reader==4.0.8
fi

echo "TWAS_RESOURCES_CHECK_COMPLETE models=$MODELS runtime=$RUNTIME"
