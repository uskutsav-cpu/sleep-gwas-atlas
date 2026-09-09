"""Shared guard for source-specific production helpers."""
import csv
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "analysis_panel.tsv"
LOCK = ROOT / "config" / "analysis_panel.lock.json"


def require_locked_traits(expected_sources: dict[str, str]) -> None:
    """Require exact trait/source mappings in the validated production panel."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "00_validate_panel.py"),
            "--manifest",
            str(MANIFEST),
            "--lock",
            str(LOCK),
        ],
        cwd=ROOT,
        check=True,
    )
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = {row["trait_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    for trait_id, source_id in expected_sources.items():
        row = rows.get(trait_id)
        if row is None:
            raise SystemExit(f"ERROR: {trait_id} is outside the locked analysis panel")
        if row["source_id"] != source_id:
            raise SystemExit(
                f"ERROR: {trait_id} source changed: expected {source_id}, got {row['source_id']}"
            )
        if row["source_status"] != "SOURCE_VERIFIED":
            raise SystemExit(
                f"ERROR: {trait_id} has source_status={row['source_status']}"
            )
