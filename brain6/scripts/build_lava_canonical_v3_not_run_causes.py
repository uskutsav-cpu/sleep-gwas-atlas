#!/usr/bin/env python3
"""Create an immutable, cell-level explanation table for canonical LAVA v3 NOT_RUNs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
ANALYSIS_ID = "brain6-lava-canonical-v3"
EXPECTED_CELLS = 17_465
EXPECTED_NOT_RUN = 3_720
CAUSES = {
    "LOW_LOCAL_H2_UNDERPOWERED": (
        "LOCAL_H2_SUPPORT_GATE",
        "The local heritability support criterion was not met; this is not a numerical failure.",
    ),
    "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS": (
        "INSUFFICIENT_SHARED_REFERENCE_VARIANTS",
        "Fewer than the minimum shared variants were available in the sealed reference.",
    ),
    "FEWER_THAN_MIN_K": (
        "INSUFFICIENT_MIN_K_COMPONENTS",
        "LAVA had fewer than the minimum components required to form the local test.",
    ),
}
COLUMNS = (
    "analysis_id", "run_id", "trait_id", "locus_id", "chromosome", "status",
    "reason", "cause_category", "cause_description", "interpretation",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def extract_cause_rows(rows: list[dict[str, str]], run_id: str = RUN_ID) -> list[dict[str, str]]:
    """Validate canonical cells and retain one auditable row per NOT_RUN cell."""
    unique: set[tuple[str, str]] = set()
    causes: list[dict[str, str]] = []
    for row in rows:
        trait, locus = row.get("phen", ""), row.get("locus_id", "")
        key = (trait, locus)
        if not trait or not locus or key in unique:
            raise ValueError(f"Canonical aggregate has a missing or duplicate cell identity: {key}")
        unique.add(key)
        if row.get("status") != "NOT_RUN":
            continue
        reason = row.get("reason", "")
        if reason not in CAUSES:
            raise ValueError(f"Unclassified canonical NOT_RUN reason: {reason!r}")
        category, description = CAUSES[reason]
        causes.append({
            "analysis_id": ANALYSIS_ID,
            "run_id": run_id,
            "trait_id": trait,
            "locus_id": locus,
            "chromosome": row.get("chromosome", ""),
            "status": "NOT_RUN",
            "reason": reason,
            "cause_category": category,
            "cause_description": description,
            "interpretation": "Not tested; this status is not evidence for a null association.",
        })
    return sorted(causes, key=lambda r: (r["trait_id"], int(r["locus_id"])))


def validate_source(run_dir: Path) -> tuple[dict, dict, Path, list[dict[str, str]]]:
    decision_path = run_dir / "canonical_family_decision.json"
    decision = read_json(decision_path)
    latest = read_json(run_dir / "latest_audit.json")
    aggregate_path = Path(decision["aggregate"]["path"])
    if not aggregate_path.is_absolute():
        aggregate_path = ROOT / aggregate_path
    if decision.get("analysis_id") != ANALYSIS_ID or decision.get("run_id") != RUN_ID:
        raise ValueError("Canonical decision is not the frozen v3 run")
    if decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or decision.get("promotion_permitted") is not False:
        raise ValueError("Canonical v3 decision no longer matches the frozen failed-QC disposition")
    if decision.get("aggregate_rows_verified") != EXPECTED_CELLS or decision["aggregate"].get("rows") != EXPECTED_CELLS:
        raise ValueError("Canonical aggregate is not the complete frozen 17,465-cell family")
    if sha256(aggregate_path) != decision["aggregate"].get("sha256"):
        raise ValueError("Canonical aggregate checksum differs from the immutable decision")
    if (latest.get("aggregate") != decision.get("aggregate") or
            latest.get("audit", {}).get("statuses") != decision.get("status_counts")):
        raise ValueError("Latest canonical audit disagrees with the immutable family decision")
    rows = read_rows(aggregate_path)
    causes = extract_cause_rows(rows)
    if len(rows) != EXPECTED_CELLS or len(causes) != EXPECTED_NOT_RUN:
        raise ValueError("Canonical aggregate dimensions differ from the frozen audited family")
    status_counts = Counter(row["status"] for row in rows)
    observed_statuses = {status: status_counts.get(status, 0) for status in ("TESTED", "NOT_RUN", "FAILED")}
    if observed_statuses != decision.get("status_counts"):
        raise ValueError("Canonical aggregate status counts differ from the decision")
    return decision, latest, aggregate_path, causes


def render_tsv(rows: list[dict[str, str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def write_immutable(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != content:
            raise FileExistsError(f"Refusing to replace existing diagnostic artifact: {path}")
        return
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(content)
    temporary.replace(path)


def build(run_dir: Path, output_path: Path, provenance_path: Path) -> dict:
    _decision, _latest, aggregate_path, causes = validate_source(run_dir)
    table_bytes = render_tsv(causes)
    script_path = Path(__file__).resolve()
    provenance = {
        "schema_version": 1,
        "status": "COMPLETE_DIAGNOSTIC_ONLY",
        "analysis_id": ANALYSIS_ID,
        "run_id": RUN_ID,
        "not_run_cells": len(causes),
        "cause_counts": dict(sorted(Counter(row["reason"] for row in causes).items())),
        "trait_counts": dict(sorted(Counter(row["trait_id"] for row in causes).items())),
        "canonical_decision_sha256": sha256(run_dir / "canonical_family_decision.json"),
        "canonical_aggregate_sha256": sha256(aggregate_path),
        "canonical_latest_audit_sha256": sha256(run_dir / "latest_audit.json"),
        "builder_script_sha256": sha256(script_path),
        "output_path": str(output_path),
        "output_sha256": hashlib.sha256(table_bytes).hexdigest(),
        "interpretation": "All 3,720 frozen NOT_RUN cells are classified from the checksum-bound canonical aggregate; no scientific threshold or run artifact is changed.",
    }
    provenance_bytes = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8")
    write_immutable(output_path, table_bytes)
    write_immutable(provenance_path, provenance_bytes)
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "work/lava-canonical-v3-production" / RUN_ID)
    parser.add_argument("--output", type=Path, default=ROOT / "brain6/results/lava/canonical_v3_not_run_cells.tsv")
    parser.add_argument("--provenance", type=Path, default=ROOT / "brain6/results/lava/canonical_v3_not_run_cells.provenance.json")
    args = parser.parse_args()
    result = build(args.run_dir.resolve(), args.output.resolve(), args.provenance.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
