#!/usr/bin/env python3
"""Summarize canonical v3 statuses and NOT_RUN causes by trait and locus."""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
from brain6.io import sha256
from build_lava_canonical_v3_not_run_causes import validate_source

RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
RUN_DIR = ROOT / "work/lava-canonical-v3-production" / RUN_ID
CAUSES = (
    "LOW_LOCAL_H2_UNDERPOWERED",
    "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS",
    "FEWER_THAN_MIN_K",
)
STATUS = ("TESTED", "NOT_RUN", "FAILED")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def render(rows: list[dict[str, Any]]) -> bytes:
    if not rows:
        raise ValueError("summary must contain rows")
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def immutable_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Refusing to replace diagnostic output: {path}")
        return
    temp = path.with_name(f".{path.name}.tmp")
    temp.write_bytes(payload)
    temp.replace(path)


def summarize(rows: list[dict[str, str]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    traits = sorted({r["phen"] for r in rows})
    loci = sorted({r["locus_id"] for r in rows}, key=int)
    if len({(r["phen"], r["locus_id"]) for r in rows}) != len(rows):
        raise ValueError("canonical aggregate contains duplicate trait-locus cells")
    if len(rows) != 17_465 or len(traits) != 7 or len(loci) != 2_495:
        raise ValueError("canonical aggregate dimensions differ from the frozen 7 x 2,495 family")

    by_trait: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    by_locus: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    chromosome: dict[str, str] = {}
    for row in rows:
        trait, locus, status = row["phen"], row["locus_id"], row["status"]
        if status not in STATUS:
            raise ValueError(f"unexpected canonical status: {status}")
        if locus in chromosome and chromosome[locus] != row["chromosome"]:
            raise ValueError(f"chromosome mismatch for locus {locus}")
        chromosome[locus] = row["chromosome"]
        by_trait[trait][status] += 1
        by_locus[locus][status] += 1
        if status == "NOT_RUN":
            if row["reason"] not in CAUSES:
                raise ValueError(f"unclassified NOT_RUN reason at {trait}/{locus}: {row['reason']}")
            by_trait[trait][row["reason"]] += 1
            by_locus[locus][row["reason"]] += 1
        elif row["reason"] in CAUSES:
            raise ValueError(f"non-NOT_RUN cell carries a NOT_RUN cause at {trait}/{locus}")

    trait_rows: list[dict[str, Any]] = []
    trait_counts: dict[str, dict[str, int]] = {}
    for trait in traits:
        c = by_trait[trait]
        if sum(c[s] for s in STATUS) != 2_495:
            raise ValueError(f"trait {trait} does not have all 2,495 locus cells")
        record: dict[str, Any] = {
            "analysis_id": "brain6-lava-canonical-v3", "run_id": RUN_ID,
            "trait_id": trait, "planned_cells": 2_495,
            "tested_cells": c["TESTED"], "not_run_cells": c["NOT_RUN"], "failed_cells": c["FAILED"],
        }
        for cause in CAUSES:
            record[cause.lower()] = c[cause]
        record["not_run_percent_of_planned"] = f"{100*c['NOT_RUN']/2495:.4f}"
        trait_rows.append(record)
        trait_counts[trait] = {"NOT_RUN": c["NOT_RUN"], **{cause: c[cause] for cause in CAUSES}}

    locus_rows: list[dict[str, Any]] = []
    for locus in loci:
        c = by_locus[locus]
        if sum(c[s] for s in STATUS) != 7:
            raise ValueError(f"locus {locus} does not have all seven trait cells")
        record = {
            "analysis_id": "brain6-lava-canonical-v3", "run_id": RUN_ID,
            "locus_id": locus, "chromosome": chromosome[locus], "planned_cells": 7,
            "tested_cells": c["TESTED"], "not_run_cells": c["NOT_RUN"], "failed_cells": c["FAILED"],
        }
        for cause in CAUSES:
            record[cause.lower()] = c[cause]
        locus_rows.append(record)

    summary = {
        "traits": len(traits), "loci": len(loci), "cells": len(rows),
        "tested_cells": sum(r["status"] == "TESTED" for r in rows),
        "not_run_cells": sum(r["status"] == "NOT_RUN" for r in rows),
        "failed_cells": sum(r["status"] == "FAILED" for r in rows),
        "cause_counts": {cause: sum(r["reason"] == cause and r["status"] == "NOT_RUN" for r in rows)
                         for cause in CAUSES},
        "not_run_loci": sum(by_locus[locus]["NOT_RUN"] > 0 for locus in loci),
        "loci_with_at_least_five_not_run": sum(by_locus[locus]["NOT_RUN"] >= 5 for locus in loci),
        "max_not_run_cells_at_one_locus": max(by_locus[locus]["NOT_RUN"] for locus in loci),
        "trait_counts": trait_counts,
    }
    if summary["not_run_cells"] != 3_720 or summary["failed_cells"] != 0:
        raise ValueError("status counts disagree with immutable canonical v3 decision")
    return trait_rows, locus_rows, summary


def build(run_dir: Path, trait_output: Path, locus_output: Path, provenance_path: Path) -> dict[str, Any]:
    decision, latest, aggregate_path, expected_causes = validate_source(run_dir)
    aggregate = read_tsv(aggregate_path)
    trait_rows, locus_rows, summary = summarize(aggregate)

    cell_path = ROOT / "brain6/results/lava/canonical_v3_not_run_cells.tsv"
    cell_prov_path = cell_path.with_suffix(".provenance.json")
    cell_prov = json.loads(cell_prov_path.read_text(encoding="utf-8"))
    if cell_prov.get("output_sha256") != sha256(cell_path) or cell_prov.get("not_run_cells") != 3_720:
        raise ValueError("source-bound canonical NOT_RUN table/provenance is invalid")
    if read_tsv(cell_path) != expected_causes:
        raise ValueError("cell-level NOT_RUN table differs from the canonical receipt-derived rows")
    if latest.get("audit", {}).get("statuses") != decision.get("status_counts"):
        raise ValueError("latest canonical audit does not match the immutable decision")

    trait_bytes, locus_bytes = render(trait_rows), render(locus_rows)
    immutable_write(trait_output, trait_bytes)
    immutable_write(locus_output, locus_bytes)
    script = Path(__file__).resolve()
    provenance = {
        "schema_version": 1,
        "status": "COMPLETE_DIAGNOSTIC_ONLY",
        "analysis_id": "brain6-lava-canonical-v3-cause-breakdown-v1",
        "run_id": RUN_ID,
        "interpretation": "Descriptive breakdown of the frozen canonical v3 receipts; these counts explain test eligibility/status, not null association evidence or numerical failures. The immutable FAILED_QC_NOT_PROMOTED decision and 5% threshold are unchanged.",
        "source_hashes": {
            "canonical_family_decision": sha256(run_dir / "canonical_family_decision.json"),
            "canonical_aggregate": sha256(aggregate_path),
            "canonical_latest_audit": sha256(run_dir / "latest_audit.json"),
            "cell_level_not_run_table": sha256(cell_path),
            "cell_level_not_run_provenance": sha256(cell_prov_path),
        },
        "builder_path": str(script.relative_to(ROOT)),
        "builder_sha256": sha256(script),
        "summary": summary,
        "trait_table": {"path": str(trait_output.relative_to(ROOT)), "rows": len(trait_rows),
                        "sha256": hashlib.sha256(trait_bytes).hexdigest()},
        "locus_table": {"path": str(locus_output.relative_to(ROOT)), "rows": len(locus_rows),
                        "sha256": hashlib.sha256(locus_bytes).hexdigest()},
    }
    immutable_write(provenance_path, (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode())
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=RUN_DIR)
    parser.add_argument("--trait-output", type=Path, default=ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv")
    parser.add_argument("--locus-output", type=Path, default=ROOT / "brain6/results/lava/canonical_v3_status_by_locus_v1.tsv")
    parser.add_argument("--provenance", type=Path, default=ROOT / "brain6/results/lava/canonical_v3_cause_breakdown_v2.provenance.json")
    args = parser.parse_args()
    result = build(args.run_dir.resolve(), args.trait_output.resolve(), args.locus_output.resolve(), args.provenance.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
