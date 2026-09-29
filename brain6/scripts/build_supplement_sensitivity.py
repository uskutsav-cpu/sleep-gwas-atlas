#!/usr/bin/env python3
"""Publish the current interim headline-claim sensitivity ledger as Table S16."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "brain6/results/sensitivity/headline_claim_sensitivity.tsv"
SCRIPT = Path(__file__).resolve()
OUTPUT = ROOT / "brain6/results/supplement/table_S16_headline_claim_sensitivity.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")
REQUIRED = {
    "claim_id", "manuscript_claim", "evidence_status", "claim_scope", "evidence_sources",
    "alternative_gwas", "independent_replication", "sample_overlap_sensitivity",
    "phenotype_definition_sensitivity", "original_396_family_significance",
    "LAVA_UKB_v1_1", "HDL_L_local_rg", "local_h2_thresholds",
    "alternate_LD_reference", "MHC_exclusion", "fine_mapping_stability",
    "coloc_prior_sensitivity", "alternate_gene_mapping",
    "objective_vs_subjective_sleep", "ancestry_sensitivity",
    "overall_sensitivity_status", "interpretation_limit",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replace_if_changed(path: Path, data: bytes) -> None:
    if path.is_file() and path.read_bytes() == data:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main() -> dict[str, object]:
    source_bytes = SOURCE.read_bytes()
    text = source_bytes.decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows or not REQUIRED.issubset(rows[0]):
        raise ValueError("Sensitivity ledger has an unexpected schema")
    claim_ids = [row["claim_id"] for row in rows]
    if len(claim_ids) != len(set(claim_ids)):
        raise ValueError("Sensitivity ledger contains duplicate claim identities")
    if any(row["overall_sensitivity_status"] in {"PASS", "ROBUST"} for row in rows):
        raise ValueError("Interim sensitivity ledger must not promote a global robustness claim")

    replace_if_changed(OUTPUT, source_bytes)
    record: dict[str, object] = {
        "schema_version": 1,
        "status": "INTERIM_PARTIAL_SENSITIVITY_COVERAGE",
        "source_path": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256(source_bytes),
        "script_path": str(SCRIPT.relative_to(ROOT)),
        "script_sha256": sha256(SCRIPT.read_bytes()),
        "output_path": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": sha256(OUTPUT.read_bytes()),
        "claim_rows": len(rows),
        "unrun_or_incomplete_cells": sum(
            1 for row in rows for field, value in row.items()
            if field not in {"claim_id", "manuscript_claim", "evidence_status", "claim_scope",
                             "evidence_sources", "overall_sensitivity_status", "interpretation_limit"}
            and (value == "NOT_RUN" or "IN_PROGRESS" in value or value == "NOT_AVAILABLE")
        ),
        "interpretation": "Interim claim coverage only. Missing and in-progress sensitivities remain explicit; this table is not a complete robustness assessment or final results freeze.",
    }
    replace_if_changed(PROVENANCE, (json.dumps(record, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return record


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
