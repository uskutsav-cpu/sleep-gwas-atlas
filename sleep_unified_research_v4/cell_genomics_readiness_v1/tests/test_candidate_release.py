#!/usr/bin/env python3
"""Check release family identity, row counts, mappings, and payload exclusions."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PACKAGE = Path(__file__).resolve().parents[1]
RELEASE = PACKAGE / "release_candidate"
EXPECTED = {
    "core_rg.tsv": 396,
    "extension_rg.tsv": 1200,
    "external_validation_rg.tsv": 41,
    "heritability.tsv": 158,
    "sensitivities.tsv": 62,
    "sleep_constructs.tsv": 12,
    "validation_217_class_counts.tsv": 4,
    "outcome_phenotypes.tsv": 158,
}


def read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    manifest = json.loads((RELEASE / "release_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "sleep-atlas-candidate-release/1.0"
    assert manifest["release_status"] == "CANDIDATE_FOR_RIGHTS_AND_METADATA_REVIEW"
    assert manifest["counts"] == {
        "heritability": 158, "rg_core": 396, "rg_extension": 1200,
        "rg_external_estimates": 41, "sensitivities": 62, "sleep_constructs": 12,
        "outcome_phenotypes": 158,
    }
    assert len(manifest["files"]) == 10
    assert digest(ROOT / manifest["builder"]["path"]) == manifest["builder"]["sha256"]
    for source in manifest["input_files"]:
        assert digest(ROOT / source["path"]) == source["sha256"], source["path"]

    schema = json.loads((RELEASE / "schema.json").read_text(encoding="utf-8"))
    for filename, record in schema["files"].items():
        fields, _ = read(RELEASE / filename)
        assert len(fields) == len(record["fields"]), filename
        assert fields == record["fields"], filename
    for filename in ("SCHEMA.md", "schema.json"):
        assert digest(RELEASE / filename) == manifest["files"][filename]["sha256"]

    all_text = []
    for filename, expected_rows in EXPECTED.items():
        fields, rows = read(RELEASE / filename)
        assert len(rows) == expected_rows, (filename, len(rows))
        assert fields and len(fields) == len(set(fields)), filename
        assert manifest["files"][filename]["rows"] == expected_rows
        assert digest(RELEASE / filename) == manifest["files"][filename]["sha256"]
        assert not ({"raw_path", "evidence_paths", "original_log", "native_full_precision_path"} & set(fields)), filename
        all_text.append((RELEASE / filename).read_text(encoding="utf-8"))

    for name, n, family in [("core_rg.tsv", 396, "core"), ("extension_rg.tsv", 1200, "extension")]:
        _, data = read(RELEASE / name)
        assert {r["analysis_family"] for r in data} == {family}
        assert {r["family_size"] for r in data} == {str(n)}
        assert all(r["adjustment_method"] == f"frozen_BH_{n}" for r in data)
        assert all(r["sleep_source_id"] for r in data)
        assert all(r["outcome_source_id"] for r in data)

    _, external = read(RELEASE / "external_validation_rg.tsv")
    assert all(r["analysis_family"] == "validation" and r["family_size"] == "217" for r in external)
    assert all(r["adjustment_method"] == "Bonferroni_0.05_over_217" for r in external)
    assert all(r["fully_independent_two_trait_replication"] == "False" for r in external)
    _, outcome_catalog = read(RELEASE / "outcome_phenotypes.tsv")
    catalog = {(r["analysis_family"], r["trait_id"]): r for r in outcome_catalog}
    extension_ids = {r["trait_id"] for r in outcome_catalog if r["analysis_family"] == "extension"}
    for r in external:
        assert r["discovery_outcome_trait_id"] in extension_ids
        meta = catalog[("validation", r["outcome_trait"])]
        assert meta["source_id"] == r["outcome_source_id"]
        assert meta["phenotype_definition"] and meta["label"]
        assert meta["phenotype_type"] == "UNKNOWN"
    for family, filename in [("core", "core_rg.tsv"), ("extension", "extension_rg.tsv")]:
        _, estimates = read(RELEASE / filename)
        for r in estimates:
            meta = catalog[(family, r["outcome_trait"])]
            assert meta["source_id"] == r["outcome_source_id"]
    _, classes = read(RELEASE / "validation_217_class_counts.tsv")
    assert sum(int(r["n_pairs"]) for r in classes) == 217
    class_counts = {r["classification"]: int(r["n_pairs"]) for r in classes}
    assert sorted(class_counts.values()) == [17, 18, 23, 159]

    _, sensitivities = read(RELEASE / "sensitivities.tsv")
    assert all(not r["difference_p_value"] for r in sensitivities)
    assert all("NOT_CALIBRATED" in r["delta_uncertainty"] for r in sensitivities)
    payload = "\n".join(all_text)
    assert not re.search(r"/Users/|/Volumes/|file://|X-Amz-Signature|token=", payload, re.I)
    print(f"CANDIDATE_RELEASE_INTEGRITY_OK files=10 rows={sum(EXPECTED.values())} validation_classes=217 independent_two_trait_replications=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
