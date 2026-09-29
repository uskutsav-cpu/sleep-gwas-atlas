#!/usr/bin/env python3
"""Build a traceable table of published, overlap-limited rg context."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6"
SOURCE = OUT / "manifests/published_context_literature.json"
DISCOVERY = OUT / "results/global/brain6_72_locked.tsv"
TABLE = OUT / "results/supplement/table_S18_published_rg_context.tsv"
PROVENANCE = OUT / "results/supplement/table_S18_published_rg_context.provenance.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    with DISCOVERY.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    matches = [row for row in rows if row["sleep_trait"] == "longsleep" and row["brain_disorder"] == "scz"]
    if len(matches) != 1:
        raise ValueError(f"Expected one locked longsleep/SCZ discovery row, found {len(matches)}")
    discovery = matches[0]
    published = source["source_reported_result"]
    if source["record_id"] != "austin_zimmerman_2023_longsleep_scz":
        raise ValueError("Unexpected literature record")
    difference = float(published["rg"]) - float(discovery["rg"])
    row = {
        "record_id": source["record_id"],
        "citation": source["citation"],
        "doi": source["doi"],
        "source_url": source["source_url"],
        "accessed_date": source["accessed_date"],
        "sleep_trait": "longsleep",
        "brain_disorder": "scz",
        "published_sleep_definition": published["sleep_case_definition"],
        "published_sleep_ancestry": published["sleep_ancestry"],
        "published_sleep_n": published["sleep_n"],
        "published_sleep_cohorts": published["sleep_cohorts"],
        "published_disorder_source": published["disorder_summary_source"],
        "published_rg": published["rg"],
        "published_se": published["se"],
        "published_p": published["p"],
        "brain6_sleep_definition": source["brain6_comparison"]["sleep_case_definition"],
        "brain6_sleep_source": source["brain6_comparison"]["sleep_source"],
        "brain6_disorder_source": source["brain6_comparison"]["disorder_source"],
        "brain6_rg": discovery["rg"],
        "brain6_se": discovery["se"],
        "brain6_p": discovery["p"],
        "published_minus_brain6_rg_point_difference": f"{difference:.4f}",
        "direction_concordant": str(float(published["rg"]) * float(discovery["rg"]) > 0).upper(),
        "classification": source["brain6_comparison"]["classification"],
        "admission_status": source["brain6_comparison"]["admission_status"],
        "overlap_status": source["brain6_comparison"]["overlap_status"],
        "interpretation": source["brain6_comparison"]["interpretation"],
    }
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    with TABLE.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    provenance = {
        "status": "PASS_CONTEXT_ONLY",
        "source_manifest": str(SOURCE.relative_to(ROOT)),
        "source_manifest_sha256": sha256(SOURCE),
        "discovery_table": str(DISCOVERY.relative_to(ROOT)),
        "discovery_table_sha256": sha256(DISCOVERY),
        "output_table": str(TABLE.relative_to(ROOT)),
        "output_table_sha256": sha256(TABLE),
        "data_rows": 1,
        "statistical_reanalysis_performed": False,
        "replication_admitted": False,
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"path": str(TABLE), "status": provenance["status"], "rows": 1}, indent=2))


if __name__ == "__main__":
    main()
