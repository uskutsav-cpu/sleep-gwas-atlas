#!/usr/bin/env python3
"""Validate the aggregate-only PGC MDD2025 global-rg sensitivity record."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "brain6/results/replication/pgc_mdd2025_insomnia_rg_sensitivity_v1.tsv"
PROVENANCE = TABLE.with_suffix(".provenance.json")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate() -> dict[str, object]:
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    builder = ROOT / provenance["builder"]
    run_path = Path(provenance["external_run_provenance"])
    run = json.loads(run_path.read_text(encoding="utf-8"))
    result = run["result"]
    log = Path(result["log_path"])
    with TABLE.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one aggregate sensitivity row, found {len(rows)}")
    row = rows[0]
    checks = {
        "status": provenance.get("status") == "PASS_AGGREGATE_ONLY_EXTERNAL_SENSITIVITY_RECORD",
        "external_run_status": run.get("status") == "PASS_LDSC_EXTERNAL_SENSITIVITY",
        "result_status": result.get("status") == "PASS_LDSC_EXTERNAL_SENSITIVITY",
        "result_values": (row.get("rg"), row.get("se"), row.get("p"),
                          row.get("cross_trait_intercept"), row.get("cross_trait_intercept_se")) ==
                         ("0.4771", "0.0228", "4.6535e-97", "0.0087", "0.0069"),
        "qualified_label": row.get("status") == "SENSITIVITY_ONLY_NOT_INDEPENDENT_REPLICATION",
        "not_ukbb": "noUKBB" in row.get("external_source", "") and "excludes UK Biobank" in row.get("overlap_status", ""),
        "source_md5": run["sources"]["mdd"].get("source_md5") == "29d3ce57cfce19ca28eb1643d9e2e428",
        "readme_md5": run["sources"]["mdd"].get("readme_md5") == "d743266ab7da04276517d3d4e033052a",
        "pgc_source_hash": provenance["source_checksums"].get("pgc_summary_stats_sha256") == run["sources"]["mdd"].get("source_sha256"),
        "run_receipt_hash": sha256(run_path) == provenance.get("external_run_provenance_sha256"),
        "rg_log_hash": sha256(log) == provenance.get("external_rg_log_sha256") == result.get("rg_log_sha256"),
        "builder_hash": sha256(builder) == provenance.get("builder_sha256"),
        "table_hash": sha256(TABLE) == provenance.get("output_sha256"),
        "external_data_only": not str(TABLE).startswith(str(Path(run["sources"]["mdd"]["source_path"]).parent)),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("PGC MDD2025 sensitivity integrity checks failed: " + ", ".join(failed))
    return {"status": "PASS_PGC_MDD2025_RG_SENSITIVITY", "checks": checks,
            "row_count": len(rows), "rg": row["rg"], "se": row["se"], "p": row["p"]}


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
