"""Bind the 35 inherited Brain6 discoveries to the dated parent-atlas audit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
GLOBAL = ROOT / "brain6/results/global/brain6_72_locked.tsv"
AUDIT = ROOT / "results/analysis/literature_novelty_audit.tsv"
SEARCH_QUERIES = ROOT / "results/analysis/literature_search_queries.tsv"
SOURCE_EVIDENCE = ROOT / "results/analysis/published_rg_source_evidence.tsv"
OUTPUT = ROOT / "brain6/results/novelty/brain6_significant_pair_novelty_crosswalk.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")

OUTPUT_FIELDS = [
    "sleep_trait", "brain_disorder", "brain_disorder_label", "rg", "se", "p",
    "original_396_family_q", "significance_under_original_396_family",
    "novelty_classification", "classification_confidence", "direct_rg_found",
    "direct_evidence_count", "numeric_direct_evidence_count", "representative_pmid",
    "representative_doi", "representative_title", "representative_published_rg",
    "phenotype_definition_match", "dataset_relation", "comparison_independence",
    "evidence_pmids", "evidence_dois", "evidence_source_ids", "search_as_of_date",
    "search_status", "direct_rg_search_result", "claim_limit",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build_rows() -> tuple[list[dict[str, str]], dict[str, Any]]:
    global_rows = read_tsv(GLOBAL)
    audit_rows = read_tsv(AUDIT)
    if len(global_rows) != 72:
        raise ValueError(f"Expected 72 locked Brain6 pairs, found {len(global_rows)}")
    audit_by_pair: dict[tuple[str, str], dict[str, str]] = {}
    for row in audit_rows:
        key = (row["sleep_trait"], row["external_trait"])
        if key in audit_by_pair:
            raise ValueError(f"Duplicate pair in parent novelty audit: {key}")
        audit_by_pair[key] = row

    discoveries = [row for row in global_rows
                   if row["significance_under_original_396_family"].lower() == "true"]
    out: list[dict[str, str]] = []
    for row in discoveries:
        key = (row["sleep_trait"], row["brain_disorder"])
        evidence = audit_by_pair.get(key)
        if evidence is None:
            raise ValueError(f"No pair-level literature audit for inherited discovery {key}")
        out.append({
            "sleep_trait": row["sleep_trait"],
            "brain_disorder": row["brain_disorder"],
            "brain_disorder_label": row["brain_disorder_label"],
            "rg": row["rg"],
            "se": row["se"],
            "p": row["p"],
            "original_396_family_q": row["original_BH_FDR_q"],
            "significance_under_original_396_family": row["significance_under_original_396_family"],
            **{field: evidence.get(field, "") for field in OUTPUT_FIELDS[8:]},
        })
    if len(out) != 35:
        raise ValueError(f"Expected 35 inherited significant pairs, found {len(out)}")
    classes = Counter(row["novelty_classification"] for row in out)
    dates = {row["search_as_of_date"] for row in out}
    if len(dates) != 1:
        raise ValueError(f"Parent audit search dates differ: {sorted(dates)}")
    summary = {
        "global_pair_count": len(global_rows),
        "inherited_significant_pair_count": len(out),
        "parent_audit_search_as_of_date": next(iter(dates)),
        "classification_counts": dict(sorted(classes.items())),
        "apparent_novelty_count": classes.get("APPARENTLY_NOVEL", 0),
        "direct_rg_count": sum(row["direct_rg_found"].lower() == "true" for row in out),
        "direct_rg_independence_counts": dict(sorted(Counter(
            row["comparison_independence"] for row in out).items())),
        "interpretation": (
            "Dated parent-atlas literature crosswalk, not a fresh current-date search; "
            "direct rg comparisons frequently reuse or overlap discovery GWAS and are not "
            "independent replication. No first-ever claim is supported."
        ),
    }
    return out, summary


def write_outputs() -> dict[str, Any]:
    rows, summary = build_rows()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".tsv.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(OUTPUT)
    inputs = [GLOBAL, AUDIT, SEARCH_QUERIES, SOURCE_EVIDENCE]
    provenance = {
        "schema_version": 1,
        "status": "PASS_DATED_LITERATURE_CROSSWALK_NO_NOVELTY_CLAIM",
        "scope": "All and only 35 Brain6 pairs significant under the inherited original 396-test family.",
        "builder": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "output": {"path": str(OUTPUT.relative_to(ROOT)), "sha256": sha256(OUTPUT), "rows": len(rows)},
        "summary": summary,
        "limitations": [
            "The parent audit is dated 2026-08-28 and was not refreshed in this build.",
            "Search-based absence is not proof that no prior report exists.",
            "Prior direct rg comparisons are mostly not fully independent of discovery GWAS.",
            "This crosswalk does not validate the original 396-family results or imply causality.",
        ],
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return provenance


def validate() -> dict[str, Any]:
    if not OUTPUT.is_file() or not PROVENANCE.is_file():
        raise ValueError("Crosswalk or provenance is missing")
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    rows = read_tsv(OUTPUT)
    expected_rows, expected_summary = build_rows()
    if rows != expected_rows:
        raise ValueError("Crosswalk differs from current locked map or parent audit")
    if (provenance.get("status") != "PASS_DATED_LITERATURE_CROSSWALK_NO_NOVELTY_CLAIM" or
            set(provenance.get("inputs", {})) != {
                str(path.relative_to(ROOT)) for path in (GLOBAL, AUDIT, SEARCH_QUERIES, SOURCE_EVIDENCE)
            } or
            provenance.get("scope") != (
                "All and only 35 Brain6 pairs significant under the inherited original 396-test family."
            ) or
            provenance.get("builder_sha256") != sha256(Path(__file__).resolve()) or
            provenance.get("output", {}).get("sha256") != sha256(OUTPUT) or
            provenance.get("summary") != expected_summary):
        raise ValueError("Crosswalk provenance or summary does not validate")
    for name, digest in provenance.get("inputs", {}).items():
        if sha256(ROOT / name) != digest:
            raise ValueError(f"Crosswalk input hash mismatch: {name}")
    return {"status": "PASS_DATED_LITERATURE_CROSSWALK_NO_NOVELTY_CLAIM", **expected_summary,
            "output_sha256": sha256(OUTPUT)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    print(json.dumps(validate() if args.validate else write_outputs(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
