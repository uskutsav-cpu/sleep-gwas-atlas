#!/usr/bin/env python3
"""Validate source-snapshotted mechanistic evidence against the locked search plan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


EVIDENCE_FIELDS = [
    "evidence_id", "search_task_id", "comparison_id", "pair_id", "locus_id", "signal1", "signal2",
    "search_task", "search_outcome", "evidence_class", "source_id", "exact_dataset_release",
    "accession", "source_url", "source_access_date", "source_license_or_terms",
    "snapshot_path", "snapshot_sha256", "primary_citation_title", "primary_citation_DOI",
    "primary_citation_PMID", "species", "biological_context", "variant_id", "regulatory_element_id",
    "gene_id", "gene_symbol", "cell_type_or_state", "pathway_id", "pathway_name", "perturbation",
    "endpoint", "effect_direction", "effect_size", "SE", "P", "FDR", "posterior_probability",
    "colocalization_interpretation", "supports_edge_from", "supports_edge_to", "contradiction_status",
    "evidence_summary", "caveat", "claim_limit", "curator", "curation_date",
]
TASK_CLASS = {
    "MOLECULAR_QTL": "COLOCALIZED_MOLECULAR_QTL",
    "VARIANT_REGULATORY_ELEMENT": "FINE_MAPPED_REGULATORY_EVIDENCE",
    "REGULATORY_ELEMENT_TARGET_GENE": "FINE_MAPPED_REGULATORY_EVIDENCE",
    "TISSUE_CELL_EXPRESSION": "TISSUE_CELL_EXPRESSION_OR_ENRICHMENT",
    "CELL_STATE_ACCESSIBILITY": "TISSUE_CELL_EXPRESSION_OR_ENRICHMENT",
    "PATHWAY": "PATHWAY_OR_NETWORK_ANNOTATION",
    "MODEL_SYSTEM_PERTURBATION": "MODEL_SYSTEM_OR_PERTURBATION_EVIDENCE",
}
EDGE_NODES = {
    "SLEEP_TRAIT", "SHARED_SIGNAL_MODEL", "CREDIBLE_VARIANT", "REGULATORY_ELEMENT",
    "TARGET_GENE", "CELL_TYPE_OR_STATE", "PATHWAY", "EXTERNAL_PHENOTYPE",
}
OUTCOMES = {"SUPPORTED", "NO_EVIDENCE_FOUND", "ACCESS_BLOCKED", "NOT_APPLICABLE"}
CONTRADICTIONS = {"NONE", "DIRECTION_CONFLICT", "CONTEXT_MISMATCH", "SOURCE_CONFLICT", "UNRESOLVED"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def present(value: str) -> bool:
    return value not in {"", "NA", ".", "PENDING"}


def optional_probability(value: str, field: str, identity: str) -> None:
    if not present(value):
        return
    try:
        number = float(value)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid {field}: {identity}") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise SystemExit(f"ERROR: invalid {field}: {identity}")


def require_supported_entities(row: dict[str, str]) -> None:
    task, identity = row["search_task"], row["evidence_id"]
    requirements = {
        "MOLECULAR_QTL": ("gene_id", "gene_symbol"),
        "VARIANT_REGULATORY_ELEMENT": ("variant_id", "regulatory_element_id"),
        "REGULATORY_ELEMENT_TARGET_GENE": ("regulatory_element_id", "gene_id", "gene_symbol"),
        "TISSUE_CELL_EXPRESSION": ("gene_id", "gene_symbol", "cell_type_or_state"),
        "CELL_STATE_ACCESSIBILITY": ("regulatory_element_id", "cell_type_or_state"),
        "PATHWAY": ("gene_id", "gene_symbol", "pathway_id", "pathway_name"),
        "MODEL_SYSTEM_PERTURBATION": ("gene_id", "gene_symbol", "perturbation", "endpoint"),
    }
    for field in requirements[task]:
        if not present(row[field]):
            raise SystemExit(f"ERROR: supported {task} evidence lacks {field}: {identity}")
    if row["supports_edge_from"] not in EDGE_NODES or row["supports_edge_to"] not in EDGE_NODES or row["supports_edge_from"] == row["supports_edge_to"]:
        raise SystemExit(f"ERROR: supported evidence lacks a valid directed chain edge: {identity}")
    if task == "MOLECULAR_QTL" and row["colocalization_interpretation"] != "SHARED_SIGNAL_MODEL_SUPPORTED":
        raise SystemExit(f"ERROR: molecular-QTL evidence is not formally colocalized: {identity}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_search_queue.tsv"))
    parser.add_argument("--family-lock", type=Path, default=Path("discovery_extension/config/mechanistic_signal_family.lock.json"))
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/mechanistic_sources.tsv"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/mechanistic_annotation_contract.json"))
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_evidence.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/mechanistic_evidence_validation.json"))
    args = parser.parse_args()

    _, queue = read_tsv(args.queue)
    lock = json.loads(args.family_lock.read_text(encoding="utf-8"))
    sources = {row["source_id"]: row for row in read_tsv(args.sources)[1]}
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    fields, evidence = read_tsv(args.evidence)
    if fields != EVIDENCE_FIELDS:
        raise SystemExit("ERROR: mechanistic evidence fields differ from the locked schema")
    if sha256(args.queue) != lock.get("queue_template_sha256"):
        raise SystemExit("ERROR: mechanistic search queue differs from its pre-search lock")
    if sha256(args.sources) != lock.get("sources_sha256") or sha256(args.contract) != lock.get("contract_sha256"):
        raise SystemExit("ERROR: mechanistic source or analysis contract differs from the signal-family lock")
    task_ids = [row["search_task_id"] for row in queue]
    if task_ids != lock.get("search_task_ids_in_locked_order") or len(task_ids) != lock.get("search_task_count"):
        raise SystemExit("ERROR: mechanistic task family/order differs from lock")
    if len({row["evidence_id"] for row in evidence}) != len(evidence):
        raise SystemExit("ERROR: duplicate mechanistic evidence identifier")
    queue_by_id = {row["search_task_id"]: row for row in queue}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence:
        identity = row["evidence_id"]
        if row["search_task_id"] not in queue_by_id:
            raise SystemExit(f"ERROR: evidence outside locked search family: {identity}")
        task = queue_by_id[row["search_task_id"]]
        for field in ("comparison_id", "pair_id", "locus_id", "signal1", "signal2", "search_task"):
            expected_field = field
            if row[field] != task[expected_field]:
                raise SystemExit(f"ERROR: evidence identity differs from locked task: {identity}/{field}")
        if task["results_accessed_before_plan_lock"] != "NO":
            raise SystemExit(f"ERROR: source results preceded mechanistic search lock: {identity}")
        if row["search_outcome"] not in OUTCOMES or row["contradiction_status"] not in CONTRADICTIONS:
            raise SystemExit(f"ERROR: invalid search outcome or contradiction status: {identity}")
        allowed_sources = set(task["allowed_source_ids"].split(";"))
        if row["source_id"] not in allowed_sources or row["source_id"] not in sources:
            raise SystemExit(f"ERROR: evidence uses an unplanned source: {identity}/{row['source_id']}")
        for field in (
            "exact_dataset_release", "accession", "source_url", "source_access_date",
            "source_license_or_terms", "snapshot_path", "snapshot_sha256", "species",
            "biological_context", "evidence_summary", "caveat", "claim_limit", "curator", "curation_date",
        ):
            if not present(row[field]):
                raise SystemExit(f"ERROR: mechanistic evidence lacks {field}: {identity}")
        if not row["source_url"].startswith(("https://", "synthetic://")):
            raise SystemExit(f"ERROR: invalid mechanistic source URL: {identity}")
        snapshot = Path(row["snapshot_path"])
        if not snapshot.is_file() or len(row["snapshot_sha256"]) != 64 or sha256(snapshot) != row["snapshot_sha256"]:
            raise SystemExit(f"ERROR: mechanistic source snapshot missing or checksum-mismatched: {identity}")
        optional_probability(row["P"], "P", identity)
        optional_probability(row["FDR"], "FDR", identity)
        optional_probability(row["posterior_probability"], "posterior_probability", identity)
        if row["search_outcome"] == "SUPPORTED":
            if row["evidence_class"] != TASK_CLASS[row["search_task"]]:
                raise SystemExit(f"ERROR: supported evidence class differs from task: {identity}")
            require_supported_entities(row)
            synthetic_fixture = row["curator"] == "synthetic_test" and row["source_url"].startswith("synthetic://")
            if not present(row["primary_citation_title"]) or (not present(row["primary_citation_DOI"]) and not present(row["primary_citation_PMID"])):
                if not synthetic_fixture:
                    raise SystemExit(f"ERROR: supported evidence lacks a primary-source citation: {identity}")
        else:
            if row["evidence_class"] != "NO_SUPPORTED_EVIDENCE":
                raise SystemExit(f"ERROR: nonsupporting search outcome has a positive evidence class: {identity}")
        grouped[row["search_task_id"]].append(row)
    if set(grouped) != set(task_ids):
        raise SystemExit("ERROR: mechanistic evidence does not cover every locked search task")
    for task_id, rows in grouped.items():
        outcomes = {row["search_outcome"] for row in rows}
        if "SUPPORTED" in outcomes and outcomes != {"SUPPORTED"}:
            raise SystemExit(f"ERROR: search task mixes supported and terminal no-evidence outcomes: {task_id}")
        if "SUPPORTED" not in outcomes and len(rows) != 1:
            raise SystemExit(f"ERROR: terminal no-evidence search outcome must have exactly one row: {task_id}")

    order = {task_id: index for index, task_id in enumerate(task_ids)}
    evidence.sort(key=lambda row: (order[row["search_task_id"]], row["evidence_id"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=EVIDENCE_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence)
    counts = Counter(row["search_outcome"] for row in evidence)
    provenance = {
        "schema_version": "1.0.0",
        "validated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "family_lock_sha256": sha256(args.family_lock),
        "queue_sha256": sha256(args.queue),
        "source_registry_sha256": sha256(args.sources),
        "contract_sha256": sha256(args.contract),
        "input_evidence_sha256": sha256(args.evidence),
        "search_task_count": len(task_ids),
        "evidence_row_count": len(evidence),
        "search_outcome_counts": dict(sorted(counts.items())),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "warning": contract["claim_limit"],
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"MECHANISTIC_EVIDENCE_VALID tasks={len(task_ids)} rows={len(evidence)} outcomes={dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
