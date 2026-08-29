#!/usr/bin/env python3
"""Freeze the exact six Phase-1 apparently-unreported validation candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_CANDIDATES = {
    "longsleep__parental_lifespan",
    "shortsleep__parental_lifespan",
    "sleep_apnea__healthspan",
    "sleep_apnea__parental_lifespan",
    "sleep_efficiency__frailty",
    "snoring__healthspan",
}
INPUTS = [
    Path("results/analysis/novel_connection_priorities.tsv"),
    Path("results/analysis/literature_novelty_audit.tsv"),
    Path("results/analysis/phase1_master_analysis.tsv"),
    Path("results/tables/rg_matrix.tsv"),
    Path("results/tables/h2_summary.tsv"),
    Path("config/analysis_panel.tsv"),
]
OUTPUT = Path("results/validation/candidate_validation_manifest.tsv")
LOCK = Path("results/validation/candidate_validation_manifest.lock.json")
FIELDS = [
    "candidate_id", "sleep_trait", "external_trait",
    "sleep_trait_label", "external_trait_label",
    "discovery_rg", "discovery_se", "discovery_p", "discovery_fdr",
    "discovery_gcov_intercept", "discovery_gcov_intercept_se",
    "discovery_SNP_overlap", "sleep_h2", "sleep_h2_se", "sleep_h2_z",
    "external_h2", "external_h2_se", "external_h2_z", "ancestry",
    "sleep_phenotype_definition", "external_phenotype_definition",
    "sleep_discovery_GWAS_source", "external_discovery_GWAS_source",
    "sleep_source_id", "external_source_id", "discovery_GWAS_sources",
    "sleep_PMID", "external_PMID", "PMIDs", "sleep_DOI", "external_DOI", "DOIs",
    "sleep_dataset_version", "external_dataset_version",
    "original_novelty_classification", "original_classification_confidence",
    "original_novelty_notes", "original_review_status", "original_claim_limit",
    "frozen_candidate_family_status",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"required input is missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unique_by(rows: list[dict[str, str]], field: str, label: str) -> dict[str, dict[str, str]]:
    output = {row[field]: row for row in rows}
    if len(output) != len(rows):
        fail(f"duplicate {field} in {label}")
    return output


def pair_id(row: dict[str, str]) -> str:
    return f"{row['sleep_trait']}__{row['external_trait']}"


def verify_lock() -> None:
    if not OUTPUT.is_file() or not LOCK.is_file():
        fail("candidate manifest or lock is missing")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    rows = read_tsv(OUTPUT)
    identities = [row["candidate_id"] for row in rows]
    if identities != lock.get("candidate_ids_in_locked_order") or set(identities) != EXPECTED_CANDIDATES:
        fail("candidate manifest membership/order differs from lock")
    if len(rows) != 6 or any(row["frozen_candidate_family_status"] != "FROZEN_BEFORE_VALIDATION_NO_REPLACEMENT" for row in rows):
        fail("candidate manifest freeze status is invalid")
    if sha256(OUTPUT) != lock.get("manifest_sha256"):
        fail("candidate manifest differs from its pre-validation lock")
    observed_inputs = {str(path): sha256(path) for path in INPUTS}
    if observed_inputs != lock.get("input_sha256"):
        fail("a canonical Phase-1 input differs from its pre-validation lock")
    if lock.get("validation_results_accessed_before_lock") is not False:
        fail("candidate lock does not certify pre-validation freezing")
    print(f"CANDIDATE_VALIDATION_MANIFEST_OK candidates={len(rows)} sha256={sha256(OUTPUT)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_lock()
        return
    if OUTPUT.exists() or LOCK.exists():
        fail("candidate family is already locked; use --verify instead of regenerating it")
    audit_rows = read_tsv(INPUTS[1])
    candidate_audit = {
        row["pair_id"]: row for row in audit_rows
        if row["novelty_classification"] == "APPARENTLY_NOVEL"
    }
    if set(candidate_audit) != EXPECTED_CANDIDATES or len(candidate_audit) != 6:
        fail(
            "canonical apparently-novel family differs from the expected six: "
            f"observed={sorted(candidate_audit)}"
        )

    master_rows = read_tsv(INPUTS[2])
    master = {pair_id(row): row for row in master_rows}
    if not EXPECTED_CANDIDATES.issubset(master):
        fail(f"Phase-1 master lacks candidates: {sorted(EXPECTED_CANDIDATES - set(master))}")
    if len([row for row in master_rows if pair_id(row) in EXPECTED_CANDIDATES]) != 6:
        fail("Phase-1 master contains duplicate candidate rows")

    rg = {
        f"{row['sleep_trait']}__{row['disease_trait']}": row
        for row in read_tsv(INPUTS[3])
    }
    h2 = unique_by(read_tsv(INPUTS[4]), "trait", str(INPUTS[4]))
    panel = unique_by(read_tsv(INPUTS[5]), "trait_id", str(INPUTS[5]))
    for identity in EXPECTED_CANDIDATES:
        if identity not in rg:
            fail(f"canonical rg matrix lacks candidate: {identity}")
        row = master[identity]
        matrix_row = rg[identity]
        for field_master, field_rg in (("rg", "rg"), ("se", "se"), ("p", "p"), ("fdr", "fdr")):
            if not math.isclose(float(row[field_master]), float(matrix_row[field_rg]), rel_tol=1e-12):
                fail(f"Phase-1 master and rg matrix disagree for {identity}/{field_master}")
        if not math.isclose(
            float(row["fdr_primary_372_sensitivity"]),
            float(matrix_row["fdr_primary_phase1"]),
            rel_tol=1e-12,
        ):
            fail(f"Phase-1 master and rg matrix disagree for {identity}/fdr_primary_372_sensitivity")
        for trait in (row["sleep_trait"], row["external_trait"]):
            if trait not in h2 or trait not in panel:
                fail(f"candidate trait lacks h2/panel metadata: {identity}/{trait}")
        if float(row["sleep_h2"]) != float(h2[row["sleep_trait"]]["h2"]):
            fail(f"sleep h2 mismatch for {identity}")
        if float(row["external_h2"]) != float(h2[row["external_trait"]]["h2"]):
            fail(f"external h2 mismatch for {identity}")

    output: list[dict[str, str]] = []
    for identity in sorted(EXPECTED_CANDIDATES):
        row = master[identity]
        audit = candidate_audit[identity]
        sleep = panel[row["sleep_trait"]]
        external = panel[row["external_trait"]]
        sleep_doi = sleep["doi"] if sleep["doi"] != "UNRESOLVED" else "10.1038/s41467-019-09576-1"
        output.append({
            "candidate_id": identity,
            "sleep_trait": row["sleep_trait"], "external_trait": row["external_trait"],
            "sleep_trait_label": row["sleep_trait_label"],
            "external_trait_label": row["external_trait_label"],
            "discovery_rg": row["rg"], "discovery_se": row["se"],
            "discovery_p": row["p"], "discovery_fdr": row["fdr"],
            "discovery_gcov_intercept": row["gcov_intercept"],
            "discovery_gcov_intercept_se": row["gcov_intercept_se"],
            "discovery_SNP_overlap": row["SNP_overlap"],
            "sleep_h2": row["sleep_h2"], "sleep_h2_se": row["sleep_h2_se"],
            "sleep_h2_z": row["sleep_h2_z"],
            "external_h2": row["external_h2"], "external_h2_se": row["external_h2_se"],
            "external_h2_z": row["external_h2_z"], "ancestry": row["ancestry"],
            "sleep_phenotype_definition": sleep["phenotype_definition"],
            "external_phenotype_definition": external["phenotype_definition"],
            "sleep_discovery_GWAS_source": sleep["source_note"],
            "external_discovery_GWAS_source": external["source_note"],
            "sleep_source_id": sleep["source_id"], "external_source_id": external["source_id"],
            "discovery_GWAS_sources": f"sleep={sleep['source_note']}; external={external['source_note']}",
            "sleep_PMID": sleep["pmid"], "external_PMID": external["pmid"],
            "PMIDs": f"{sleep['pmid']};{external['pmid']}",
            "sleep_DOI": sleep_doi, "external_DOI": external["doi"],
            "DOIs": f"{sleep_doi};{external['doi']}",
            "sleep_dataset_version": sleep["dataset_version"],
            "external_dataset_version": external["dataset_version"],
            "original_novelty_classification": audit["novelty_classification"],
            "original_classification_confidence": audit["classification_confidence"],
            "original_novelty_notes": audit["review_rationale"],
            "original_review_status": audit["review_status"],
            "original_claim_limit": audit["claim_limit"],
            "frozen_candidate_family_status": "FROZEN_BEFORE_VALIDATION_NO_REPLACEMENT",
        })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    lock = {
        "schema_version": "1.0.0",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "candidate_count": len(output),
        "candidate_ids_in_locked_order": [row["candidate_id"] for row in output],
        "candidate_selection_rule": "All and only Phase-1 literature_novelty_audit rows classified APPARENTLY_NOVEL",
        "candidate_replacement_policy": "FORBIDDEN; failed and unavailable candidates remain in the family",
        "validation_results_accessed_before_lock": False,
        "input_sha256": {str(path): sha256(path) for path in INPUTS},
        "manifest_path": str(OUTPUT),
        "manifest_sha256": sha256(OUTPUT),
    }
    LOCK.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"CANDIDATE_VALIDATION_MANIFEST_FROZEN candidates={len(output)} sha256={lock['manifest_sha256']}")


if __name__ == "__main__":
    main()
