#!/usr/bin/env python3
"""Fail closed on incomplete, unsupported, or overclaimed novelty decisions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALLOWED_CLASSES = {
    "KNOWN_REPLICATION", "KNOWN_BUT_NEW_DATASET", "PARTIAL_EXTENSION",
    "NO_DIRECT_RG_FOUND", "APPARENTLY_NOVEL", "UNCERTAIN",
}
ALLOWED_STRENGTHS = {"STRONG", "MODERATE", "WEAK"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rg", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
    )
    parser.add_argument(
        "--audit", type=Path,
        default=Path("discovery_extension/results/novelty/extension_novelty_audit.tsv"),
    )
    args = parser.parse_args()
    hits = {
        (row["sleep_trait"], row["extension_trait_id"])
        for row in read_tsv(args.rg) if float(row["extension_fdr"]) < 0.05
    }
    audit = read_tsv(args.audit)
    observed = {(row["sleep_trait"], row["extension_trait_id"]) for row in audit}
    if observed != hits or len(observed) != len(audit):
        raise SystemExit(
            f"ERROR: novelty audit family mismatch hits={len(hits)} audit={len(audit)} "
            f"missing={sorted(hits-observed)[:5]} extra={sorted(observed-hits)[:5]}"
        )
    for row in audit:
        pair = row["pair_id"]
        if row["audit_status"] != "COMPLETE":
            raise SystemExit(f"ERROR: incomplete pair-level audit: {pair}")
        if row["novelty_class"] not in ALLOWED_CLASSES:
            raise SystemExit(f"ERROR: invalid novelty class for {pair}")
        if row["novelty_strength"] not in ALLOWED_STRENGTHS:
            raise SystemExit(f"ERROR: invalid novelty strength for {pair}")
        for field in (
            "exact_prior_rg_found", "closest_prior_result", "prior_method", "prior_effect",
            "prior_publication", "prior_DOI", "prior_PMID", "search_databases",
            "search_queries_used", "search_date", "evidence_PMIDs_DOIs_URLs",
            "decision_rationale", "reviewer_notes", "reviewer",
        ):
            if not row[field].strip() or row[field].startswith("PENDING"):
                raise SystemExit(f"ERROR: completed audit lacks {field}: {pair}")
        if row["novelty_class"] == "APPARENTLY_NOVEL":
            targeted_queries = [
                query.strip() for query in row["search_queries_used"].split(" || ")
                if query.strip()
            ]
            if len(targeted_queries) < 3:
                raise SystemExit(f"ERROR: APPARENTLY_NOVEL lacks several targeted searches: {pair}")
            if row["exact_prior_rg_found"] != "NO" or row["direct_prior_same_pair"] != "NO":
                raise SystemExit(f"ERROR: APPARENTLY_NOVEL conflicts with direct prior rg evidence: {pair}")
        if row["novelty_strength"] == "STRONG":
            if row["novelty_class"] != "APPARENTLY_NOVEL" or row["exact_prior_rg_found"] != "NO":
                raise SystemExit(f"ERROR: STRONG novelty lacks APPARENTLY_NOVEL/no-direct-prior evidence: {pair}")
            if row["independent_replication_status"] != "INDEPENDENT_REPLICATION_PASS":
                raise SystemExit(f"ERROR: STRONG novelty lacks successful independent replication: {pair}")
    print(f"PAIR_NOVELTY_AUDIT_VALID pairs={len(audit)} complete={len(audit)}")


if __name__ == "__main__":
    main()
