#!/usr/bin/env python3
"""Fail closed on incomplete, unsupported, or overclaimed novelty decisions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALLOWED_DECISIONS = {"STRONG", "MODERATE", "WEAK", "PREVIOUSLY_KNOWN"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rg", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_rg_primary.tsv"),
    )
    parser.add_argument(
        "--audit", type=Path,
        default=Path("discovery_extension/results/novelty/pair_level_novelty_audit.tsv"),
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
        if row["novelty_decision"] not in ALLOWED_DECISIONS:
            raise SystemExit(f"ERROR: invalid novelty decision for {pair}")
        for field in ("search_databases", "search_queries", "search_date", "evidence_PMIDs_DOIs_URLs", "decision_rationale", "reviewer"):
            if not row[field].strip() or row[field].startswith("PENDING"):
                raise SystemExit(f"ERROR: completed audit lacks {field}: {pair}")
        if row["novelty_decision"] == "STRONG":
            if row["direct_prior_same_pair"] != "NO":
                raise SystemExit(f"ERROR: STRONG novelty has direct prior same-pair evidence: {pair}")
            if row["independent_replication_status"] != "INDEPENDENT_REPLICATION_PASS":
                raise SystemExit(f"ERROR: STRONG novelty lacks successful independent replication: {pair}")
    print(f"PAIR_NOVELTY_AUDIT_VALID pairs={len(audit)} complete={len(audit)}")


if __name__ == "__main__":
    main()
