#!/usr/bin/env python3
"""Audit pair-context univariate outputs against canonical trait-only LAVA v3.

Receipt-verified descriptive audit only; it never changes or aggregates runs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
from brain6.io import sha256

_runner_path = ROOT / "brain6/scripts/run_lava_family.py"
_runner_spec = importlib.util.spec_from_file_location("lava_pair_audit_runner", _runner_path)
assert _runner_spec is not None and _runner_spec.loader is not None
_runner = importlib.util.module_from_spec(_runner_spec)
_runner_spec.loader.exec_module(_runner)
EXPECTED_PAIRS = _runner.EXPECTED_PAIRS
verify_unit = _runner.verify_unit


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def numeric(raw: str) -> float | None:
    if raw in {"", "NA", "NaN", "nan"}:
        return None
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError(f"non-finite value: {raw}")
    return value


def stable_inventory_hash(rows: list[dict[str, str]]) -> str:
    payload = "\n".join(f"{row['path']}\t{row['sha256']}" for row in rows).encode()
    return hashlib.sha256(payload).hexdigest()


def summarize_context(contexts: list[dict[str, str]], canonical: dict[str, str],
                      tolerance: float = 1e-12) -> dict[str, Any]:
    ps = [numeric(row["p"]) for row in contexts]
    valid = [value for value in ps if value is not None]
    spread = max(valid) - min(valid) if len(valid) > 1 else 0.0
    canonical_p = numeric(canonical["p"])
    deltas = [abs(value - canonical_p) for value in valid] if canonical_p is not None else []
    return {
        "n_contexts": len(contexts),
        "within_context_p_range": spread,
        "contexts_disagree": spread > tolerance,
        "canonical_p": canonical_p,
        "max_abs_context_minus_canonical": max(deltas, default=0.0),
        "differs_from_canonical": bool(deltas and max(deltas) > tolerance),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()

    audit = json.loads(args.audit.read_text())
    if audit.get("audit_status") != "PASS_RECEIPT_INTEGRITY":
        raise ValueError("latest production audit does not pass receipt integrity")
    manifests = {
        "baseline": ROOT / "brain6/manifests/lava_family_v2_run.json",
        "roundoff": ROOT / "brain6/manifests/lava_roundoff_v1_run.json",
    }
    ids: dict[str, str] = {}
    roots: dict[str, Path] = {}
    for label, manifest_path in manifests.items():
        manifest = json.loads(manifest_path.read_text())
        ids[label] = manifest["run_id"]
        roots[label] = Path(manifest["output_root"])
        audit_key = "lava" if label == "baseline" else "lava_roundoff"
        check = audit[audit_key]
        if (check.get("run_id") != ids[label] or check.get("receipt_errors") != 0
                or check.get("complete") is not True or check.get("verified_loci") != 2495):
            raise ValueError(f"{label} receipt audit does not match its run manifest")

    canonical = read_tsv(args.canonical)
    canon = {(r["phen"], r["locus_id"]): r for r in canonical}
    if len(canon) != 17465 or len(canonical) != 17465:
        raise ValueError("canonical aggregate must contain exactly 17,465 unique cells")

    contexts: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    file_inventory: dict[str, list[dict[str, str]]] = {k: [] for k in manifests}
    verified_receipts: dict[str, int] = defaultdict(int)
    for label, root in roots.items():
        for locus in range(1, 2496):
            locus_id = str(locus)
            unit = root / "loci" / locus_id
            receipt = verify_unit(unit, locus_id, set(EXPECTED_PAIRS), ids[label])
            if receipt is None:
                raise ValueError(f"missing {label} receipt for locus {locus_id}")
            verified_receipts[label] += 1
            path = unit / "worker_output/univariate.tsv"
            file_inventory[label].append({"path": str(path), "sha256": sha256(path)})
            for row in read_tsv(path):
                trait, row_locus, pair = row["phen"].strip('"'), row["locus_id"].strip('"'), row["pair_id"].strip('"')
                if row_locus != locus_id or pair not in EXPECTED_PAIRS:
                    raise ValueError(f"unexpected univariate row identity in {label}/{locus_id}")
                key = (trait, row_locus)
                if key not in canon:
                    raise ValueError(f"pair-context cell is absent from canonical family: {key}")
                contexts[(label, trait, row_locus)].append({
                    "pair_id": pair,
                    "p": row["p"].strip('"'),
                    "h2_obs": row["h2.obs"].strip('"'),
                    "h2_latent": row["h2.latent"].strip('"'),
                })
    if any(n != 2495 for n in verified_receipts.values()):
        raise ValueError("not all prior-run receipts were independently verified")

    output_rows: list[dict[str, str]] = []
    summary: list[dict[str, str]] = []
    tol = 1e-12
    for label in manifests:
        for trait in sorted({k[0] for k in canon}):
            keys = [(trait, str(locus)) for locus in range(1, 2496)]
            observed = 0
            repeated = 0
            conflicting = 0
            canonical_tested = 0
            context_vs_canonical_conflicts = 0
            max_within = 0.0
            max_vs_canonical = 0.0
            for key in keys:
                base = canon[key]
                ctx = contexts.get((label, *key), [])
                if not ctx:
                    continue
                observed += 1
                comparison = summarize_context(ctx, base, tol)
                spread = comparison["within_context_p_range"]
                is_repeated = comparison["n_contexts"] > 1
                is_conflict = comparison["contexts_disagree"]
                if is_repeated:
                    repeated += 1
                if is_conflict:
                    conflicting += 1
                    max_within = max(max_within, spread)
                canonical_p = comparison["canonical_p"]
                vs_conflict = comparison["differs_from_canonical"]
                if base["status"] == "TESTED":
                    canonical_tested += 1
                    if vs_conflict:
                        context_vs_canonical_conflicts += 1
                        max_vs_canonical = max(max_vs_canonical, comparison["max_abs_context_minus_canonical"])
                for context in ctx:
                    p = numeric(context["p"])
                    delta = None if canonical_p is None or p is None else p - canonical_p
                    output_rows.append({
                        "run": label, "run_id": ids[label], "trait": trait, "locus_id": key[1],
                        "pair_id": context["pair_id"], "pair_context_p": "NA" if p is None else repr(p),
                        "pair_context_h2_obs": context["h2_obs"], "pair_context_h2_latent": context["h2_latent"],
                        "n_pair_contexts": str(len(ctx)), "within_context_p_range": repr(spread),
                        "contexts_disagree_gt_1e-12": str(is_conflict).upper(),
                        "canonical_v3_status": base["status"],
                        "canonical_v3_p": "NA" if canonical_p is None else repr(canonical_p),
                        "delta_context_minus_canonical": "NA" if delta is None else repr(delta),
                        "canonical_comparison_admissible": str(canonical_p is not None).upper(),
                    })
            summary.append({
                "run": label, "trait": trait, "observed_unique_trait_locus": str(observed),
                "missing_unique_trait_locus": str(2495 - observed), "repeated_keys": str(repeated),
                "repeated_keys_disagree_gt_1e-12": str(conflicting),
                "canonical_tested_keys_with_pair_context": str(canonical_tested),
                "canonical_tested_keys_context_differs_gt_1e-12": str(context_vs_canonical_conflicts),
                "max_within_context_p_range": repr(max_within),
                "max_abs_context_minus_canonical_p": repr(max_vs_canonical),
            })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    for path, rows in ((args.output, output_rows), (args.summary_output, summary)):
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    provenance: dict[str, Any] = {
        "schema_version": 1,
        "analysis_id": "brain6-lava-pair-context-univariate-audit-v1",
        "status": "DESCRIPTIVE_COMPLETE_NO_FAMILY_PROMOTION",
        "interpretation": "Compares receipt-verified pair-context univariate values with the independently computed canonical trait-only v3 family. This audit does not combine runs, repair old receipts, or authorize downstream inference.",
        "run_ids": ids,
        "receipt_verified_loci": dict(verified_receipts),
        "intended_unique_trait_locus_tests": 17465,
        "pair_context_rows": len(output_rows),
        "canonical_aggregate_sha256": sha256(args.canonical),
        "canonical_decision_sha256": sha256(ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json"),
        "checkpoint_audit_sha256": sha256(args.audit),
        "source_manifests": {str(p): sha256(p) for p in manifests.values()},
        "univariate_inventory_sha256": {k: stable_inventory_hash(v) for k, v in file_inventory.items()},
        "comparison_tolerance": 1e-12,
        "output_sha256": sha256(args.output),
        "summary_sha256": sha256(args.summary_output),
        "builder": str(Path(__file__)),
        "builder_sha256": sha256(Path(__file__)),
    }
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"pair_context_rows": len(output_rows), "summary_rows": len(summary), "verified_loci": dict(verified_receipts), "output_sha256": provenance["output_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
