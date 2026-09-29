#!/usr/bin/env python3
"""Compare receipt-verified LAVA baseline and roundoff pair-by-locus outputs.

This is a descriptive comparison only. It does not aggregate shared
univariate tests, calculate family-level correction, or promote either run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys_path = ROOT / "extensions/brain6"
sys.path.insert(0, str(sys_path))
from brain6.io import sha256


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(value: str) -> float | None:
    if value in {"", "NA", "NaN", "nan"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite numeric value: {value}")
    return result


def stable_file_inventory_hash(records: list[dict[str, str]]) -> str:
    payload = "\n".join(f"{r['path']}\t{r['sha256']}" for r in records).encode()
    return hashlib.sha256(payload).hexdigest()


def validate_audit(audit_path: Path, baseline_id: str, roundoff_id: str) -> dict[str, Any]:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("audit_status") != "PASS_RECEIPT_INTEGRITY":
        raise ValueError("checkpoint audit did not pass receipt integrity")
    for key, run_id in (("lava", baseline_id), ("lava_roundoff", roundoff_id)):
        run = audit.get(key, {})
        if (run.get("run_id") != run_id or run.get("receipt_errors") != 0 or
                run.get("complete") is not True or run.get("verified_loci") != run.get("planned_loci")):
            raise ValueError(f"{key} is not a complete receipt-verified run for {run_id}")
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()

    baseline_manifest_path = ROOT / "brain6/manifests/lava_family_v2_run.json"
    roundoff_manifest_path = ROOT / "brain6/manifests/lava_roundoff_v1_run.json"
    baseline_manifest = json.loads(baseline_manifest_path.read_text(encoding="utf-8"))
    roundoff_manifest = json.loads(roundoff_manifest_path.read_text(encoding="utf-8"))
    baseline_id = baseline_manifest["run_id"]
    roundoff_id = roundoff_manifest["run_id"]
    audit = validate_audit(args.audit, baseline_id, roundoff_id)
    baseline_root = Path(baseline_manifest["output_root"])
    roundoff_root = Path(roundoff_manifest["output_root"])
    family = json.loads((ROOT / "brain6/config/lava_family_v2.json").read_text(encoding="utf-8"))
    family_path = ROOT / "brain6/config/lava_family_v2.json"
    loci_path = ROOT / family["locus_definition"]["path"]
    loci = [line.split()[0] for line in loci_path.read_text().splitlines()[1:] if line.strip()]
    expected_pairs = tuple(family["pairs"])
    if len(loci) != 2495 or len(set(loci)) != 2495 or len(expected_pairs) != 5:
        raise ValueError("frozen LAVA locus/pair inventory does not match the intended five-pair family")

    rows: list[dict[str, str]] = []
    inventories: dict[str, list[dict[str, str]]] = {"baseline": [], "roundoff": []}
    for locus_id in loci:
        run_rows: dict[str, dict[str, dict[str, str]]] = {}
        for label, run_root in (("baseline", baseline_root), ("roundoff", roundoff_root)):
            path = run_root / "loci" / locus_id / "worker_output" / "pair_results.tsv"
            if not path.is_file():
                raise ValueError(f"missing pair results for receipt-verified locus {locus_id}: {path}")
            file_hash = sha256(path)
            inventories[label].append({"path": str(path), "sha256": file_hash})
            local_rows = read_tsv(path)
            by_pair = {row["pair_id"]: row for row in local_rows}
            if len(local_rows) != len(expected_pairs) or set(by_pair) != set(expected_pairs):
                raise ValueError(f"pair result inventory differs at {label} locus {locus_id}")
            run_rows[label] = by_pair
        for pair_id in expected_pairs:
            before, after = run_rows["baseline"][pair_id], run_rows["roundoff"][pair_id]
            if before["locus_id"] != locus_id or after["locus_id"] != locus_id:
                raise ValueError(f"locus key mismatch for {pair_id}/{locus_id}")
            bp, ap = number(before["p"]), number(after["p"])
            br, ar = number(before["local_rg"]), number(after["local_rg"])
            rows.append({
                "pair_id": pair_id, "locus_id": locus_id,
                "baseline_status": before["status"], "roundoff_status": after["status"],
                "status_transition": f"{before['status']}->{after['status']}",
                "baseline_p": "NA" if bp is None else repr(bp),
                "roundoff_p": "NA" if ap is None else repr(ap),
                "delta_p_roundoff_minus_baseline": "NA" if bp is None or ap is None else repr(ap - bp),
                "baseline_local_rg": "NA" if br is None else repr(br),
                "roundoff_local_rg": "NA" if ar is None else repr(ar),
                "delta_local_rg_roundoff_minus_baseline": "NA" if br is None or ar is None else repr(ar - br),
                "baseline_reason": before["reason"], "roundoff_reason": after["reason"],
            })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    transitions = Counter(row["status_transition"] for row in rows)
    both_tested = [row for row in rows if row["baseline_status"] == row["roundoff_status"] == "TESTED"]
    numeric_pairs = [row for row in both_tested if row["delta_p_roundoff_minus_baseline"] != "NA"
                     and row["delta_local_rg_roundoff_minus_baseline"] != "NA"]
    summary_rows: list[dict[str, str]] = []
    all_statuses = ("TESTED", "UNIVARIATE_UNDERPOWERED", "NO_OVERLAP", "FAILED")
    for pair_id in expected_pairs:
        pair_rows = [row for row in rows if row["pair_id"] == pair_id]
        pair_record = {"pair_id": pair_id, "loci": str(len(pair_rows))}
        for label in ("baseline", "roundoff"):
            for status in all_statuses:
                pair_record[f"{label}_{status.lower()}_loci"] = str(sum(row[f"{label}_status"] == status for row in pair_rows))
        pair_record["status_changed_loci"] = str(sum(row["baseline_status"] != row["roundoff_status"] for row in pair_rows))
        pair_record["both_tested_loci"] = str(sum(row["baseline_status"] == row["roundoff_status"] == "TESTED" for row in pair_rows))
        summary_rows.append(pair_record)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    with args.summary_output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    retry_provenance_path = (roundoff_root / "targeted_retries/envfix-20260924-001/retry_provenance.json")
    retry_provenance = json.loads(retry_provenance_path.read_text(encoding="utf-8"))
    retry_receipt_path = roundoff_root / "targeted_retries/envfix-20260924-001/locus_740/receipt.json"
    original_receipt_path = roundoff_root / "loci/740/receipt.json"
    if (retry_provenance.get("status") != "PASS_RETRY_ARTIFACT" or
            retry_provenance.get("retry_pair_result", {}).get("status") != "UNIVARIATE_UNDERPOWERED" or
            retry_provenance.get("original_receipt_sha256") != sha256(original_receipt_path) or
            retry_provenance.get("retry_receipt_sha256") != sha256(retry_receipt_path)):
        raise ValueError("isolated environmental retry provenance/receipts failed validation")

    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6-lava-baseline-roundoff-descriptive-comparison-v1",
        "status": "DESCRIPTIVE_COMPLETE_FAMILY_QC_UNRESOLVED",
        "interpretation": "Receipt-verified pair-by-locus comparison only. It does not resolve the roundoff aggregator's shared-univariate disagreement, apply family correction, or promote local-rg inference.",
        "baseline_run_id": baseline_id,
        "roundoff_run_id": roundoff_id,
        "audit_path": str(args.audit),
        "audit_sha256": sha256(args.audit),
        "audit_time_utc": audit["audited_at_utc"],
        "planned_pair_locus_slots": len(rows),
        "pair_ids": list(expected_pairs),
        "loci": len(loci),
        "status_transitions": dict(sorted(transitions.items())),
        "both_tested_slots": len(both_tested),
        "both_tested_with_numeric_p_and_rg": len(numeric_pairs),
        "numeric_difference_summary": {
            "max_abs_delta_p": max((abs(float(r["delta_p_roundoff_minus_baseline"])) for r in numeric_pairs), default=None),
            "max_abs_delta_local_rg": max((abs(float(r["delta_local_rg_roundoff_minus_baseline"])) for r in numeric_pairs), default=None),
            "significance_status_changes_at_p_0_05": sum(
                (float(r["baseline_p"]) < 0.05) != (float(r["roundoff_p"]) < 0.05) for r in numeric_pairs
            ),
        },
        "source_manifests": {
            str(baseline_manifest_path): sha256(baseline_manifest_path),
            str(roundoff_manifest_path): sha256(roundoff_manifest_path),
            str(family_path): sha256(family_path),
            str(loci_path): sha256(loci_path),
        },
        "pair_result_file_inventory_sha256": {
            label: stable_file_inventory_hash(records) for label, records in inventories.items()
        },
        "comparison_table_path": str(args.output),
        "comparison_table_sha256": sha256(args.output),
        "summary_table_path": str(args.summary_output),
        "summary_table_sha256": sha256(args.summary_output),
        "builder_path": str(Path(__file__)),
        "builder_sha256": sha256(Path(__file__)),
        "roundoff_aggregation_error": "Pair-specific inputs disagree on the shared univariate test: insomnia/2",
        "environmental_retry": {
            "retry_id": "envfix-20260924-001",
            "status": "PASS_RETRY_ARTIFACT",
            "artifact_path": str(retry_provenance_path),
            "artifact_sha256": sha256(retry_provenance_path),
            "receipt_sha256": sha256(retry_receipt_path),
            "original_receipt_sha256_unchanged": sha256(original_receipt_path),
            "pair_result": retry_provenance["retry_pair_result"],
            "interpretation": "Retry changed the preserved environmental failure to UNIVARIATE_UNDERPOWERED; it remains isolated and is not merged into the family.",
        },
    }
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: provenance[k] for k in ("status", "planned_pair_locus_slots", "status_transitions", "both_tested_slots", "numeric_difference_summary")}, sort_keys=True))


if __name__ == "__main__":
    main()
