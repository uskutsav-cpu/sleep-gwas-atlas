#!/usr/bin/env python3
"""Collate the complete locked TWAS family and apply family-wise BH FDR."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "twas_id", "trait_id", "model_family", "modality", "context", "gene_id", "gene_symbol",
    "zscore", "uncalibrated_zscore", "effect_size", "p_value", "uncalibrated_p_value", "fdr",
    "n_snps_used", "n_snps_in_model", "gwas_N", "gwas_h2", "gwas_h2_scale",
    "coverage_fraction", "status", "model_id", "provenance_id",
]
COVERAGE_FIELDS = [
    "run_id", "trait_id", "model_id", "model_family", "context", "analysis_status", "reason",
    "result_rows", "provenance_path", "provenance_sha256",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def bh(rows: list[dict[str, object]]) -> None:
    ordered = sorted(range(len(rows)), key=lambda index: (float(rows[index]["p_value"]), str(rows[index]["twas_id"])))
    adjusted = [1.0] * len(rows); running = 1.0; total = len(rows)
    for rank_index in range(total - 1, -1, -1):
        index = ordered[rank_index]; rank = rank_index + 1
        running = min(running, float(rows[index]["p_value"]) * total / rank)
        adjusted[index] = min(1.0, running)
    for row, value in zip(rows, adjusted):
        row["fdr"] = format(value, ".15g")
        row["status"] = "SUPPORTED" if value <= 0.05 else "NO_EVIDENCE_FOUND"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/twas_run_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/twas_run_manifest.lock.json")
    parser.add_argument("--eligibility", default="results/tables/twas_trait_eligibility.tsv")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--out", default="results/tables/twas.tsv")
    parser.add_argument("--coverage-out", default="results/tables/twas_coverage.tsv")
    parser.add_argument("--provenance-out", default="results/tables/twas.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest_path, lock_path, policy_path = root / args.manifest, root / args.manifest_lock, root / args.policy
    eligibility_path = root / args.eligibility
    manifest = read_tsv(manifest_path)
    eligibility = read_tsv(eligibility_path)
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if lock.get("run_manifest_sha256") != sha256(manifest_path) or lock.get("policy_sha256") != sha256(policy_path):
        fail("TWAS manifest differs from lock")
    if lock.get("run_ids_in_locked_order") != [row["run_id"] for row in manifest]:
        fail("TWAS run family/order differs from lock")
    if lock.get("eligibility_sha256") != sha256(eligibility_path) or len(eligibility) != 45:
        fail("TWAS trait eligibility family differs from lock")
    output: list[dict[str, object]] = []
    coverage: list[dict[str, object]] = []
    for run in manifest:
        result_path = root / run["output_path"]
        provenance_path = result_path.parent / "provenance.json"
        if not result_path.is_file() or not provenance_path.is_file():
            fail(f"TWAS run is incomplete: {run['run_id']}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if provenance.get("run_id") != run["run_id"] or provenance.get("output_sha256") != sha256(result_path) or provenance.get("manifest_sha256") != sha256(manifest_path):
            fail(f"TWAS run differs from provenance: {run['run_id']}")
        calibration = provenance.get("variance_control", {})
        if (
            calibration.get("status") != "APPLIED" or calibration.get("phi_required") is not True
            or int(calibration.get("gwas_N", 0)) != int(run["gwas_N"])
            or float(calibration.get("gwas_h2", float("nan"))) != float(run["gwas_h2"])
        ):
            fail(f"TWAS variance control differs from lock: {run['run_id']}")
        with result_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != provenance.get("row_count"):
            fail(f"TWAS row count differs from provenance: {run['run_id']}")
        provenance_id = "TWAS:" + sha256(provenance_path)
        for row in rows:
            pvalue, raw_p = float(row["pvalue"]), float(row["uncalibrated_pvalue"])
            phi = float(row["phi"]); used = int(float(row["n_snps_used"])); in_model = int(float(row["n_snps_in_model"]))
            if not all(math.isfinite(value) for value in (pvalue, raw_p, phi)) or not 0 <= pvalue <= 1 or not 0 <= raw_p <= 1 or phi < 0 or used < 1 or in_model < used:
                fail(f"invalid TWAS result: {run['run_id']}/{row['gene']}")
            gene_id = row["gene"].split(".")[0]
            output.append({
                "twas_id": f"{run['run_id']}__{gene_id}", "trait_id": run["trait_id"],
                "model_family": run["model_family"], "modality": run["modality"], "context": run["context"],
                "gene_id": gene_id, "gene_symbol": row.get("gene_name", "NA") or "NA",
                "zscore": row["zscore"], "uncalibrated_zscore": row["uncalibrated_zscore"],
                "effect_size": row.get("effect_size", "NA") or "NA",
                "p_value": row["pvalue"], "uncalibrated_p_value": row["uncalibrated_pvalue"],
                "fdr": "NA", "n_snps_used": used, "n_snps_in_model": in_model,
                "gwas_N": run["gwas_N"], "gwas_h2": run["gwas_h2"], "gwas_h2_scale": run["gwas_h2_scale"],
                "coverage_fraction": format(used / in_model, ".15g"), "status": "PENDING_FDR",
                "model_id": run["model_id"], "provenance_id": provenance_id,
            })
        coverage.append({
            "run_id": run["run_id"], "trait_id": run["trait_id"], "model_id": run["model_id"],
            "model_family": run["model_family"], "context": run["context"], "result_rows": len(rows),
            "analysis_status": "COMPLETED", "reason": "PASS",
            "provenance_path": str(provenance_path.relative_to(root)), "provenance_sha256": sha256(provenance_path),
        })
    for row in eligibility:
        if row["analysis_status"] == "ELIGIBLE":
            continue
        if row["analysis_status"] != "NOT_APPLICABLE" or row["reason"] == "PASS":
            fail(f"invalid TWAS terminal eligibility outcome: {row['trait_id']}")
        coverage.append({
            "run_id": "NA", "trait_id": row["trait_id"], "model_id": "NA", "model_family": "ALL",
            "context": "ALL", "analysis_status": "NOT_APPLICABLE", "reason": row["reason"],
            "result_rows": 0, "provenance_path": str(eligibility_path.relative_to(root)),
            "provenance_sha256": sha256(eligibility_path),
        })
    by_family: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in output:
        by_family[(str(row["trait_id"]), str(row["model_family"]))].append(row)
    for rows in by_family.values(): bh(rows)
    out_text, coverage_text = table_text(FIELDS, output), table_text(COVERAGE_FIELDS, coverage)
    out, coverage_out, provenance_out = root / args.out, root / args.coverage_out, root / args.provenance_out
    provenance = {
        "analysis_id": policy["analysis_id"], "run_count": len(manifest), "result_row_count": len(output),
        "manifest_sha256": sha256(manifest_path), "manifest_lock_sha256": sha256(lock_path),
        "eligibility_sha256": sha256(eligibility_path),
        "policy_sha256": sha256(policy_path), "outputs": {
            args.out: hashlib.sha256(out_text.encode()).hexdigest(), args.coverage_out: hashlib.sha256(coverage_text.encode()).hexdigest(),
        }, "multiple_testing": policy["twas"]["multiple_testing"], "claim_limit": policy["claim_limit"],
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    if args.validate_only:
        if not out.is_file() or out.read_text(encoding="utf-8") != out_text or not coverage_out.is_file() or coverage_out.read_text(encoding="utf-8") != coverage_text or not provenance_out.is_file() or provenance_out.read_text(encoding="utf-8") != provenance_text:
            fail("TWAS aggregate differs from deterministic recomputation")
    else:
        for path, text in ((out, out_text), (coverage_out, coverage_text), (provenance_out, provenance_text)):
            if path.exists(): fail(f"TWAS aggregate already exists: {path}")
            path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text, encoding="utf-8")
    if not args.quiet: print(f"TWAS_COLLATION_OK runs={len(manifest)} rows={len(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
