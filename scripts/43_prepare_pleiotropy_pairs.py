#!/usr/bin/env python3
"""Freeze the exact 396-pair PLACO+/conjFDR input-readiness family."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT_SPEC = importlib.util.spec_from_file_location(
    "pleiotropy_preflight", ROOT / "scripts/42_pleiotropy_preflight.py"
)
if PREFLIGHT_SPEC is None or PREFLIGHT_SPEC.loader is None:
    raise RuntimeError("could not load pleiotropy preflight helpers")
preflight = importlib.util.module_from_spec(PREFLIGHT_SPEC)
PREFLIGHT_SPEC.loader.exec_module(preflight)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--out", default="results/tables/pleiotropy_pair_manifest.tsv")
    parser.add_argument("--lock-out", default="results/tables/pleiotropy_pair_manifest.lock.json")
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads(
        (root / "config/pleiotropy_analysis_policy.json").read_text(encoding="utf-8")
    )
    panel = preflight.mixer.read_tsv(root / "config/analysis_panel.tsv")
    traits = [row["trait_id"] for row in panel]
    if len(traits) != 45 or len(set(traits)) != 45:
        raise SystemExit("ERROR: panel is not the exact locked 45-trait set")
    sleeps = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
    diseases = [row["trait_id"] for row in panel if row["domain"] != "sleep"]
    if len(sleeps) != 12 or len(diseases) != 33:
        raise SystemExit("ERROR: expected 12 sleep and 33 non-sleep traits")
    rg_rows = preflight.mixer.read_tsv(root / "results/tables/rg_matrix.tsv")
    rg_by_pair = {(row["sleep_trait"], row["disease_trait"]): row for row in rg_rows}
    expected = [(sleep, disease) for sleep in sleeps for disease in diseases]
    if len(rg_rows) != 396 or set(rg_by_pair) != set(expected):
        raise SystemExit("ERROR: rg table is not the exact 396-pair family")

    planned_prefilters = {
        row["trait_id"]
        for row in preflight.mixer.read_tsv(root / "config/hm3_prefilter_plans.tsv")
    }
    inputs: dict[str, dict[str, object]] = {}
    for trait in traits:
        data, qc, strategy = preflight.mixer.choose_harmonized(root, trait)
        if not strategy and trait not in planned_prefilters:
            strategy = "not supplied"
        ready = data.is_file() and qc.is_file() and strategy == "not supplied"
        inputs[trait] = {
            "path": data,
            "qc": qc,
            "strategy": strategy or "UNKNOWN",
            "ready": ready,
            "sha256": sha256(data) if data.is_file() else "",
            "bytes": data.stat().st_size if data.is_file() else 0,
        }

    rows = []
    for sleep, disease in expected:
        first, second = inputs[sleep], inputs[disease]
        ready = bool(first["ready"] and second["ready"])
        rg = rg_by_pair[(sleep, disease)]
        rows.append({
            "pair_id": f"{sleep}__{disease}",
            "sleep_trait": sleep,
            "non_sleep_trait": disease,
            "analysis_tier": rg["analysis_tier"],
            "interpretation_status": rg["interpretation_status"],
            "global_rg": rg["rg"],
            "global_rg_se": rg["se"],
            "global_rg_p": rg["p"],
            "global_rg_fdr": rg["fdr"],
            "sleep_sumstats": str(Path(first["path"]).relative_to(root)),
            "sleep_sumstats_bytes": first["bytes"],
            "sleep_sumstats_sha256": first["sha256"],
            "sleep_prefilter_strategy": first["strategy"],
            "non_sleep_sumstats": str(Path(second["path"]).relative_to(root)),
            "non_sleep_sumstats_bytes": second["bytes"],
            "non_sleep_sumstats_sha256": second["sha256"],
            "non_sleep_prefilter_strategy": second["strategy"],
            "placo_pair_threshold": policy["placo_locked_pair_family_threshold"],
            "conjfdr_threshold": policy["pleiofdr_conjfdr_threshold"],
            "input_status": "READY_FULL_SUMSTATS" if ready else "BLOCKED_FULL_SUMSTATS",
            "blocker": "" if ready else "one or both traits lack full non-HapMap3 post-QC summary statistics",
        })
    output = root / args.out
    write_tsv(output, rows)
    lock = {
        "analysis_id": policy["analysis_id"],
        "policy": "config/pleiotropy_analysis_policy.json",
        "policy_sha256": sha256(root / "config/pleiotropy_analysis_policy.json"),
        "manifest": str(output.relative_to(root)),
        "manifest_sha256": sha256(output),
        "panel_order_sha256": hashlib.sha256("".join(f"{trait}\n" for trait in traits).encode()).hexdigest(),
        "pair_ids_in_locked_order": [row["pair_id"] for row in rows],
        "pair_count": len(rows),
        "ready_pair_count": sum(row["input_status"] == "READY_FULL_SUMSTATS" for row in rows),
        "blocked_pair_count": sum(row["input_status"] != "READY_FULL_SUMSTATS" for row in rows),
        "global_rg_filters_pair_eligibility": False,
        "placo_source_sha256": policy["placo_source_sha256"],
        "pleiofdr_commit": policy["pleiofdr_commit"],
        "locus_definition_sha256": policy["locus_definition_sha256"],
    }
    lock_path = root / args.lock_out
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = lock_path.with_suffix(lock_path.suffix + ".tmp")
    temporary.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(lock_path)
    print(
        f"Pleiotropy pair family: {len(rows)} locked; "
        f"{lock['ready_pair_count']} full-input ready; {lock['blocked_pair_count']} blocked"
    )
    return 0 if lock["blocked_pair_count"] == 0 or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
