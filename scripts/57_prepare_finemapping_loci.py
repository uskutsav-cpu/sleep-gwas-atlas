#!/usr/bin/env python3
"""Freeze every primary cross-method shared locus before fine-mapping results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


FIELDS = [
    "comparison_id", "shared_locus_id", "pair_id", "sleep_trait", "non_sleep_trait",
    "analysis_tier", "lava_block_id", "chromosome", "start_bp", "end_bp",
    "placo_lead_snp", "placo_lead_p", "conjfdr_lead_snp", "conjfdr_lead_fdr", "effect_direction",
    "sleep_full_input", "non_sleep_full_input", "reference_prefix",
    "summary1_path", "summary2_path", "variant_order_path", "ld_path", "task_path",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def full_input(root: Path, trait: str) -> str:
    candidates = [
        (
            root / "data/harmonized_mixer_full" / f"{trait}.harmonized.tsv.gz",
            root / "data/harmonized_mixer_full" / f"{trait}.qc.txt",
        ),
        (
            root / "data/harmonized" / f"{trait}.harmonized.tsv.gz",
            root / "data/harmonized" / f"{trait}.qc.txt",
        ),
    ]
    for harmonized, qc in candidates:
        if harmonized.is_file() and qc.is_file() and not any(
            line.startswith("prefilter_strategy\t") for line in qc.read_text(encoding="utf-8").splitlines()
        ):
            return str(harmonized.relative_to(root))
    fail(f"full non-HapMap3 harmonized input is absent for {trait}")
    raise AssertionError


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--shared", default="results/atlas/shared_loci.tsv")
    parser.add_argument("--policy", default="config/fine_mapping_analysis_policy.json")
    parser.add_argument("--preflight", default="results/tables/fine_mapping_preflight.json")
    parser.add_argument("--out", default="results/tables/fine_mapping_locus_manifest.tsv")
    parser.add_argument("--lock", default="results/tables/fine_mapping_locus_manifest.lock.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    shared_path = root / args.shared
    policy_path = root / args.policy
    preflight_path = root / args.preflight
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    if not preflight.get("ready") or preflight.get("policy_sha256") != sha256(policy_path):
        fail("fine-mapping preflight is absent, blocked, or policy-drifted")
    fields, shared = read_tsv(shared_path)
    required = {
        "shared_locus_id", "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
        "locus_id", "CHR", "START", "STOP", "placo_lead_snp", "placo_lead_p",
        "conjfdr_lead_snp", "conjfdr_lead_fdr", "evidence_status",
        "effect_direction",
    }
    if not required.issubset(fields):
        fail("shared-locus table lacks required cross-method fields")
    primary = [row for row in shared if row["analysis_tier"] == "PRIMARY_PHASE1"]
    if not primary:
        fail("no primary cross-method shared locus is available")
    if any(row["evidence_status"] != "PLACO_PLUS_AND_CONJFDR_SAME_LOCKED_LD_BLOCK" for row in primary):
        fail("primary fine-mapping family contains a non-consensus locus")
    if len({row["shared_locus_id"] for row in primary}) != len(primary):
        fail("duplicate primary shared-locus identifier")
    rows = []
    for row in primary:
        identity = row["shared_locus_id"]
        chromosome, start, end = int(row["CHR"]), int(row["START"]), int(row["STOP"])
        if not 1 <= chromosome <= 22 or not 0 < start <= end:
            fail(f"invalid locus interval for {identity}")
        base = f"data/fine_mapping/{identity}"
        rows.append({
            "comparison_id": f"{identity}__TRAIT_TRAIT",
            "shared_locus_id": identity,
            "pair_id": row["pair_id"],
            "sleep_trait": row["sleep_trait"],
            "non_sleep_trait": row["non_sleep_trait"],
            "analysis_tier": row["analysis_tier"],
            "lava_block_id": row["locus_id"],
            "chromosome": str(chromosome),
            "start_bp": str(start),
            "end_bp": str(end),
            "placo_lead_snp": row["placo_lead_snp"],
            "placo_lead_p": row["placo_lead_p"],
            "conjfdr_lead_snp": row["conjfdr_lead_snp"],
            "conjfdr_lead_fdr": row["conjfdr_lead_fdr"],
            "effect_direction": row["effect_direction"],
            "sleep_full_input": full_input(root, row["sleep_trait"]),
            "non_sleep_full_input": full_input(root, row["non_sleep_trait"]),
            "reference_prefix": policy["reference"]["prefix"],
            "summary1_path": f"{base}/sleep.tsv.gz",
            "summary2_path": f"{base}/non_sleep.tsv.gz",
            "variant_order_path": f"{base}/variants.tsv",
            "ld_path": f"{base}/ld.tsv.gz",
            "task_path": f"results/fine_mapping/tasks/{identity}.tsv",
        })
    out = root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = out.with_suffix(out.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(out)
    lock = {
        "schema_version": "atlas-v1.0-finemapping-locus-manifest.1",
        "locked_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "results_accessed_before_lock": False,
        "selection_rule": policy["entry_rule"],
        "locus_count": len(rows),
        "comparison_ids_in_locked_order": [row["comparison_id"] for row in rows],
        "manifest_sha256": sha256(out),
        "shared_loci_sha256": sha256(shared_path),
        "policy_sha256": sha256(policy_path),
        "preflight_sha256": sha256(preflight_path),
        "claim_limit": policy["claim_limit"],
    }
    lock_path = root / args.lock
    lock_path.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FINEMAPPING_LOCUS_FAMILY_LOCKED loci={len(rows)} result_free=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
