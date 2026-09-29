#!/usr/bin/env python3
"""Read-only audit of available live PLACO and LAVA production checkpoints."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "extensions/brain6"
sys.path.insert(0, str(EXTENSION))
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from brain6.artifacts import verify_artifact  # noqa: E402
from materialize_lava_family_details import resolve_run_variant  # noqa: E402
from run_lava_family import verify_unit  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_placo_code_lock(repo: Path, lock: dict) -> dict:
    """Verify every implementation hash frozen into the v3 family lock."""
    checked = []
    entries = [*lock.get("extension_code_inventory", [])]
    entries.extend([
        {"path": lock["adapter_path"], "sha256": lock["adapter_sha256"]},
        {"path": lock["bundled_adapter_path"], "sha256": lock["bundled_adapter_sha256"]},
        {"path": lock["placo_source_path"], "sha256": lock["placo_source_sha256"]},
    ])
    for entry in entries:
        path = Path(entry["path"])
        if not path.is_absolute():
            path = repo / path
            try:
                path.resolve().relative_to(repo.resolve())
            except ValueError as exc:
                raise ValueError(f"PLACO lock path escapes repository: {entry['path']}") from exc
        if not path.is_file():
            raise ValueError(f"PLACO locked implementation file is missing: {path}")
        observed = sha256_file(path)
        if observed != entry["sha256"]:
            raise ValueError(
                f"PLACO locked implementation hash mismatch: {entry['path']} "
                f"expected={entry['sha256']} observed={observed}"
            )
        checked.append(str(entry["path"]))
    return {"status": "PASS", "files_checked": len(checked)}


def audit_placo_pair(root: Path, pair_id: str, expected_chunks: int) -> dict:
    chunk_ids = []
    for path in root.glob("chunk_[0-9]*"):
        receipt = path / "receipt.json"
        if not receipt.is_file():
            continue
        record = verify_artifact(path)
        if record.get("stage") != "native_job" or record.get("scientific_status") != "PASS":
            raise ValueError(f"PLACO chunk failed its receipt status: {path}")
        chunk_id = int(path.name.removeprefix("chunk_"))
        chunk_ids.append(chunk_id)
    chunk_ids.sort()
    if chunk_ids != list(range(len(chunk_ids))):
        raise ValueError(f"PLACO {pair_id} receipts are not a contiguous prefix: {chunk_ids[:5]}")
    return {
        "verified_chunks": len(chunk_ids),
        "expected_chunks": expected_chunks,
        "complete": len(chunk_ids) == expected_chunks,
        "first_missing_chunk": len(chunk_ids) if len(chunk_ids) < expected_chunks else None,
    }


def audit_placo_collation(pair_root: Path) -> dict:
    """Report pair-collation state without opening or touching live databases."""
    final = pair_root / "collated"
    lock = pair_root / "collated.lock"
    partials = sorted(p.name for p in pair_root.glob(".collated.*.partial") if p.is_dir())
    if final.is_dir():
        receipt_path = final / "receipt.json"
        if receipt_path.is_file():
            record = verify_artifact(final)
            if record.get("stage") != "collate_placo":
                raise ValueError(f"Unexpected PLACO collation receipt stage: {final}")
            status_path = final / "status.json"
            if not status_path.is_file() or json.loads(status_path.read_text(encoding="utf-8")).get("status") != "PASS":
                raise ValueError(f"PLACO collation receipt lacks passing pair QC status: {final}")
            return {"status": "COMPLETE_QC_PASS", "receipt_verified": True,
                    "lock_present": lock.exists(), "partial_directories": partials}
        return {"status": "FINAL_DIRECTORY_UNVERIFIED", "receipt_verified": False,
                "lock_present": lock.exists(), "partial_directories": partials}
    if lock.exists() and partials:
        status = "COLLATION_IN_PROGRESS"
    elif lock.exists():
        status = "LOCK_PRESENT_NO_PARTIAL_DIRECTORY"
    elif partials:
        status = "PARTIAL_DIRECTORY_NO_LOCK"
    else:
        status = "NOT_STARTED"
    return {"status": status, "receipt_verified": False,
            "lock_present": lock.exists(), "partial_directories": partials}


def audit_lava(run_root: Path, run_id: str, loci_path: Path, expected_pairs: set[str],
               maximum_failure_rate: float = 0.01,
               expected_slot_count: int | None = None,
               max_workers: int = 8) -> dict:
    with loci_path.open(newline="", encoding="utf-8") as stream:
        loci = [row["LOC"] for row in csv.DictReader(stream, delimiter=" ", skipinitialspace=True)]
    status_counts: dict[str, int] = {}
    reason_counts: dict[str, int] = {}
    no_overlap_by_pair: dict[str, int] = {}
    no_overlap_loci: set[str] = set()
    failures: list[dict[str, str]] = []
    completed = 0
    if max_workers < 1:
        raise ValueError("LAVA receipt audit max_workers must be positive")

    def verify_locus(locus_id: str) -> tuple[str, dict | None]:
        unit = run_root / "loci" / str(locus_id)
        return locus_id, verify_unit(unit, str(locus_id), expected_pairs, run_id)

    # Receipt verification is I/O-bound, especially on the external exFAT
    # results volume. A small thread pool preserves deterministic row order
    # while overlapping independent immutable file/hash checks.
    worker_count = max(1, min(max_workers, len(loci)))
    with ThreadPoolExecutor(max_workers=worker_count) as pool:
        verified_units = pool.map(verify_locus, loci)
        for locus_id, record in verified_units:
            if record is None:
                continue
            completed += 1
            with (run_root / "loci" / str(locus_id) / "worker_output/pair_results.tsv").open(
                    newline="", encoding="utf-8") as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
                    reason = row.get("reason", "")
                    if reason:
                        key = f"{row['status']}:{reason}"
                        reason_counts[key] = reason_counts.get(key, 0) + 1
                    if row["status"] == "NO_OVERLAP":
                        no_overlap_loci.add(str(locus_id))
                        pair_id = row["pair_id"]
                        no_overlap_by_pair[pair_id] = no_overlap_by_pair.get(pair_id, 0) + 1
                    if row["status"] == "FAILED":
                        failures.append({"locus_id": str(locus_id), "pair_id": row["pair_id"],
                                         "reason": reason})
    planned_slots = len(loci) * len(expected_pairs)
    if expected_slot_count is not None and planned_slots != expected_slot_count:
        raise ValueError(
            f"LAVA planned slot denominator drifted: loci × pairs={planned_slots}, "
            f"frozen family={expected_slot_count}"
        )
    if not 0 <= maximum_failure_rate < 1:
        raise ValueError("LAVA maximum failure rate must be in [0, 1)")
    failed_fraction = len(failures) / planned_slots if planned_slots else 0.0
    complete = completed == len(loci)
    ceiling_status = (
        "PASS" if complete and failed_fraction <= maximum_failure_rate else
        "FAIL" if complete else
        "PROVISIONALLY_WITHIN_CEILING" if failed_fraction <= maximum_failure_rate else
        "PROVISIONALLY_OVER_CEILING"
    )
    return {
        "verified_loci": completed,
        "planned_loci": len(loci),
        "complete": complete,
        "pair_locus_status_counts": status_counts,
        "pair_locus_reason_counts": reason_counts,
        "failed_slots": failures,
        "observed_failed_slot_fraction_of_planned_family": failed_fraction,
        "predeclared_failure_rate_ceiling": maximum_failure_rate,
        "failure_rate_ceiling_status": ceiling_status,
        "no_overlap_unique_loci": len(no_overlap_loci),
        "no_overlap_slots_by_pair": no_overlap_by_pair,
        "receipt_errors": 0,
    }


def audit_lava_manifest(repo: Path, manifest_name: str, expected_variant: str,
                        family: dict) -> dict:
    """Audit one LAVA run only after binding its root to a supported code identity."""
    manifest = json.loads(
        (repo / "brain6/manifests" / manifest_name).read_text(encoding="utf-8")
    )
    run_id = str(manifest["run_id"])
    run_root = Path(manifest["output_root"])
    if run_root.name != run_id:
        raise ValueError(
            f"LAVA output root does not match manifest run ID: {run_root} != {run_id}"
        )
    identity_path = run_root / "run_identity.json"
    if not identity_path.is_file():
        raise ValueError(f"LAVA run identity is missing: {identity_path}")
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    variant = resolve_run_variant(identity)
    if variant["name"] != expected_variant:
        raise ValueError(
            f"LAVA manifest {manifest_name} expected {expected_variant}, "
            f"but its checksum-bound identity resolves to {variant['name']}"
        )
    result = audit_lava(
        run_root, run_id, repo / family["locus_definition"]["path"],
        set(family["pairs"]),
        float(family["execution"]["maximum_locus_failure_fraction"]),
        int(family["bivariate_family"]["n_slots"]),
    )
    result["run_id"] = run_id
    result["run_variant"] = variant["name"]
    if "completed_loci" in manifest:
        result["manifest_completed_loci"] = int(manifest["completed_loci"])
        result["manifest_matches_receipts"] = (
            result["manifest_completed_loci"] == result["verified_loci"]
        )
    return result


def audit(repo: Path) -> dict:
    placo_run = json.loads((repo / "brain6/manifests/placo_family_v3_run.json").read_text(encoding="utf-8"))
    placo_lock = json.loads((repo / "brain6/config/placo_family_v3/family_lock.json").read_text(encoding="utf-8"))
    placo_code_lock = audit_placo_code_lock(repo, placo_lock)
    placo_pairs = {}
    for pair_id, spec in placo_lock["pairs"].items():
        pair_root = Path(placo_run["output_storage"]) / pair_id
        pair = audit_placo_pair(pair_root,
                                pair_id, int(spec["n_chunks"]))
        pair["collation"] = audit_placo_collation(pair_root)
        pair["manifest_completed_chunks"] = int(placo_run["pairs"][pair_id]["completed_chunks"])
        pair["manifest_matches_receipts"] = pair["manifest_completed_chunks"] == pair["verified_chunks"]
        placo_pairs[pair_id] = pair

    lava_family = json.loads((repo / "brain6/config/lava_family_v2.json").read_text(encoding="utf-8"))
    lava = audit_lava_manifest(
        repo, "lava_family_v2_run.json", "baseline", lava_family
    )
    lava_roundoff = audit_lava_manifest(
        repo, "lava_roundoff_v1_run.json", "roundoff_v1", lava_family
    )
    return {
        "audited_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "audit_status": "PASS_RECEIPT_INTEGRITY",
        "placo": {"code_lock": placo_code_lock,
                  "verified_chunks": sum(x["verified_chunks"] for x in placo_pairs.values()),
                  "planned_chunks": sum(x["expected_chunks"] for x in placo_pairs.values()),
                  "manifest_counts_match": all(x["manifest_matches_receipts"] for x in placo_pairs.values()),
                  "pairs": placo_pairs},
        "lava": lava,
        "lava_roundoff": lava_roundoff,
        "interpretation": "Receipt/collation-state audit only; baseline and corrected LAVA runs remain separate, unfinished families are not final results, and no statistical correction is performed.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument("--out", type=Path, help="Optional JSON report path; omit to print to stdout")
    args = parser.parse_args()
    result = audit(args.repo.resolve())
    rendered = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
