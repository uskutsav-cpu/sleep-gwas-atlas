#!/usr/bin/env python3
"""Resume-safe streaming acquisition, harmonization, and LDSC munging per locked trait."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MIN_FREE_BYTES = 700 * 1024**2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def existing_receipt_complete(
    receipt_path: Path,
    munged: Path,
    log: Path,
    *,
    policy_sha256: str,
    code_sha256: dict[str, str],
    snapshot_sha256: str,
    reference_sha256: str,
) -> bool:
    if not receipt_path.is_file() or not munged.is_file() or not log.is_file():
        return False
    try:
        receipt = json.loads(receipt_path.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    return (
        receipt.get("pipeline_status") == "STREAM_HARMONIZE_MUNGE_PASS"
        and receipt.get("munged_output_sha256") == sha256(munged)
        and receipt.get("munging_log_sha256") == sha256(log)
        and receipt.get("source_verification", {}).get("verification_status") == "PASS"
        and receipt.get("harmonization_policy_sha256") == policy_sha256
        and receipt.get("pipeline_code_sha256") == code_sha256
        and receipt.get("remote_snapshot_sha256") == snapshot_sha256
        and receipt.get("reference_sha256") == reference_sha256
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--trait-id", action="append")
    selection.add_argument("--all", action="store_true")
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/remote_object_snapshot.tsv"),
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz"),
    )
    parser.add_argument(
        "--reference-provenance",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/hm3_variant_reference_build.json"),
    )
    parser.add_argument(
        "--snapshot-provenance",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/remote_object_snapshot.json"),
    )
    parser.add_argument(
        "--streaming-contract",
        type=Path,
        default=Path("discovery_extension/config/streaming_acquisition_contract.json"),
    )
    parser.add_argument(
        "--harmonization-policy",
        type=Path,
        default=Path("discovery_extension/config/extension_harmonization_policy.json"),
    )
    parser.add_argument(
        "--acknowledge-network-gib",
        type=float,
        required=True,
        help="Must be at least the exact compressed GiB selected for this run.",
    )
    parser.add_argument("--keep-harmonized", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    subprocess.run(
        [str(ROOT / ".venv/bin/python"), "discovery_extension/scripts/00_verify_core_checkpoint.py"],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [str(ROOT / ".venv/bin/python"), "discovery_extension/scripts/05_validate_extension_panel.py"],
        cwd=ROOT,
        check=True,
    )
    required_files = (
        args.reference,
        args.reference_provenance,
        args.snapshot,
        args.snapshot_provenance,
        args.streaming_contract,
        args.harmonization_policy,
    )
    for required in required_files:
        if not required.is_file():
            raise SystemExit(f"ERROR: streaming prerequisite is absent: {required}")
    reference_provenance = json.loads(args.reference_provenance.read_text())
    reference_receipt = reference_provenance.get("variant_manifest_stream_receipt", {})
    if (
        reference_provenance.get("output_sha256") != sha256(args.reference)
        or reference_receipt.get("verification_status") != "PASS"
        or reference_receipt.get("observed_size_bytes")
        != reference_receipt.get("expected_size_bytes")
        or reference_receipt.get("observed_md5") != reference_receipt.get("expected_md5")
    ):
        raise SystemExit("ERROR: streamed HapMap3 identity reference provenance does not validate")
    snapshot_provenance = json.loads(args.snapshot_provenance.read_text())
    if (
        snapshot_provenance.get("failure_count") != 0
        or snapshot_provenance.get("output_sha256") != sha256(args.snapshot)
    ):
        raise SystemExit("ERROR: remote object snapshot provenance does not validate")
    streaming_contract = json.loads(args.streaming_contract.read_text())
    if (
        streaming_contract.get("defined_before_extension_h2_or_rg") is not True
        or streaming_contract.get("panel_membership_changed") is not False
        or streaming_contract.get("statistical_thresholds_changed") is not False
    ):
        raise SystemExit("ERROR: streaming acquisition contract lacks the pre-result boundary")
    harmonization_policy = json.loads(args.harmonization_policy.read_text())
    if harmonization_policy.get("panel_sha256") != sha256(args.panel):
        raise SystemExit("ERROR: harmonization policy is not bound to the locked candidate panel")
    policy_sha256 = sha256(args.harmonization_policy)
    snapshot_sha256 = sha256(args.snapshot)
    reference_sha256 = sha256(args.reference)
    scientific_code_paths = [
        Path("discovery_extension/scripts/streaming_io.py"),
        Path("discovery_extension/scripts/10_harmonize_panukbb.py"),
        Path("discovery_extension/scripts/11_munge_extension.sh"),
    ]
    code_sha256 = {str(path): sha256(path) for path in scientific_code_paths}
    orchestration_code_sha256 = sha256(
        Path("discovery_extension/scripts/10_stream_and_munge_panukbb.py")
    )

    panel_rows = read_tsv(args.panel)
    panel = {row["extension_trait_id"]: row for row in panel_rows}
    selected = [row["extension_trait_id"] for row in panel_rows] if args.all else args.trait_id
    assert selected is not None
    missing = [trait_id for trait_id in selected if trait_id not in panel]
    if missing:
        raise SystemExit(f"ERROR: selected traits are not in the locked panel: {missing}")

    snapshots = read_tsv(args.snapshot)
    source_snapshot = {
        row["extension_trait_id"]: row
        for row in snapshots
        if row["object_role"] == "phenotype_sumstats"
    }
    absent_snapshot = [trait_id for trait_id in selected if trait_id not in source_snapshot]
    if absent_snapshot:
        raise SystemExit(f"ERROR: selected traits lack a remote source snapshot: {absent_snapshot}")
    invalid_snapshot = [
        trait_id
        for trait_id in selected
        if not source_snapshot[trait_id]["head_verification_status"].startswith("HEAD_SIZE_VERSION")
        or "versionId=" not in source_snapshot[trait_id]["versioned_url"]
    ]
    if invalid_snapshot:
        raise SystemExit(f"ERROR: selected traits lack a passing version-pinned HEAD: {invalid_snapshot}")
    shared_reference_rows = [
        row
        for row in snapshots
        if row["extension_trait_id"] == "SHARED_VARIANT_REFERENCE"
        and row["filename"] == "full_variant_qc_metrics.txt.bgz"
    ]
    if (
        len(shared_reference_rows) != 1
        or shared_reference_rows[0]["versioned_url"] != reference_provenance.get("variant_manifest")
    ):
        raise SystemExit("ERROR: reference provenance is not bound to the current remote snapshot")

    required_bytes = sum(int(panel[trait_id]["source_file_size_bytes"]) for trait_id in selected)
    required_gib = required_bytes / 1024**3
    if args.acknowledge_network_gib + 1e-9 < required_gib:
        raise SystemExit(
            f"ERROR: selected sources total {required_gib:.6f} GiB; "
            f"--acknowledge-network-gib was {args.acknowledge_network_gib:.6f}"
        )
    print(
        f"STREAMING_FAMILY_ACKNOWLEDGED traits={len(selected)} compressed_gib={required_gib:.6f} "
        f"keep_harmonized={str(args.keep_harmonized).lower()}"
    )
    if args.preflight_only:
        print(
            f"STREAMING_PREFLIGHT_PASS traits={len(selected)} compressed_gib={required_gib:.6f} "
            f"free_gib={shutil.disk_usage(ROOT).free / 1024**3:.3f}"
        )
        return

    for index, trait_id in enumerate(selected, start=1):
        free = shutil.disk_usage(ROOT).free
        if free < MIN_FREE_BYTES:
            raise SystemExit(
                f"ERROR: operational free space fell below {MIN_FREE_BYTES / 1024**2:.0f} MiB "
                f"before {trait_id}; completed {index - 1}/{len(selected)}"
            )
        harmonized_rel = Path(f"discovery_extension/data/harmonized/{trait_id}.txt.gz")
        qc_rel = Path(f"discovery_extension/results/qc/harmonization/{trait_id}.tsv")
        receipt_rel = Path(f"discovery_extension/provenance/streaming_receipts/{trait_id}.json")
        munged_rel = Path(f"discovery_extension/data/munged/{trait_id}.sumstats.gz")
        munging_log_rel = Path(f"discovery_extension/data/munged/{trait_id}.log")
        harmonized = ROOT / harmonized_rel
        qc = ROOT / qc_rel
        receipt_path = ROOT / receipt_rel
        munged = ROOT / munged_rel
        munging_log = ROOT / munging_log_rel
        if existing_receipt_complete(
            receipt_path,
            munged,
            munging_log,
            policy_sha256=policy_sha256,
            code_sha256=code_sha256,
            snapshot_sha256=snapshot_sha256,
            reference_sha256=reference_sha256,
        ):
            print(f"STREAMING_TRAIT_SKIP_VERIFIED index={index}/{len(selected)} trait={trait_id}")
            continue

        harmonized.parent.mkdir(parents=True, exist_ok=True)
        qc.parent.mkdir(parents=True, exist_ok=True)
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        remote = source_snapshot[trait_id]
        print(
            f"STREAMING_TRAIT_START index={index}/{len(selected)} trait={trait_id} "
            f"compressed_gib={int(panel[trait_id]['source_file_size_bytes']) / 1024**3:.6f}"
        )
        subprocess.run(
            [
                str(ROOT / ".venv/bin/python"),
                "discovery_extension/scripts/10_harmonize_panukbb.py",
                "--trait-id", trait_id,
                "--source-url", remote["versioned_url"],
                "--panel", str(args.panel),
                "--reference", str(args.reference),
                "--out", str(harmonized_rel),
                "--qc-out", str(qc_rel),
                "--receipt-out", str(receipt_rel),
            ],
            cwd=ROOT,
            check=True,
        )
        subprocess.run(
            ["bash", "discovery_extension/scripts/11_munge_extension.sh", trait_id],
            cwd=ROOT,
            check=True,
        )
        if not munged.is_file() or not munging_log.is_file():
            raise SystemExit(f"ERROR: LDSC munging did not produce both outputs for {trait_id}")
        receipt = json.loads(receipt_path.read_text())
        receipt.update(
            {
                "pipeline_status": "STREAM_HARMONIZE_MUNGE_PASS",
                "munged_output": str(munged.relative_to(ROOT)),
                "munged_output_sha256": sha256(munged),
                "munging_log": str(munging_log.relative_to(ROOT)),
                "munging_log_sha256": sha256(munging_log),
                "remote_snapshot": str(args.snapshot),
                "remote_snapshot_sha256": snapshot_sha256,
                "reference": str(args.reference),
                "reference_sha256": reference_sha256,
                "harmonization_policy": str(args.harmonization_policy),
                "harmonization_policy_sha256": policy_sha256,
                "pipeline_code_sha256": code_sha256,
                "orchestration_code_sha256": orchestration_code_sha256,
                "sealed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "harmonized_intermediate_retained": args.keep_harmonized,
            }
        )
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        if not args.keep_harmonized:
            harmonized.unlink()
            receipt["harmonized_intermediate_removed_utc"] = datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            )
            receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        print(
            f"STREAMING_TRAIT_PASS index={index}/{len(selected)} trait={trait_id} "
            f"munged_sha256={receipt['munged_output_sha256']}"
        )

    print(f"STREAMING_FAMILY_PASS traits={len(selected)} compressed_gib={required_gib:.6f}")


if __name__ == "__main__":
    main()
