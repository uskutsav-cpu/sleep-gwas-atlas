#!/usr/bin/env python3
"""Fail closed if the locked atlas-v1.0 core differs from its checkpoint."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
CHECKPOINT = ROOT / "discovery_extension" / "core_checkpoint.json"


def fail(message: str) -> None:
    raise SystemExit(f"CORE_CHECKPOINT_FAIL: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", "-C", str(ROOT), *args],
        check=check,
        capture_output=True,
        text=True,
    )


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    checkpoint_commit = checkpoint["git_commit"]
    if git("cat-file", "-e", f"{checkpoint_commit}^{{commit}}", check=False).returncode:
        fail(f"checkpoint commit is unavailable: {checkpoint_commit}")
    if git("merge-base", "--is-ancestor", checkpoint_commit, "HEAD", check=False).returncode:
        fail(f"HEAD does not descend from checkpoint commit {checkpoint_commit}")

    mismatches: list[str] = []
    for relative, expected in checkpoint["artifact_hashes_sha256"].items():
        path = ROOT / relative
        if not path.is_file():
            mismatches.append(f"missing {relative}")
            continue
        observed = sha256(path)
        if observed != expected:
            mismatches.append(f"hash {relative}: expected {expected}, observed {observed}")
    if mismatches:
        fail("; ".join(mismatches))

    panel = load_tsv(ROOT / "config" / "analysis_panel.tsv")
    results = load_tsv(ROOT / checkpoint["primary_result_path"])
    sleep = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
    external = [row["trait_id"] for row in panel if row["domain"] != "sleep"]
    expected_pairs = {(sleep_id, external_id) for sleep_id in sleep for external_id in external}
    observed_pairs = {(row["sleep_trait"], row["disease_trait"]) for row in results}

    if len(panel) != checkpoint["core_trait_count"]:
        fail(f"core trait count changed: {len(panel)}")
    if len(sleep) != checkpoint["sleep_trait_count"]:
        fail(f"sleep trait count changed: {len(sleep)}")
    if len(external) != checkpoint["non_sleep_trait_count"]:
        fail(f"non-sleep trait count changed: {len(external)}")
    if len(results) != checkpoint["primary_pair_count"]:
        fail(f"primary result row count changed: {len(results)}")
    if len(observed_pairs) != len(results):
        fail("primary result contains duplicate sleep-by-trait pairs")
    if observed_pairs != expected_pairs:
        fail("primary result is not the exact locked sleep x non-sleep Cartesian product")

    tiers: dict[str, int] = {}
    for row in results:
        tier = row["analysis_tier"]
        tiers[tier] = tiers.get(tier, 0) + 1
    expected_tiers = {
        "PRIMARY_PHASE1": checkpoint["primary_phase1_pair_count"],
        "QC_FAILED_SENSITIVITY": checkpoint["qc_failed_sensitivity_pair_count"],
    }
    if tiers != expected_tiers:
        fail(f"analysis-tier counts changed: {tiers}")

    print(
        "CORE_CHECKPOINT_OK "
        f"commit={checkpoint_commit} "
        f"traits={len(panel)} pairs={len(results)} "
        f"result_sha256={checkpoint['primary_result_sha256']}"
    )


if __name__ == "__main__":
    try:
        main()
    except (KeyError, ValueError, json.JSONDecodeError) as exc:
        fail(f"invalid checkpoint or core schema: {exc}")

