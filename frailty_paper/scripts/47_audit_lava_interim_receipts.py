#!/usr/bin/env python3
"""Reproduce a bounded, read-only diagnostic of interim LAVA receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
RUN_ROOT = Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1")
LOCK = REPO / "frailty_paper/config/lava_frailty_sensitivity_v1.yaml"
WORKER = REPO / "frailty_paper/scripts/44_run_lava_sensitivity_locus.R"
COLLATOR = REPO / "frailty_paper/scripts/46_collate_lava_sensitivity.py"
LOCUS_LIMIT = 2495
FAILURE_LIMIT = 0.01


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_receipt(value: dict, log_text: str, trait: str) -> tuple[bool, dict[str, int], str | None]:
    """Apply the collator predicate and distinguish known from unexpected failures."""
    process_status = value.get("process_status")
    if process_status not in {"PROCESSED", "PROCESS_FAILED"}:
        raise ValueError(f"unknown process status: {process_status}")
    for key in ("fi_status", "sleep_status"):
        if key not in value:
            raise ValueError(f"receipt lacks {key}")

    process_failed = process_status == "PROCESS_FAILED"
    fi_failed = value["fi_status"] != "TESTED"
    sleep_failed = value["sleep_status"] != "TESTED"
    counts: Counter[str] = Counter()
    unexpected_process_message = None
    if process_failed:
        counts["process_failed"] += 1
        if ("process.locus returned NULL" in str(value.get("error", ""))
                and "Negative variance estimate for all phenotypes" in log_text):
            counts["all_phenotypes_negative_variance"] += 1
        elif ("process.locus returned NULL" in str(value.get("error", ""))
              and "none of specified SNP IDs are present in reference data" in log_text):
            counts["no_specified_snps_in_reference"] += 1
        else:
            counts["unexpected_process_failure"] += 1
            error_lines = [line.strip() for line in log_text.splitlines()
                           if "Error:" in line or "error:" in line]
            receipt_error = str(value.get("error", "")).strip()
            log_error = error_lines[-1] if error_lines else "no explicit error line in log"
            unexpected_process_message = f"receipt: {receipt_error or 'missing error'}; log: {log_error}"
    else:
        if fi_failed:
            counts["frailty_univariate_not_tested"] += 1
            if "phenotype(s) 'frailty'" not in log_text:
                raise ValueError("frailty status lacks negative-variance diagnostic")
            counts["frailty_negative_variance"] += 1
        if sleep_failed:
            counts["sleep_univariate_not_tested"] += 1
            if f"phenotype(s) '{trait}'" not in log_text:
                raise ValueError("sleep status lacks negative-variance diagnostic")
            counts["sleep_negative_variance"] += 1
    return process_failed or fi_failed or sleep_failed, dict(counts), unexpected_process_message


def latest_complete_prefix(pair_root: Path, trait: str, locus_limit: int = LOCUS_LIMIT) -> int:
    """Return the last contiguous locus with an identity-matched receipt and nonempty log."""
    latest = 0
    for index in range(1, locus_limit + 1):
        receipt = pair_root / "loci" / f"locus_{index:04d}.json"
        log = pair_root / "logs" / f"locus_{index:04d}.log"
        if not receipt.is_file() or not log.is_file() or log.stat().st_size == 0:
            break
        value = json.loads(receipt.read_text(encoding="utf-8"))
        if value.get("trait") != trait or value.get("locus_index") != index:
            raise ValueError(f"receipt identity mismatch at locus index {index}")
        if value.get("process_status") not in {"PROCESSED", "PROCESS_FAILED"}:
            raise ValueError(f"receipt is not complete at locus index {index}")
        if "fi_status" not in value or "sleep_status" not in value:
            raise ValueError(f"receipt lacks per-trait status at locus index {index}")
        latest = index
    return latest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trait", default="accel_sleep_duration")
    parser.add_argument("--through-index", required=True,
                        help="inclusive prefix index, or 'latest' to detect contiguous completed receipts")
    parser.add_argument("--run-root", type=Path, default=RUN_ROOT)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    pair_root = args.run_root / "pairs" / args.trait
    selection_method = "explicit"
    if args.through_index == "latest":
        through_index = latest_complete_prefix(pair_root, args.trait)
        if through_index == 0:
            parser.error("no contiguous completed receipt/log prefix found")
        selection_method = "latest_contiguous_nonempty_receipt_log_pair"
    else:
        try:
            through_index = int(args.through_index)
        except ValueError:
            parser.error("--through-index must be an integer or 'latest'")
        if not 1 <= through_index <= LOCUS_LIMIT:
            parser.error(f"--through-index must be from 1 to {LOCUS_LIMIT}")

    counts: Counter[str] = Counter()
    unexpected_process_messages: Counter[str] = Counter()
    unexpected_process_failures: list[dict[str, str | int]] = []
    flagged: set[int] = set()
    manifest = hashlib.sha256()
    for index in range(1, through_index + 1):
        receipt = pair_root / "loci" / f"locus_{index:04d}.json"
        log = pair_root / "logs" / f"locus_{index:04d}.log"
        if not receipt.is_file() or not log.is_file():
            raise FileNotFoundError(f"receipt or log missing at locus index {index}")
        value = json.loads(receipt.read_text(encoding="utf-8"))
        if value.get("trait") != args.trait or value.get("locus_index") != index:
            raise ValueError(f"receipt identity mismatch at index {index}")
        log_text = log.read_text(encoding="utf-8", errors="replace")
        for path in (receipt, log):
            rel = path.relative_to(args.run_root).as_posix()
            manifest.update(f"{rel}\t{sha256(path)}\n".encode())

        is_flagged, locus_counts, unexpected_process_message = classify_receipt(value, log_text, args.trait)
        counts.update(locus_counts)
        if unexpected_process_message:
            unexpected_process_messages[unexpected_process_message] += 1
            unexpected_process_failures.append({
                "locus_index": index,
                "LOC": value.get("LOC"),
                "message": unexpected_process_message,
            })
        if is_flagged:
            flagged.add(index)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "trait": args.trait,
        "through_locus_index_inclusive": through_index,
        "receipts_audited": through_index,
        "checkpoint_selection": selection_method,
        "loci_flagged_by_frozen_collator_rule": len(flagged),
        "locus_family_denominator": LOCUS_LIMIT,
        "locus_failure_fraction": len(flagged) / LOCUS_LIMIT,
        "locked_locus_failure_fraction_limit": FAILURE_LIMIT,
        "locked_maximum_failed_loci": int(LOCUS_LIMIT * FAILURE_LIMIT),
        "gate_already_exceeded": len(flagged) / LOCUS_LIMIT > FAILURE_LIMIT,
        "counts": dict(sorted(counts.items())),
        "unexpected_process_failure_messages": dict(sorted(unexpected_process_messages.items())),
        "unexpected_process_failures": unexpected_process_failures,
        "receipt_and_log_manifest_sha256": manifest.hexdigest(),
        "analysis_lock_sha256": sha256(LOCK),
        "worker_script_sha256": sha256(WORKER),
        "collator_script_sha256": sha256(COLLATOR),
        "auditor_script_sha256": sha256(Path(__file__)),
        "external_run_root": str(args.run_root),
        "mutates_run_or_receipts": False,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
