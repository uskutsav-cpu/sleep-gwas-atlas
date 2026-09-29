#!/usr/bin/env python3
"""Read-only, receipt-backed comparison of the locked longsleep/sleepdur pilot."""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_RUN_ROOT = Path(
    "/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/"
    "frailty_v1/lava_sensitivity_v1"
)
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOCUS_FILE = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
TRAITS = ("longsleep", "sleepdur")
EXPECTED_LOCI = 2495


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_trait(run_root: Path, trait: str) -> tuple[dict[int, dict], dict[str, int], str]:
    pair_root = run_root / "pairs" / trait
    loci_dir = pair_root / "loci"
    receipts: dict[int, dict] = {}
    statuses: Counter[str] = Counter()
    manifest = hashlib.sha256()
    for receipt_path in sorted(loci_dir.glob("locus_*.json")):
        # Finder metadata files on macOS start with ._; they are not receipts.
        if receipt_path.name.startswith("._"):
            continue
        value = json.loads(receipt_path.read_text(encoding="utf-8"))
        index = value.get("locus_index")
        if not isinstance(index, int) or not 1 <= index <= EXPECTED_LOCI:
            raise ValueError(f"invalid locus index in {receipt_path}")
        if value.get("trait") != trait or receipt_path.stem != f"locus_{index:04d}":
            raise ValueError(f"receipt identity mismatch: {receipt_path}")
        if index in receipts:
            raise ValueError(f"duplicate receipt for {trait} locus {index}")
        log_path = pair_root / "logs" / f"locus_{index:04d}.log"
        if not log_path.is_file() or log_path.stat().st_size == 0:
            raise ValueError(f"missing/empty paired log for {trait} locus {index}")
        if value.get("process_status") not in {"PROCESSED", "PROCESS_FAILED"}:
            raise ValueError(f"incomplete process status for {trait} locus {index}")
        sleep_status = value.get("sleep_status")
        if sleep_status not in {"TESTED", "PHENOTYPE_DROPPED", "PROCESS_FAILED"}:
            raise ValueError(f"unexpected sleep status for {trait} locus {index}: {sleep_status}")
        receipts[index] = value
        statuses[sleep_status] += 1
        for path in (receipt_path, log_path):
            rel = path.relative_to(run_root).as_posix()
            manifest.update(f"{rel}\t{sha256(path)}\n".encode())
    return receipts, dict(sorted(statuses.items())), manifest.hexdigest()


def audit(run_root: Path, locus_file: Path = DEFAULT_LOCUS_FILE) -> dict:
    data = {trait: collect_trait(run_root, trait) for trait in TRAITS}
    overlap = sorted(set(data["longsleep"][0]) & set(data["sleepdur"][0]))
    if len(overlap) < 100:
        raise ValueError(f"pilot has only {len(overlap)} matched loci; expected at least 100")
    locus_chromosomes: dict[int, int] = {}
    with locus_file.open(encoding="utf-8") as stream:
        header = stream.readline().split()
        if header != ["LOC", "CHR", "START", "STOP"]:
            raise ValueError(f"unexpected locus-file header: {header}")
        for line in stream:
            fields = line.split()
            if len(fields) != 4:
                raise ValueError("malformed row in frozen locus file")
            loc, chrom = int(fields[0]), int(fields[1])
            if loc in locus_chromosomes:
                raise ValueError(f"duplicate locus {loc} in frozen locus file")
            locus_chromosomes[loc] = chrom
    if len(locus_chromosomes) != EXPECTED_LOCI:
        raise ValueError(f"frozen locus file contains {len(locus_chromosomes)} loci")
    try:
        chromosome_counts = Counter(locus_chromosomes[index] for index in overlap)
    except KeyError as error:
        raise ValueError(f"receipt locus is absent from frozen locus file: {error}") from error
    counts = {}
    for trait in TRAITS:
        receipts, _, _ = data[trait]
        tested = sum(receipts[i].get("sleep_status") == "TESTED" for i in overlap)
        reference_k = [receipts[i]["K"] for i in overlap
                       if isinstance(receipts[i].get("K"), int)]
        counts[trait] = {
            "matched_loci": len(overlap),
            "eligible_tested": tested,
            "not_eligible": len(overlap) - tested,
            "eligibility_rate": tested / len(overlap),
            "tested_status_counts_on_matched_loci": dict(sorted(
                Counter(receipts[i].get("sleep_status") for i in overlap).items()
            )),
            "shared_reference_K": {
                "observed_loci": len(reference_k),
                "minimum": min(reference_k) if reference_k else None,
                "median": statistics.median(reference_k) if reference_k else None,
                "maximum": max(reference_k) if reference_k else None,
                "below_locked_minimum_2": sum(value < 2 for value in reference_k),
            },
            "receipt_and_log_manifest_sha256": data[trait][2],
        }
    long_rate = counts["longsleep"]["eligibility_rate"]
    duration_rate = counts["sleepdur"]["eligibility_rate"]
    long_fail = counts["longsleep"]["not_eligible"]
    duration_fail = counts["sleepdur"]["not_eligible"]
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": str(run_root),
        "analysis_id": "frailty_fi_sleep_lava_sensitivity_v1",
        "comparison": "longsleep_vs_sleepdur",
        "matched_locus_indices": overlap,
        "matched_locus_count": len(overlap),
        "locus_file_sha256": sha256(locus_file),
        "chromosome_counts": dict(sorted(chromosome_counts.items())),
        "chromosomes_covered": len(chromosome_counts),
        "autosomes_covered": f"{len(chromosome_counts)}/22",
        "traits": counts,
        "eligibility_rate_difference_percentage_points": (duration_rate - long_rate) * 100,
        "relative_reduction_in_not_eligible": (
            (long_fail - duration_fail) / long_fail if long_fail else None
        ),
        "interpretation": "interim matched-receipt screen; assess chromosome coverage before treating as a formal pilot; not bivariate evidence",
        "mutates_run_or_receipts": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--locus-file", type=Path, default=DEFAULT_LOCUS_FILE)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.run_root, args.locus_file)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
