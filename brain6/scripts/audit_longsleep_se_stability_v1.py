#!/usr/bin/env python3
"""Audit genome-wide SE scale without reading association effects or P values."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/longsumstats.txt.zip")
README = ROOT / "brain6/results/lava_confirmatory_source_triage_v1/source_docs/Saxena_fullUKBB_Longsleep_summary_stats_README"
DECISION = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json"
OUTPUT = ROOT / "brain6/results/lava_confirmatory_source_triage_v1/longsleep_se_stability_v1.json"
EXPECTED = {
    SOURCE: "0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885",
    README: "ec79d6599ad8ec0353c2d8476bd182e512b5a333519c5b6d9263bf2a9e9df0aa",
    DECISION: "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22",
}
CASES, CONTROLS = 34_184, 305_742
MIN_INFO = 0.95
MIN_MAF = 0.01
SAMPLE_STRIDE = 100
AUTOSOMES = {str(i) for i in range(1, 23)}


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quantiles(values: list[float]) -> dict[str, float | int]:
    if not values:
        return {"n": 0}
    ordered = sorted(values)
    def at(fraction: float) -> float:
        index = (len(ordered) - 1) * fraction
        left = math.floor(index)
        right = math.ceil(index)
        return ordered[left] + (ordered[right] - ordered[left]) * (index - left)
    return {"n": len(ordered), "q05": at(0.05), "q25": at(0.25),
            "median": statistics.median(ordered), "q75": at(0.75), "q95": at(0.95)}


def main() -> None:
    for path, expected in EXPECTED.items():
        if sha(path) != expected:
            raise ValueError(f"Frozen input changed: {path}")
    total_n = CASES + CONTROLS
    case_fraction = CASES / total_n
    counts: Counter[str] = Counter()
    ratios: list[float] = []
    by_chromosome: dict[str, list[float]] = defaultdict(list)
    with zipfile.ZipFile(SOURCE) as archive:
        members = [item for item in archive.infolist() if not item.is_dir()]
        if len(members) != 1 or members[0].filename != "longsumstats.txt":
            raise ValueError("Unexpected long-sleep archive members")
        with archive.open(members[0]) as raw:
            reader = csv.reader((line.decode("utf-8") for line in raw), delimiter="\t")
            header = next(reader)
            expected_header = ["SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ",
                               "INFO", "BETA_LONGSLEEP", "SE_LONGSLEEP", "P_LONGSLEEP"]
            if header != expected_header:
                raise ValueError("Unexpected original source schema")
            chrom, freq, info, se = (header.index(name) for name in
                ("CHR", "A1FREQ", "INFO", "SE_LONGSLEEP"))
            for row in reader:
                counts["source_rows"] += 1
                if len(row) != len(header):
                    counts["malformed_rows"] += 1
                    continue
                try:
                    chromosome = str(int(row[chrom]))
                except (ValueError, OverflowError):
                    counts["non_autosomal_rows"] += 1
                    continue
                try:
                    frequency = float(row[freq])
                    quality = float(row[info])
                    observed_se = float(row[se])
                except (ValueError, OverflowError):
                    counts["non_numeric_rows"] += 1
                    continue
                if not (chromosome in AUTOSOMES
                        and all(math.isfinite(x) for x in (frequency, quality, observed_se))
                        and MIN_MAF <= frequency <= 1 - MIN_MAF
                        and MIN_INFO <= quality <= 1 and observed_se > 0):
                    counts["ineligible_rows"] += 1
                    continue
                counts["eligible_rows"] += 1
                if counts["eligible_rows"] % SAMPLE_STRIDE:
                    continue
                genotype_variance = 2 * frequency * (1 - frequency) * quality
                simple_se = math.sqrt(case_fraction * (1 - case_fraction)
                                      / (total_n * genotype_variance))
                ratio = observed_se / simple_se
                if not math.isfinite(ratio) or ratio <= 0:
                    raise ValueError("Invalid SE ratio after eligibility filtering")
                ratios.append(ratio)
                by_chromosome[chromosome].append(ratio)
    if counts["malformed_rows"]:
        raise ValueError("Original source contains malformed rows")
    if len(by_chromosome) != 22 or min(map(len, by_chromosome.values())) < 100:
        raise ValueError("Insufficient sampled source coverage across chromosomes")
    report = {
        "analysis_id": "brain6_longsleep_se_stability_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": {str(path): digest for path, digest in EXPECTED.items()},
        "source_rows": counts["source_rows"],
        "eligible_rows": counts["eligible_rows"],
        "excluded_row_counts": {key: counts[key] for key in
                                ("malformed_rows", "non_autosomal_rows", "non_numeric_rows", "ineligible_rows")},
        "sampled_rows": len(ratios),
        "selection": {"min_info": MIN_INFO, "min_maf": MIN_MAF,
                      "sample_stride_among_eligible_rows": SAMPLE_STRIDE,
                      "association_effect_or_p_fields_read": False},
        "published_cases": CASES,
        "published_controls": CONTROLS,
        "observed_to_simple_linear_se_ratio": quantiles(ratios),
        "ratio_outside_0_9_to_1_1": sum(not 0.9 <= value <= 1.1 for value in ratios),
        "chromosome_ratio": {key: quantiles(value) for key, value in
                             sorted(by_chromosome.items(), key=lambda item: int(item[0]))},
        "interpretation": "Genome-wide SE scale diagnostic only. A near-constant ratio supports a stable linear-model sample-size scale but cannot recover per-variant analyzed N or validate LAVA binary reconstruction; mixed-model relatedness, covariates and residual variance remain unmodeled.",
        "canonical_modified": False,
        "lava_run_started": False,
    }
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"source_rows": counts["source_rows"],
                      "eligible_rows": counts["eligible_rows"],
                      "sampled_rows": len(ratios),
                      "median_ratio": report["observed_to_simple_linear_se_ratio"]["median"],
                      "output_sha256": sha(OUTPUT)}, sort_keys=True))


if __name__ == "__main__":
    main()
