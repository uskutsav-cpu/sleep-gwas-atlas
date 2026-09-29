#!/usr/bin/env python3
"""Source-only SE-scale diagnostic for the locked Dashti long-sleep release."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/lava_confirmatory_source_triage_v1/longsleep_effect_scale_v1"
ARCHIVE = Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/longsumstats.txt.zip")
README = ROOT / ("brain6/results/lava_confirmatory_source_triage_v1/source_docs/"
                 "Saxena_fullUKBB_Longsleep_summary_stats_README")
DECISION = ROOT / ("work/lava-canonical-v3-production/"
    "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/"
    "canonical_family_decision.json")
EXPECTED = {
    ARCHIVE: "0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885",
    README: "ec79d6599ad8ec0353c2d8476bd182e512b5a333519c5b6d9263bf2a9e9df0aa",
    DECISION: "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22",
}
CASES, CONTROLS = 34_184, 305_742
TARGET = 1_000
MIN_INFO = 0.9


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    for path, digest in EXPECTED.items():
        if sha(path) != digest:
            raise ValueError(f"Source hash mismatch: {path}")
    n = CASES + CONTROLS
    p = CASES / n
    linear_ratios: list[float] = []
    logistic_ratios: list[float] = []
    source_rows_seen = 0
    with zipfile.ZipFile(ARCHIVE) as z:
        members = [x for x in z.namelist() if not x.startswith("__MACOSX/")]
        if len(members) != 1:
            raise ValueError("Expected one long-sleep source member")
        with z.open(members[0]) as raw:
            reader = csv.DictReader((line.decode("utf-8") for line in raw), delimiter="\t")
            required = {"A1FREQ", "INFO", "SE_LONGSLEEP"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError("Original long-sleep source schema changed")
            for row in reader:
                source_rows_seen += 1
                try:
                    frequency = float(row["A1FREQ"])
                    info = float(row["INFO"])
                    observed_se = float(row["SE_LONGSLEEP"])
                except ValueError:
                    continue
                if not all(math.isfinite(x) for x in (frequency, info, observed_se)):
                    continue
                if not 0 < frequency < 1 or not MIN_INFO <= info <= 1 or observed_se <= 0:
                    continue
                genotype_variance = 2 * frequency * (1 - frequency) * info
                expected_linear_se = math.sqrt(p * (1 - p) / (n * genotype_variance))
                expected_logistic_se = math.sqrt(1 / (n * p * (1 - p) * genotype_variance))
                linear_ratios.append(observed_se / expected_linear_se)
                logistic_ratios.append(observed_se / expected_logistic_se)
                if len(linear_ratios) == TARGET:
                    break
    if len(linear_ratios) != TARGET:
        raise ValueError("Insufficient high-INFO source rows")
    summary = {
        "schema_version": 1,
        "analysis_id": "brain6_longsleep_effect_scale_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sampling": "first 1000 source rows with finite positive SE, 0<EAF<1, and 0.9<=INFO<=1",
        "source_rows_seen": source_rows_seen,
        "eligible_rows": TARGET,
        "cases": CASES,
        "controls": CONTROLS,
        "linear_ratio_median": statistics.median(linear_ratios),
        "linear_ratio_q25": statistics.quantiles(linear_ratios, n=4)[0],
        "linear_ratio_q75": statistics.quantiles(linear_ratios, n=4)[2],
        "logistic_ratio_median": statistics.median(logistic_ratios),
        "logistic_ratio_q25": statistics.quantiles(logistic_ratios, n=4)[0],
        "logistic_ratio_q75": statistics.quantiles(logistic_ratios, n=4)[2],
        "source_sha256": {str(path): digest for path, digest in EXPECTED.items()},
        "effect_or_pvalue_fields_read": False,
        "interpretation": "SE scale supports a linear-probability association release, consistent with the study's BOLT-LMM primary model. Simple independent-sample formulas are diagnostic only; relatedness, covariates, and imputation can alter exact SE. No LAVA input model or source admission follows from this check.",
        "canonical_modified": False,
        "new_lava_run_started": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    output = OUT / "effect_scale_summary.json"
    output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"rows": TARGET, "source_rows_seen": source_rows_seen,
                      "linear_ratio_median": summary["linear_ratio_median"],
                      "logistic_ratio_median": summary["logistic_ratio_median"],
                      "output_sha256": sha(output)}, sort_keys=True))


if __name__ == "__main__":
    main()
