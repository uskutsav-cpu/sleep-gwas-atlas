#!/usr/bin/env python3
"""Report same-family Bonferroni sensitivity for the frozen 84-pair latent RG family."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


SOURCE = "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
OUTPUT = "frailty_paper/analysis/latent_factor_correction_sensitivity.tsv"
MANIFEST = "frailty_paper/analysis/latent_factor_correction_sensitivity.manifest.json"
EXPECTED_SLEEP = {
    "insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype", "sleepiness",
    "napping", "snoring", "sleep_apnea", "sleep_efficiency", "accel_sleep_duration",
    "sleep_timing",
}
EXPECTED_FACTORS = {
    "frailty_general", "frailty_factor_1", "frailty_factor_2", "frailty_factor_3",
    "frailty_factor_4", "frailty_factor_5", "frailty_factor_6",
}
DENOMINATOR = 84
ALPHA = 0.05
OUTPUT_COLUMNS = [
    "sleep_trait", "frailty_factor", "rg", "se", "raw_p", "bh_q_secondary_84",
    "bonferroni_family_size", "bonferroni_p_secondary_84", "bh_significant_at_0.05",
    "bonferroni_significant_at_0.05", "comparison_interpretation", "analysis_family",
    "cohort_overlap_status", "exact_participant_overlap", "source_table",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def derive(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    keys = [(r.get("sleep_trait", ""), r.get("disease_trait", "")) for r in rows]
    expected = {(sleep, factor) for sleep in EXPECTED_SLEEP for factor in EXPECTED_FACTORS}
    if len(rows) != DENOMINATOR or len(set(keys)) != DENOMINATOR or set(keys) != expected:
        raise ValueError("Latent-factor source must contain exactly the locked 12 x 7 unique pair family")

    result: list[dict[str, str]] = []
    for row in rows:
        if row.get("analysis_family") != "secondary_sleep_x_latent_frailty":
            raise ValueError("Unexpected analysis family in latent-factor source")
        if int(row.get("family_denominator", "0")) != DENOMINATOR:
            raise ValueError("Latent-factor source denominator differs from the locked 84-pair family")
        p = float(row["p_value"])
        q = float(row["q_value"])
        rg = float(row["rg"])
        se = float(row["se"])
        if not (math.isfinite(p) and 0.0 <= p <= 1.0 and math.isfinite(q) and 0.0 <= q <= 1.0):
            raise ValueError("Non-finite or out-of-domain latent-factor p/q value")
        if not (math.isfinite(rg) and math.isfinite(se) and se > 0):
            raise ValueError("Non-finite or out-of-domain latent-factor association value")
        bonferroni_p = min(p * DENOMINATOR, 1.0)
        bh_sig = q <= ALPHA
        bonf_sig = bonferroni_p <= ALPHA
        if bh_sig and bonf_sig:
            interpretation = "RETAINS_BH_AND_BONFERRONI"
        elif bh_sig:
            interpretation = "BH_ONLY_NOT_BONFERRONI"
        elif bonf_sig:
            interpretation = "BONFERRONI_ONLY_NOT_BH"
        else:
            interpretation = "NEITHER_BH_NOR_BONFERRONI"
        result.append({
            "sleep_trait": row["sleep_trait"],
            "frailty_factor": row["disease_trait"],
            "rg": row["rg"],
            "se": row["se"],
            "raw_p": row["p_value"],
            "bh_q_secondary_84": row["q_value"],
            "bonferroni_family_size": str(DENOMINATOR),
            "bonferroni_p_secondary_84": format(bonferroni_p, ".12g"),
            "bh_significant_at_0.05": str(bh_sig).upper(),
            "bonferroni_significant_at_0.05": str(bonf_sig).upper(),
            "comparison_interpretation": interpretation,
            "analysis_family": row["analysis_family"],
            "cohort_overlap_status": row["cohort_overlap_status"],
            "exact_participant_overlap": "UNKNOWN",
            "source_table": SOURCE,
        })
    return result


def encode(rows: list[dict[str, str]]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def build(repo: Path) -> dict[str, object]:
    source = repo / SOURCE
    script = Path(__file__).resolve()
    source_rows = list(csv.DictReader(source.open(encoding="utf-8", newline=""), delimiter="\t"))
    rows = derive(source_rows)
    output = repo / OUTPUT
    manifest_path = repo / MANIFEST
    payload = encode(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    manifest = {
        "schema_version": 1,
        "analysis": "same-family Bonferroni sensitivity for secondary latent-factor LDSC",
        "command": "python frailty_paper/scripts/40_build_latent_correction_sensitivity.py --repo <repo-root>",
        "source_file": SOURCE,
        "source_sha256": sha256(source),
        "script_file": str(script.relative_to(repo)),
        "script_sha256": sha256(script),
        "output_file": OUTPUT,
        "output_sha256": hashlib.sha256(payload).hexdigest(),
        "family_denominator": DENOMINATOR,
        "alpha": ALPHA,
        "rows": len(rows),
        "bh_significant": sum(r["bh_significant_at_0.05"] == "TRUE" for r in rows),
        "bonferroni_significant": sum(r["bonferroni_significant_at_0.05"] == "TRUE" for r in rows),
        "interpretation": "Bonferroni is an alternate correction sensitivity; the predeclared BH q-values remain primary.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=".")
    args = parser.parse_args()
    manifest = build(Path(args.repo).resolve())
    print(json.dumps({k: manifest[k] for k in ("rows", "bh_significant", "bonferroni_significant", "output_file", "output_sha256")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
