#!/usr/bin/env python3
"""Compare frozen all-396 BH FI results with a stricter same-family Bonferroni check."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import sys
from pathlib import Path


ALPHA = 0.05
EXPECTED_TRAITS = {
    "insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype",
    "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency",
    "accel_sleep_duration", "sleep_timing",
}
REQUIRED = {
    "pair_id", "sleep_trait", "non_sleep_trait", "global_rg", "global_rg_se",
    "global_rg_p", "global_rg_fdr_all_396",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_rows(source_rows: list[dict[str, str]], denominator: int) -> list[dict[str, str]]:
    if denominator != 396:
        raise ValueError(f"Expected inherited all-396 denominator, got {denominator}")
    if len(source_rows) != 12 or {row.get("sleep_trait") for row in source_rows} != EXPECTED_TRAITS:
        raise ValueError("Frozen FI extract does not contain exactly the 12 locked sleep traits")
    output = []
    seen_pairs = set()
    for row in source_rows:
        missing = REQUIRED - row.keys()
        if missing:
            raise ValueError(f"Missing required input columns: {sorted(missing)}")
        if row["non_sleep_trait"] != "frailty":
            raise ValueError(f"Unexpected non-sleep trait in FI extract: {row['non_sleep_trait']}")
        if row["pair_id"] in seen_pairs:
            raise ValueError(f"Duplicate pair_id: {row['pair_id']}")
        seen_pairs.add(row["pair_id"])
        p_value = float(row["global_rg_p"])
        q_value = float(row["global_rg_fdr_all_396"])
        rg = float(row["global_rg"])
        se = float(row["global_rg_se"])
        if not (0 <= p_value <= 1 and 0 <= q_value <= 1 and se >= 0):
            raise ValueError(f"Invalid P/q/SE for {row['pair_id']}")
        adjusted = min(1.0, p_value * denominator)
        output.append({
            "pair_id": row["pair_id"],
            "sleep_trait": row["sleep_trait"],
            "non_sleep_trait": row["non_sleep_trait"],
            "global_rg": row["global_rg"],
            "global_rg_se": row["global_rg_se"],
            "global_rg_p_frozen": row["global_rg_p"],
            "bh_q_all_396_frozen": row["global_rg_fdr_all_396"],
            "bh_all_396_significant_at_0.05": str(q_value <= ALPHA).upper(),
            "bonferroni_family_size": str(denominator),
            "bonferroni_p_all_396": format(adjusted, ".12g"),
            "bonferroni_all_396_significant_at_0.05": str(adjusted <= ALPHA).upper(),
            "comparison_interpretation": "SAME_FROZEN_P_VALUES; stricter correction sensitivity only; primary BH q-values unchanged",
        })
    return output


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    repo = args.repo.resolve()
    source_manifest_path = repo / "frailty_paper/manifests/frozen_atlas_frailty_manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    if source_manifest.get("artifact_role") != "READ_ONLY_EXTRACT_OF_FROZEN_CORE_ATLAS_RESULTS_NOT_NEW_ANALYSIS":
        raise SystemExit("Unexpected source manifest role")
    denominator_text = source_manifest.get("fdr_denominator", "")
    match = re.fullmatch(r"(\d+) locked sleep-by-non-sleep pairs", denominator_text)
    if not match:
        raise SystemExit(f"Cannot parse frozen correction denominator: {denominator_text!r}")
    denominator = int(match.group(1))
    source = repo / source_manifest["extracted_table"]
    source_hash = sha256(source)
    if source_hash != source_manifest.get("extracted_table_sha256"):
        raise SystemExit("Frozen FI result table does not match its audit-manifest hash")
    rows = build_rows(read_tsv(source), denominator)
    output = repo / "frailty_paper/analysis/frozen_fi_correction_sensitivity.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "artifact_role": "DERIVED_CORRECTION_SENSITIVITY_FROM_FROZEN_P_VALUES_NOT_NEW_GWAS_ANALYSIS",
        "input": source.relative_to(repo).as_posix(),
        "input_sha256": source_hash,
        "input_audit_manifest": source_manifest_path.relative_to(repo).as_posix(),
        "input_audit_manifest_sha256": sha256(source_manifest_path),
        "locked_family_size": denominator,
        "primary_correction": "BH q-values over all 396 frozen sleep-by-non-sleep tests, inherited unchanged",
        "sensitivity_correction": "Bonferroni adjusted P = min(raw P * 396, 1), alpha 0.05",
        "rows": len(rows),
        "bh_significant_count": sum(r["bh_all_396_significant_at_0.05"] == "TRUE" for r in rows),
        "bonferroni_significant_count": sum(r["bonferroni_all_396_significant_at_0.05"] == "TRUE" for r in rows),
        "script": Path(__file__).resolve().relative_to(repo).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "command": "python " + Path(__file__).resolve().relative_to(repo).as_posix() + (" " + " ".join(sys.argv[1:]) if sys.argv[1:] else ""),
        "python_version": platform.python_version(),
        "output_sha256": sha256(output),
        "interpretation_limit": "The Bonferroni column is an explicitly secondary robustness description. It does not replace inherited BH q-values, add independent evidence, or remedy absent original harmonized/munged inputs.",
    }
    sidecar = output.with_suffix(".manifest.json")
    sidecar.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "rows": len(rows), "BH": manifest["bh_significant_count"], "Bonferroni": manifest["bonferroni_significant_count"]}))


if __name__ == "__main__":
    main()
