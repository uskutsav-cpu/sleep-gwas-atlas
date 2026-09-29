#!/usr/bin/env python3
"""Build pair-level cohort-overlap annotations for frozen aging-context RG rows."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv"
LEDGER = ROOT / "frailty_analysis/manifests/cohort_overlap.tsv"
OUTPUT = ROOT / "frailty_paper/analysis/frozen_atlas_aging_context_overlap.tsv"
MANIFEST = ROOT / "frailty_paper/manifests/frozen_atlas_aging_context_overlap_manifest_v2.json"
COHORTS = ("uk_biobank", "finngen", "23andme", "charge")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pair_status(a: str, b: str) -> str:
    a, b = a.strip().upper(), b.strip().upper()
    if a not in {"YES", "NO", "UNKNOWN"} or b not in {"YES", "NO", "UNKNOWN"}:
        raise ValueError(f"unexpected cohort flag pair: {a!r}, {b!r}")
    if a == b == "YES":
        return "SHARED_COHORT_EXPECTED"
    if a == b == "NO":
        return "NO_SHARED_COHORT_REPORTED"
    return "POSSIBLE_OR_UNKNOWN"


def build_rows(result_rows: list[dict[str, str]], ledger_rows: list[dict[str, str]],
               expected_count: int = 72) -> list[dict[str, str]]:
    by_trait = {row["trait_id"]: row for row in ledger_rows}
    if len(by_trait) != len(ledger_rows):
        raise ValueError("cohort ledger has duplicate trait_id rows")
    output = []
    for result in result_rows:
        sleep, context = result["sleep_trait"], result["non_sleep_trait"]
        if sleep not in by_trait or context not in by_trait:
            raise ValueError(f"cohort ledger missing source row for {sleep} × {context}")
        srow, crow = by_trait[sleep], by_trait[context]
        row = {
            "pair_id": result["pair_id"],
            "sleep_trait": sleep,
            "context_trait": context,
            "sleep_source_id": srow["source_id"],
            "context_source_id": crow["source_id"],
            "exact_participant_intersection": "UNKNOWN",
            "interpretation_status": "COHORT_LEVEL_ASSESSMENT_ONLY",
        }
        shared = []
        all_no = True
        for cohort in COHORTS:
            a, b = srow[cohort], crow[cohort]
            status = pair_status(a, b)
            row[f"sleep_{cohort}"] = a
            row[f"context_{cohort}"] = b
            row[f"{cohort}_pair_status"] = status
            if status == "SHARED_COHORT_EXPECTED":
                shared.append(cohort)
            if status != "NO_SHARED_COHORT_REPORTED":
                all_no = False
        row["declared_shared_cohorts"] = ",".join(shared) if shared else "NONE_IDENTIFIED"
        row["overall_cohort_status"] = (
            "SHARED_COHORT_EXPECTED" if shared else
            "NO_COMMON_COHORT_DECLARED; EXACT OVERLAP STILL UNKNOWN" if all_no else
            "POSSIBLE_OR_UNKNOWN; EXACT OVERLAP UNKNOWN"
        )
        row["source_overlap_notes"] = (
            "Flags reflect source-level cohort declarations only; UNKNOWN never means zero overlap. "
            "See the cited cohort ledger evidence for cohort composition."
        )
        output.append(row)
    if len(output) != expected_count or len({row["pair_id"] for row in output}) != expected_count:
        raise ValueError(f"expected {expected_count} unique pair rows, received {len(output)}")
    return output


def publish_idempotently(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise SystemExit(f"REFUSING_NONIDENTICAL_OVERWRITE: {path.relative_to(ROOT)}")
        return
    with path.open("xb") as handle:
        handle.write(payload)


def main() -> int:
    result_rows, ledger_rows = read_tsv(RESULTS), read_tsv(LEDGER)
    rows = build_rows(result_rows, ledger_rows)
    fields = list(rows[0])
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    payload = output.getvalue().encode()
    publish_idempotently(OUTPUT, payload)
    manifest = {
        "artifact_role": "COHORT_LEVEL_OVERLAP_CROSSWALK_FOR_READ_ONLY_FROZEN_AGING_CONTEXT_RG",
        "result_extract_sha256": sha256(RESULTS),
        "cohort_ledger_sha256": sha256(LEDGER),
        "pair_count": len(rows),
        "cohort_axes": list(COHORTS),
        "exact_participant_intersections": "UNKNOWN for every pair",
        "shared_cohort_expected_pair_count": sum(r["overall_cohort_status"] == "SHARED_COHORT_EXPECTED" for r in rows),
        "possible_or_unknown_pair_count": sum(r["overall_cohort_status"].startswith("POSSIBLE_OR_UNKNOWN") for r in rows),
        "no_common_cohort_declared_pair_count": sum(r["overall_cohort_status"].startswith("NO_COMMON_COHORT") for r in rows),
        "output": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(payload).hexdigest(),
        "claim_limit": "Cohort-level source declarations do not establish exact participant overlap or independence.",
        "reproducibility": {
            "generator": str(Path(__file__).resolve().relative_to(ROOT)),
            "generator_sha256": sha256(Path(__file__).resolve()),
            "command": "python3 frailty_paper/scripts/30_build_aging_context_overlap.py",
            "python_version": sys.version.split()[0],
        },
    }
    publish_idempotently(MANIFEST, (json.dumps(manifest, indent=2) + "\n").encode())
    print("AGING_CONTEXT_OVERLAP_OK pairs=72 shared_expected=" + str(manifest["shared_cohort_expected_pair_count"]) +
          " possible_or_unknown=" + str(manifest["possible_or_unknown_pair_count"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
