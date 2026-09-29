#!/usr/bin/env python3
"""Describe MHC-excluded counts in the frozen partial PLACO candidate set."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "brain6/config/placo_candidate_mhc_exclusion_v1.json"
SCRIPT = Path(__file__).resolve()
OUT = ROOT / "brain6/results/sensitivity/placo_candidate_mhc_exclusion.tsv"
PROVENANCE = OUT.with_suffix(".provenance.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_mhc(chromosome: str, position: int, *, chr6: str, start: int, end: int) -> bool:
    return chromosome.removeprefix("chr") == chr6 and start <= position < end


def atomic_write(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def analyze() -> dict:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    paths = {key: ROOT / value for key, value in plan["source_paths"].items()}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    family = json.loads(paths["family_lock"].read_text(encoding="utf-8"))
    provenance = json.loads(paths["candidate_provenance"].read_text(encoding="utf-8"))
    if family.get("n_selected_tracks") != 5:
        raise ValueError("The frozen PLACO family no longer specifies five tracks")
    if family.get("across_pair_adjustment") is None:
        raise ValueError("Frozen PLACO across-track correction is absent")
    if provenance.get("pairs_included") != plan["expected_pairs"] or provenance.get("family_complete") is not False:
        raise ValueError("Candidate provenance does not match the frozen partial-family scope")
    if provenance.get("candidate_threshold") != 1e-8:
        raise ValueError("Candidate threshold differs from the frozen PLACO sensitivity plan")

    variant_path = paths["candidate_variants"]
    locus_path = paths["candidate_loci"]
    expected_outputs = {variant_path: None, locus_path: None}
    for item in provenance.get("outputs", []):
        recorded = ROOT / item["path"]
        if recorded in expected_outputs:
            expected_outputs[recorded] = item["sha256"]
    for path, expected in expected_outputs.items():
        if expected is None or sha256(path) != expected:
            raise ValueError(f"Candidate source does not match frozen provenance: {path}")

    coords = plan["excluded_interval"]
    chr6, start, end = coords["chromosome"], int(coords["start_inclusive"]), int(coords["end_exclusive"])
    pair_stats: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    variant_rows = read_tsv(variant_path)
    expected_pairs = set(plan["expected_pairs"])
    for row in variant_rows:
        pair = row["pair_id"]
        if pair not in expected_pairs:
            raise ValueError(f"Unexpected PLACO pair in candidate variants: {pair}")
        stat = pair_stats[pair]
        stat["candidate_variants_all"] += 1
        if row["reference_status"] != "EXACT_MATCH":
            stat["candidate_variants_unpositioned_or_mismatched"] += 1
            continue
        stat["candidate_variants_position_validated"] += 1
        if is_mhc(row["CHR"], int(row["BP"]), chr6=chr6, start=start, end=end):
            stat["candidate_variants_inside_mhc"] += 1
        else:
            stat["candidate_variants_outside_mhc"] += 1

    locus_rows = read_tsv(locus_path)
    interval_stats: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in locus_rows:
        pair = row["pair_id"]
        if pair not in expected_pairs:
            raise ValueError(f"Unexpected PLACO pair in candidate loci: {pair}")
        stat = interval_stats[pair]
        stat["candidate_intervals_all"] += 1
        overlaps = (row["CHR"].removeprefix("chr") == chr6 and
                    int(row["START"]) < end and int(row["STOP"]) > start)
        stat["candidate_intervals_overlapping_mhc"] += int(overlaps)
        stat["candidate_intervals_wholly_outside_mhc"] += int(not overlaps)

    fields = ["pair_id", "candidate_variants_all", "candidate_variants_position_validated",
              "candidate_variants_inside_mhc", "candidate_variants_outside_mhc",
              "candidate_variants_unpositioned_or_mismatched", "candidate_intervals_all",
              "candidate_intervals_overlapping_mhc", "candidate_intervals_wholly_outside_mhc"]
    output_rows = []
    for pair in plan["expected_pairs"]:
        values = {**pair_stats[pair], **interval_stats[pair]}
        output_rows.append({field: pair if field == "pair_id" else str(values.get(field, 0)) for field in fields})
    all_values = defaultdict(int)
    for row in output_rows:
        if row["pair_id"] != "pair_id":
            for field in fields[1:]:
                all_values[field] += int(row[field])
    output_rows.append({field: "ALL_PARTIAL_PAIRS" if field == "pair_id" else str(all_values[field])
                        for field in fields})
    import io
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
    payload = buffer.getvalue()
    atomic_write(OUT, payload)
    details = {
        "analysis_id": plan["analysis_id"],
        "status": "PASS_PARTIAL_FAMILY_DESCRIPTIVE_COUNTS",
        "plan_path": str(PLAN.relative_to(ROOT)),
        "plan_sha256": sha256(PLAN),
        "script_path": str(SCRIPT.relative_to(ROOT)),
        "script_sha256": sha256(SCRIPT),
        "python_version": sys.version.split()[0],
        "coordinate_build": plan["coordinate_build"],
        "excluded_interval": coords,
        "interval_reference": plan["interval_reference"],
        "candidate_threshold_unchanged": provenance["candidate_threshold"],
        "pair_scope": plan["expected_pairs"],
        "family_complete": False,
        "protected_track_b_present": False,
        "source_hashes": {str(path.relative_to(ROOT)): sha256(path) for path in paths.values()},
        "output_path": str(OUT.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(payload.encode("utf-8")).hexdigest(),
        "interpretation": plan["inference_policy"],
    }
    atomic_write(PROVENANCE, json.dumps(details, indent=2, sort_keys=True) + "\n")
    return details


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
