#!/usr/bin/env python3
"""Audit a complete frozen SUPERGNOVA family and join protected candidates."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import pathlib
import re
import sys


REPO = pathlib.Path(__file__).resolve().parents[2]
RESULTS = REPO / "brain6/results/brain6_alternative_local_validation_v1"
sys.path.insert(0, str(REPO / "brain6/scripts"))
from run_supergnova_family_v1 import audit_raw_output, load_protocol, sha256  # noqa: E402


def read_tsv(path: pathlib.Path) -> list[dict]:
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv_new(path: pathlib.Path, rows: list[dict], fieldnames: list[str]) -> None:
    if path.exists():
        raise FileExistsError(path)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fieldnames,
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def parse_region(region: str) -> tuple[int, int, int]:
    chrom, interval = region.split(":", 1)
    start, end = interval.split("-")
    return int(chrom.removeprefix("chr")), int(start), int(end)


def parse_candidate(candidate_locus_id: str) -> tuple[int, int, int]:
    match = re.search(r"_chr(\d+)_(\d+)_(\d+)$", candidate_locus_id)
    if not match:
        raise ValueError(f"Unrecognized frozen candidate locus ID: {candidate_locus_id}")
    return tuple(map(int, match.groups()))


def to_float(s: str) -> float | None:
    if s in ("", "NA", "nan", "NaN"):
        return None
    x = float(s)
    return x if math.isfinite(x) else None


def make_tables() -> dict:
    protocol = load_protocol()
    threshold = protocol["alpha_block"]
    root = pathlib.Path(protocol["ssd_root"])
    lock_path = RESULTS / "source_eligibility_lock_v1.json"
    lock = json.loads(lock_path.read_text())
    if lock["parent_protocol_sha256"] != sha256(RESULTS / "protocol_freeze.json"):
        raise RuntimeError("Eligibility lock/protocol hash mismatch")
    if set(lock["pair_status"]) != set(protocol["pairs"]):
        raise RuntimeError("Eligibility lock does not cover exact frozen family")
    if (lock["family_denominator_pairs"] != len(protocol["pairs"]) or
            lock["blocks_per_pair"] != protocol["n_blocks_nonmhc"] or
            lock["alpha_block"] != threshold):
        raise RuntimeError("Eligibility lock changed the frozen correction")
    partitions = read_tsv(pathlib.Path(protocol["partition_file"]))
    all_blocks: list[dict] = []
    pair_blocks: dict[str, list[dict]] = {}
    receipts = []
    for pair in protocol["pairs"]:
        raw = root / "runs_v1" / f"{pair}.supergnova.txt"
        log = root / "runs_v1" / f"{pair}.log"
        receipt_path = root / "runs_v1" / f"{pair}.receipt.json"
        decision = lock["pair_status"][pair]
        if decision == "METHOD_INAPPLICABLE":
            if raw.exists() or log.exists() or receipt_path.exists():
                raise RuntimeError(f"Inapplicable pair has method artifacts: {pair}")
            rows = [{"pair_id": pair, "chr": r["chr"], "start": r["start"],
                     "end": r["end"], "rho": "NA", "corr": "NA",
                     "h2_1": "NA", "h2_2": "NA", "var": "NA", "p": "NA",
                     "m": "NA", "analysis_status": "METHOD_INAPPLICABLE",
                     "p_fwer": "", "covariance_sign": "",
                     "candidate_support_status": "METHOD_INAPPLICABLE",
                     "source_limitations": lock["inapplicable_reason"],
                     "source_eligibility": decision}
                    for r in partitions]
            pair_blocks[pair] = rows
            all_blocks.extend(rows)
            continue
        if decision != "ELIGIBLE_SECONDARY":
            raise RuntimeError(f"Unknown eligibility decision: {pair}={decision}")
        if not (raw.is_file() and log.is_file() and receipt_path.is_file()):
            raise RuntimeError(f"Incomplete pair: {pair}")
        receipt = json.loads(receipt_path.read_text())
        if (receipt["protocol_sha256"] != sha256(RESULTS / "protocol_freeze.json") or
                receipt["output_sha256"] != sha256(raw) or
                receipt["log_sha256"] != sha256(log)):
            raise RuntimeError(f"Broken pair receipt: {pair}")
        if sum(audit_raw_output(raw, protocol).values()) != protocol["n_blocks_nonmhc"]:
            raise RuntimeError(f"Incomplete block family: {pair}")
        receipts.append({"path": str(receipt_path), "sha256": sha256(receipt_path)})
        rows = []
        with raw.open() as stream:
            for source in csv.DictReader(stream, delimiter=" "):
                p = to_float(source["p"])
                rho = to_float(source["rho"])
                status = source["analysis_status"]
                if status == "ESTIMATED":
                    classification = "SUPPORTED_FWER" if p is not None and p <= threshold else "NOT_SUPPORTED"
                else:
                    classification = "NOT_TESTED_INPUT_QC"
                row = {"pair_id": pair, **source,
                       "p_fwer": "" if p is None else format(min(1.0, p * len(protocol["pairs"]) * protocol["n_blocks_nonmhc"]), ".12g"),
                       "covariance_sign": "" if rho is None else ("POSITIVE" if rho > 0 else "NEGATIVE" if rho < 0 else "ZERO"),
                       "candidate_support_status": classification,
                       "source_limitations": lock["eligible_limitations"][pair],
                       "source_eligibility": decision}
                rows.append(row)
        pair_blocks[pair] = rows
        all_blocks.extend(rows)

    block_fields = ["pair_id", "chr", "start", "end", "rho", "corr", "h2_1", "h2_2",
                    "var", "p", "m", "analysis_status", "p_fwer", "covariance_sign",
                    "candidate_support_status", "source_limitations", "source_eligibility"]
    write_tsv_new(RESULTS / "genome_wide_blocks.tsv", all_blocks, block_fields)

    candidates = read_tsv(REPO / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv")
    if len(candidates) != 25 or len({r["candidate_locus_id"] for r in candidates}) != 25:
        raise RuntimeError("Protected pair candidate cardinality changed")
    candidate_rows = []
    for source in candidates:
        pair = source["pair_id"]
        chrom, start, end = parse_candidate(source["candidate_locus_id"])
        overlaps = [r for r in pair_blocks[pair]
                    if int(r["chr"]) == chrom and int(r["start"]) <= end and int(r["end"]) >= start]
        if not overlaps:
            raise RuntimeError(f"Candidate outside frozen blocks: {source['candidate_locus_id']}")
        estimated = [r for r in overlaps if r["analysis_status"] == "ESTIMATED"]
        hits = [r for r in estimated if r["candidate_support_status"] == "SUPPORTED_FWER"]
        best = min(estimated, key=lambda r: float(r["p"])) if estimated else None
        if lock["pair_status"][pair] == "METHOD_INAPPLICABLE":
            classification = "METHOD_INAPPLICABLE"
        elif hits:
            classification = "SUPPORTED_FWER"
        elif estimated:
            classification = "NOT_SUPPORTED"
        else:
            classification = "NOT_TESTED_INPUT_QC"
        candidate_rows.append({
            "candidate_locus_id": source["candidate_locus_id"],
            "pair_id": pair, "region_group": source["region_group"],
            "lead_variants": source["lead_variants"],
            "overlapping_blocks": len(overlaps), "estimated_blocks": len(estimated),
            "supported_blocks": len(hits), "best_block": "" if best is None else
                f'{best["chr"]}:{best["start"]}-{best["end"]}',
            "best_p": "" if best is None else best["p"],
            "best_p_fwer": "" if best is None else best["p_fwer"],
            "covariance_sign": "" if best is None else best["covariance_sign"],
            "status": classification,
            "source_limitations": lock["inapplicable_reason"] if classification == "METHOD_INAPPLICABLE" else lock["eligible_limitations"][pair],
            "source_eligibility": lock["pair_status"][pair],
            "primary_lava_status": "FAILED_QC_NOT_PROMOTED",
        })
    candidate_fields = list(candidate_rows[0])
    write_tsv_new(RESULTS / "candidate_crosswalk.tsv", candidate_rows, candidate_fields)

    regions = read_tsv(REPO / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/geographic_region_groups.tsv")
    if len(regions) != 20:
        raise RuntimeError("Protected geographic region cardinality changed")
    region_rows = []
    for region in regions:
        members = [r for r in candidate_rows if r["region_group"] == region["region_group"]]
        if len(members) != int(region["n_pair_specific_candidates"]):
            raise RuntimeError(f"Candidate/region mismatch: {region['region_group']}")
        statuses = {r["status"] for r in members}
        status = ("SUPPORTED_FWER" if "SUPPORTED_FWER" in statuses else
                  "NOT_SUPPORTED" if "NOT_SUPPORTED" in statuses else
                  "NOT_TESTED_INPUT_QC" if "NOT_TESTED_INPUT_QC" in statuses else
                  "METHOD_INAPPLICABLE")
        region_rows.append({"region_group": region["region_group"],
                            "pairs": region["pairs"],
                            "candidate_locus_ids": region["candidate_locus_ids"],
                            "n_pair_specific_candidates": region["n_pair_specific_candidates"],
                            "n_supported_candidates": sum(r["status"] == "SUPPORTED_FWER" for r in members),
                            "n_method_inapplicable_candidates": sum(r["status"] == "METHOD_INAPPLICABLE" for r in members),
                            "status": status, "primary_lava_status": "FAILED_QC_NOT_PROMOTED"})
    write_tsv_new(RESULTS / "region_crosswalk.tsv", region_rows, list(region_rows[0]))

    report = {
        "protocol_sha256": sha256(RESULTS / "protocol_freeze.json"),
        "method": "SUPERGNOVA secondary observed-scale local genetic covariance",
        "primary_lava_status": protocol["primary_lava_status"],
        "n_pairs": len(protocol["pairs"]),
        "n_blocks_per_pair": protocol["n_blocks_nonmhc"],
        "n_family_cells": len(all_blocks),
        "n_candidates": len(candidate_rows), "n_geographic_regions": len(region_rows),
        "alpha_block": threshold, "pair_receipts": receipts,
        "source_eligibility_lock_sha256": sha256(lock_path),
        "pair_eligibility_status": lock["pair_status"],
        "inapplicable_reason": lock["inapplicable_reason"],
        "output_sha256": {name: sha256(RESULTS / name) for name in
                          ("genome_wide_blocks.tsv", "candidate_crosswalk.tsv", "region_crosswalk.tsv")},
        "block_status_counts": {status: sum(r["analysis_status"] == status for r in all_blocks)
                                for status in ("ESTIMATED", "NOT_TESTED_LOW_SNP", "NOT_TESTED_LOW_RANK", "METHOD_INAPPLICABLE")},
        "candidate_status_counts": {status: sum(r["status"] == status for r in candidate_rows)
                                    for status in ("SUPPORTED_FWER", "NOT_SUPPORTED", "NOT_TESTED_INPUT_QC", "METHOD_INAPPLICABLE")},
        "region_status_counts": {status: sum(r["status"] == status for r in region_rows)
                                 for status in ("SUPPORTED_FWER", "NOT_SUPPORTED", "NOT_TESTED_INPUT_QC", "METHOD_INAPPLICABLE")},
    }
    path = RESULTS / "run_provenance.json"
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    print(json.dumps(make_tables(), indent=2, sort_keys=True))
