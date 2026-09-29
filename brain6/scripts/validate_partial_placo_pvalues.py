#!/usr/bin/env python3
"""Replay candidate PLACO+ p-values from frozen Z scores and parameter receipts.

The replay only validates the four published partial-family pair outputs. It
does not complete protected Track B, revise the multiple-testing family, or
promote any candidate locus.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
MASTER = ROOT / "brain6/results/placo/placo_master.tsv"
CANDIDATES = ROOT / "brain6/results/loci/placo_candidate_variants_partial.tsv"
FAMILY_LOCK = ROOT / "brain6/config/placo_family_v3/family_lock.json"
OUTPUT = ROOT / "brain6/results/loci/placo_candidate_pvalue_validation_v1.tsv"
SUMMARY_OUTPUT = ROOT / "brain6/results/supplement/table_S39_partial_placo_pvalue_replay.tsv"
PROVENANCE = ROOT / "brain6/results/loci/placo_candidate_pvalue_validation_v1.provenance.json"
R_HELPER = ROOT / "brain6/scripts/recompute_partial_placo_candidate_pvalues.R"
DEFAULT_RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
DEFAULT_PLACO_SOURCE = Path("/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R")
DEFAULT_OUTPUT_ROOT = Path("/Volumes/Extreme SSD/brain6-work/placo-family-v3-results")
RELATIVE_TOLERANCE = 1e-8
PARAMETER_FIELDS = ["pair_id", "VarZ1", "VarZ2", "CorZ", "AbsTol"]
OUTPUT_FIELDS = [
    "pair_id", "SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO_REPORTED",
    "P_PLACO_RECOMPUTED", "relative_error", "validation_status",
]
SUMMARY_FIELDS = ["pair_id", "candidate_variants", "matching_p_values", "maximum_relative_error",
                  "relative_error_tolerance", "interpretation"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def rel_error(reported: float, recomputed: float) -> float:
    if reported == recomputed:
        return 0.0
    return abs(reported - recomputed) / max(abs(reported), abs(recomputed), 1e-300)


def is_match(reported: float, recomputed: float, tolerance: float = RELATIVE_TOLERANCE) -> bool:
    return math.isfinite(reported) and math.isfinite(recomputed) and rel_error(reported, recomputed) <= tolerance


def parameter_record(pair_id: str, output_root: Path, placo_source: Path,
                     rscript: Path) -> tuple[dict[str, str], dict[str, str]]:
    pair_root = output_root / pair_id
    parameter_path = pair_root / "parameters/parameters.json"
    receipt_path = pair_root / "parameters/receipt.json"
    document = json.loads(receipt_path.read_text(encoding="utf-8"))
    if document.get("status") != "COMPLETE" or document.get("scientific_status") not in {"PASS", "SUCCESS"}:
        raise ValueError(f"parameter receipt is not successful for {pair_id}")
    recorded = {Path(row["path"]).name: row for row in document.get("outputs", [])}
    record = recorded.get("parameters.json")
    if not record or record.get("sha256") != sha256(parameter_path):
        raise ValueError(f"parameter receipt does not bind parameters for {pair_id}")
    parameters = json.loads(parameter_path.read_text(encoding="utf-8"))
    if parameters.get("method") != "PLACO_PLUS":
        raise ValueError(f"unexpected method in {pair_id} parameter file")
    receipt_inputs = {Path(item["path"]).name: item for item in document.get("inputs", [])}
    for name, path in (("PLACO_v0.2.0.R", placo_source), (rscript.name, rscript)):
        item = receipt_inputs.get(name)
        if not item or item.get("sha256") != sha256(path):
            raise ValueError(f"parameter receipt does not bind {name} for {pair_id}")
    row = {
        "pair_id": pair_id,
        "VarZ1": repr(float(parameters["VarZ"][0])),
        "VarZ2": repr(float(parameters["VarZ"][1])),
        "CorZ": repr(float(parameters["CorZ"])),
        "AbsTol": repr(float(1e-13)),
    }
    provenance = {
        "parameters": {"path": str(parameter_path), "sha256": sha256(parameter_path)},
        "receipt": {"path": str(receipt_path), "sha256": sha256(receipt_path)},
    }
    return row, provenance


def collect_candidate_zscores(master_rows: list[dict[str, str]], candidate_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, str]]:
    wanted: dict[str, dict[tuple[str, str, str], dict[str, str]]] = {}
    for row in candidate_rows:
        wanted.setdefault(row["pair_id"], {})[(row["SNP"], row["CHR"], row["BP"])] = row
    selected: list[dict[str, str]] = []
    source_hashes: dict[str, str] = {}
    for master in master_rows:
        pair_id = master["pair_id"]
        if pair_id not in wanted:
            continue
        path = ROOT / master["output_path"]
        if not path.is_file() or sha256(path) != master["output_sha256"]:
            raise ValueError(f"published PLACO output hash mismatch for {pair_id}")
        source_hashes[master["output_path"]] = master["output_sha256"]
        remaining = dict(wanted[pair_id])
        with gzip.open(path, "rb") as stream:
            header = stream.readline().rstrip(b"\r\n").split(b"\t")
            columns = {name.decode("ascii"): index for index, name in enumerate(header)}
            required = {"SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO", "status"}
            if not required.issubset(columns):
                raise ValueError(f"unexpected PLACO output schema for {pair_id}")
            while line := stream.readline():
                row = line.rstrip(b"\r\n").split(b"\t")
                key = (row[columns["SNP"]].decode(), row[columns["CHR"]].decode(), row[columns["BP"]].decode())
                if key not in remaining:
                    continue
                candidate = remaining.pop(key)
                if row[columns["status"]] != b"TESTED":
                    raise ValueError(f"candidate {pair_id}/{key[0]} is not TESTED")
                selected.append({
                    "pair_id": pair_id, "SNP": key[0], "CHR": key[1], "BP": key[2],
                    "Z1": row[columns["Z1"]].decode(), "Z2": row[columns["Z2"]].decode(),
                    "P_PLACO_REPORTED": row[columns["P_PLACO"]].decode(),
                })
        if remaining:
            raise ValueError(f"{len(remaining)} candidate variants were not found in published PLACO output for {pair_id}")
    if len(selected) != len(candidate_rows):
        raise ValueError(f"candidate coverage mismatch: expected {len(candidate_rows)}, found {len(selected)}")
    return selected, source_hashes


def validate(rscript: Path = DEFAULT_RSCRIPT, placo_source: Path = DEFAULT_PLACO_SOURCE,
             output_root: Path = DEFAULT_OUTPUT_ROOT) -> dict[str, Any]:
    if not rscript.is_file() or not placo_source.is_file():
        raise FileNotFoundError("Pinned PLACO runtime or source is unavailable")
    master_rows = read_tsv(MASTER)
    candidate_rows = read_tsv(CANDIDATES)
    selected, pair_hashes = collect_candidate_zscores(master_rows, candidate_rows)
    params, parameter_sources = [], {}
    for pair_id in sorted({row["pair_id"] for row in selected}):
        parameter, sources = parameter_record(pair_id, output_root, placo_source, rscript)
        params.append(parameter)
        parameter_sources[pair_id] = sources

    with tempfile.TemporaryDirectory(prefix="brain6-placo-p-replay-") as tmp:
        tmp_path = Path(tmp)
        input_path, params_path, output_path = tmp_path / "candidates.tsv", tmp_path / "parameters.tsv", tmp_path / "recomputed.tsv"
        for path, fields, rows in (
            (input_path, ["pair_id", "SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO_REPORTED"], selected),
            (params_path, PARAMETER_FIELDS, params),
        ):
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        proc = subprocess.run([str(rscript), str(R_HELPER), str(input_path), str(params_path),
                               str(placo_source), str(output_path)], capture_output=True, text=True, check=False)
        if proc.returncode:
            raise RuntimeError("Pinned PLACO candidate p-value replay failed: " + (proc.stdout + proc.stderr)[-5000:])
        recomputed_rows = read_tsv(output_path)

    if len(recomputed_rows) != len(selected):
        raise ValueError("R replay row count differs from candidate input")
    output_rows, pair_summary = [], {}
    for row in recomputed_rows:
        reported, recomputed = float(row["P_PLACO_REPORTED"]), float(row["P_PLACO_RECOMPUTED"])
        passed = is_match(reported, recomputed)
        output_rows.append({
            "pair_id": row["pair_id"], "SNP": row["SNP"], "CHR": row["CHR"], "BP": row["BP"],
            "Z1": row["Z1"], "Z2": row["Z2"], "P_PLACO_REPORTED": row["P_PLACO_REPORTED"],
            "P_PLACO_RECOMPUTED": row["P_PLACO_RECOMPUTED"],
            "relative_error": f"{rel_error(reported, recomputed):.12g}",
            "validation_status": "PASS" if passed else "FAIL",
        })
        summary = pair_summary.setdefault(row["pair_id"], {"variants": 0, "passed": 0, "max_relative_error": 0.0})
        summary["variants"] += 1
        summary["passed"] += int(passed)
        summary["max_relative_error"] = max(summary["max_relative_error"], rel_error(reported, recomputed))
    if any(row["validation_status"] != "PASS" for row in output_rows):
        raise ValueError("One or more candidate PLACO p-values differ from pinned-source replay")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    summary_rows = [{
        "pair_id": pair_id,
        "candidate_variants": summary["variants"],
        "matching_p_values": summary["passed"],
        "maximum_relative_error": f"{summary['max_relative_error']:.12g}",
        "relative_error_tolerance": f"{RELATIVE_TOLERANCE:.12g}",
        "interpretation": "Exact pinned-source p-value replay; partial PLACO evidence only, no locus promotion",
    } for pair_id, summary in sorted(pair_summary.items())]
    SUMMARY_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY_OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_partial_placo_candidate_pvalue_replay_v1",
        "status": "PASS_PINNED_SOURCE_CANDIDATE_PVALUES_RECOMPUTED",
        "scope": {
            "candidate_variants": len(output_rows), "pairs": sorted(pair_summary),
            "full_five_track_family": False, "protected_track_b_recomputed": False,
            "candidate_locus_promotion": False,
        },
        "method": {
            "replay": "Exact PLACO_PLUS p-value function from hash-pinned upstream PLACO_v0.2.0.R",
            "parameter_estimation": "Receipt-bound genome-wide VarZ and CorZ for each pair; not re-estimated",
            "relative_error_tolerance": RELATIVE_TOLERANCE,
        },
        "pair_summary": pair_summary,
        "sources": {
            "candidate_variants": {"path": str(CANDIDATES.relative_to(ROOT)), "sha256": sha256(CANDIDATES)},
            "placo_master": {"path": str(MASTER.relative_to(ROOT)), "sha256": sha256(MASTER)},
            "family_lock": {"path": str(FAMILY_LOCK.relative_to(ROOT)), "sha256": sha256(FAMILY_LOCK)},
            "pair_outputs": pair_hashes, "parameter_receipts": parameter_sources,
            "placo_source": {"path": str(placo_source), "sha256": sha256(placo_source)},
            "rscript": {"path": str(rscript), "sha256": sha256(rscript)},
            "replay_helper": {"path": str(R_HELPER.relative_to(ROOT)), "sha256": sha256(R_HELPER)},
            "builder": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256(Path(__file__).resolve())},
        },
        "output": {"path": str(OUTPUT.relative_to(ROOT)), "sha256": sha256(OUTPUT), "rows": len(output_rows)},
        "summary_output": {"path": str(SUMMARY_OUTPUT.relative_to(ROOT)),
                           "sha256": sha256(SUMMARY_OUTPUT), "rows": len(summary_rows)},
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance["provenance_sha256"] = sha256(PROVENANCE)
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    parser.add_argument("--placo-source", type=Path, default=DEFAULT_PLACO_SOURCE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    print(json.dumps(validate(args.rscript, args.placo_source, args.output_root), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
