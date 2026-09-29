"""Validate the isolated, diagnostic PLACO factor-normalization sensitivity."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2"
EXPECTED_PAIRS = {
    "insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz"
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_legacy_v1_integrity(root: Path) -> dict[str, str]:
    """Check frozen v1 bytes while v2 becomes the current code-bound output."""
    legacy = root / "brain6/results/loci/placo_factor_normalized_sensitivity_v1"
    provenance = json.loads((legacy / "provenance.json").read_text(encoding="utf-8"))
    if (provenance.get("analysis_id") != "brain6_placo_factor_normalized_sensitivity_v1" or
            provenance.get("status") != "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED"):
        raise ValueError("archived factor-normalized v1 identity/status changed")
    checksums: dict[str, str] = {}
    for record in provenance.get("outputs", []):
        path = root / record["path"]
        if (not path.is_file() or path.stat().st_size != record["bytes"] or
                sha256(path) != record["sha256"]):
            raise ValueError(f"archived factor-normalized v1 output changed: {path}")
        checksums[record["name"]] = record["sha256"]
    if set(checksums) != {"variants", "loci", "replaced_edges", "summary"}:
        raise ValueError("archived factor-normalized v1 output manifest is incomplete")
    return checksums


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_range_edges(rows: list[dict[str, str]]) -> int:
    """Validate diagnostic replacements without allowing raw invalid LD through."""
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        key = row["CHR"], row["SNP1"], row["SNP2"]
        if key in seen:
            raise ValueError(f"duplicate normalized-range edge: {key}")
        seen.add(key)
        raw = float(row["lava_r_unscaled"])
        normalized = float(row["factor_r_normalized"])
        normalized_r2 = float(row["factor_r2_normalized"])
        if not math.isfinite(raw) or abs(raw) <= 1.0:
            raise ValueError(f"range-exception row is not raw out-of-range LD: {key}")
        if not math.isfinite(normalized) or abs(normalized) > 1.0:
            raise ValueError(f"normalized diagnostic LD is outside [-1,1]: {key}")
        if not math.isclose(normalized_r2, normalized * normalized, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"normalized r-squared is inconsistent: {key}")
        if row["diagnostic_status"] != "NORMALIZED_FACTOR_SENSITIVITY_ONLY_NOT_GENOTYPE_VALIDATED":
            raise ValueError(f"unexpected promotion/interpretation status: {key}")
    return len(rows)


def validate_loci_match_baseline(original: list[dict[str, str]],
                                 sensitivity: list[dict[str, str]]) -> dict[str, int]:
    """Report whether factor normalization preserves baseline leads and intervals."""
    def interval_keys(rows: list[dict[str, str]]) -> set[tuple[str, str, str, str]]:
        return {(r["pair_id"], r["CHR"], r["START"], r["STOP"]) for r in rows}

    def lead_keys(rows: list[dict[str, str]]) -> set[tuple[str, str]]:
        return {
            (row["pair_id"], lead)
            for row in rows
            for lead in row["lead_variants"].split(";")
            if lead
        }

    original_intervals, sensitivity_intervals = interval_keys(original), interval_keys(sensitivity)
    original_leads, sensitivity_leads = lead_keys(original), lead_keys(sensitivity)
    return {
        "baseline_interval_count": len(original_intervals),
        "sensitivity_interval_count": len(sensitivity_intervals),
        "matching_interval_count": len(original_intervals & sensitivity_intervals),
        "baseline_lead_count": len(original_leads),
        "sensitivity_lead_count": len(sensitivity_leads),
        "matching_lead_count": len(original_leads & sensitivity_leads),
        "interval_sets_match": int(original_intervals == sensitivity_intervals),
        "lead_sets_match": int(original_leads == sensitivity_leads),
    }


def validate(root: Path = ROOT) -> dict[str, Any]:
    legacy_v1_output_sha256 = validate_legacy_v1_integrity(root)
    out = root / "brain6/results/loci/placo_factor_normalized_sensitivity_v2"
    provenance_path = out / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    scope = provenance["scope"]
    if provenance.get("analysis_id") != "brain6_placo_factor_normalized_sensitivity_v2":
        raise ValueError("unexpected PLACO factor-normalized sensitivity identity")
    if provenance.get("status") != "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED" or provenance["method"].get("promotion") != "NONE":
        raise ValueError("diagnostic sensitivity cannot be promoted by this validator")
    if provenance["method"].get("external_genotype_validation") != "NOT_AVAILABLE":
        raise ValueError("external genotype-validation status must remain explicit")
    if set(scope["pairs"]) != EXPECTED_PAIRS or scope["full_five_pair_family"] is not False:
        raise ValueError("unexpected PLACO pair scope")
    if (scope["candidate_rows"], scope["candidate_variants"], scope["candidate_intervals"],
            scope["range_exception_edges_normalized"], scope["reference_N"]) != (2297, 2297, 19, 12317, 100000):
        raise ValueError("PLACO factor-normalized scope differs from the audited four-pair family")
    if (scope["invalid_cross_block_edges"] != 0 or scope["same_block_edges_replaced"] != 157403
            or scope["primary_blocks_touched"] != 45 or scope["roundoff_boundary_values_normalized"] != 1922):
        raise ValueError("factor-normalized block/roundoff diagnostics differ from the audited scope")

    output_rows: dict[str, dict[str, str]] = {}
    for record in provenance["outputs"]:
        path = root / record["path"]
        if not path.is_file() or path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"sensitivity output checksum mismatch: {path}")
        output_rows[record["name"]] = {"path": str(path), "sha256": record["sha256"]}
    required = {"variants", "loci", "replaced_edges", "summary"}
    if set(output_rows) != required:
        raise ValueError("sensitivity output manifest is incomplete")
    if legacy_v1_output_sha256 != {
            name: record["sha256"] for name, record in output_rows.items()}:
        raise ValueError("v2 result tables differ from the preserved v1 result bytes")

    sources = provenance["sources"]
    for name in ("master", "builder", "family_lock", "run_manifest", "locus_rule",
                 "prior_valid_ld_provenance", "reference_provenance",
                 "factor_reader_audit_script", "factor_psd_audit"):
        record = sources[name]
        path = Path(record["path"])
        if not path.is_absolute():
            path = root / path
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ValueError(f"sensitivity source checksum mismatch ({name}): {path}")
    master_rows = read_tsv(root / "brain6/results/placo/placo_master.tsv")
    master_by_pair = {row["pair_id"]: row for row in master_rows}
    if set(sources["candidate_input_hashes"]) != EXPECTED_PAIRS:
        raise ValueError("candidate-input checksum ledger has an unexpected pair set")
    for pair, expected_hash in sources["candidate_input_hashes"].items():
        if master_by_pair[pair]["output_sha256"] != expected_hash:
            raise ValueError(f"candidate-input hash differs from PLACO master for {pair}")
        input_path = root / master_by_pair[pair]["output_path"]
        if not input_path.is_file() or sha256(input_path) != expected_hash:
            raise ValueError(f"candidate-input checksum mismatch for {pair}: {input_path}")
    for path_string, record in sources["reference_files"].items():
        path = Path(path_string)
        if (not path.is_file() or path.stat().st_size != record["bytes"]
                or len(record["sha256"]) != 64):
            raise ValueError(f"reference file is unavailable or differs in size: {path}")

    variants = read_tsv(Path(output_rows["variants"]["path"]))
    loci = read_tsv(Path(output_rows["loci"]["path"]))
    exceptions = read_tsv(Path(output_rows["replaced_edges"]["path"]))
    summary = read_tsv(Path(output_rows["summary"]["path"]))
    if len(variants) != 2297 or len({(r["pair_id"], r["SNP"]) for r in variants}) != len(variants):
        raise ValueError("candidate variant rows are incomplete or duplicated")
    reference_status = {row["reference_status"] for row in variants}
    if (sum(row["reference_status"] == "EXACT_MATCH" for row in variants) != 2273
            or sum(row["reference_status"] == "NO_REFERENCE_VARIANT" for row in variants) != 24
            or reference_status != {"EXACT_MATCH", "NO_REFERENCE_VARIANT"}):
        raise ValueError("exact-position reference coverage differs from the candidate input audit")
    if len(loci) != 19 or {r["pair_id"] for r in loci} != EXPECTED_PAIRS:
        raise ValueError("candidate locus output is incomplete or has unexpected pairs")
    if any(r["locus_status"] != "PLACO_ONLY_CANDIDATE_LOCUS" for r in loci):
        raise ValueError("candidate loci have an unexpected evidence status")
    if validate_range_edges(exceptions) != 12317:
        raise ValueError("normalized raw-range exception count is incomplete")
    if {r["pair_id"] for r in summary} != EXPECTED_PAIRS or sum(int(r["candidate_rows"]) for r in summary) != 2297:
        raise ValueError("pair summary does not account for the four-pair candidate set")
    if sum(int(r["factor_normalized_sensitivity_intervals"]) for r in summary) != 19:
        raise ValueError("pair summary does not account for all sensitivity intervals")

    baseline = read_tsv(root / "brain6/results/loci/placo_candidate_loci_partial.tsv")
    match = validate_loci_match_baseline(baseline, loci)
    if match["interval_sets_match"] != 1 or match["lead_sets_match"] != 1:
        raise ValueError("factor-normalized sensitivity changed the baseline lead or interval set")
    return {
        "status": "PASS_DIAGNOSTIC_ONLY_NO_PROMOTION",
        "candidate_variants": len(variants),
        "candidate_intervals": len(loci),
        "normalized_range_exception_edges": len(exceptions),
        "reference_blocks_touched": scope["primary_blocks_touched"],
        "roundoff_boundary_values_normalized": scope["roundoff_boundary_values_normalized"],
        "baseline_comparison": match,
        "legacy_v1_output_sha256": legacy_v1_output_sha256,
        "output_sha256": {name: record["sha256"] for name, record in output_rows.items()},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="print machine-readable validation result")
    args = parser.parse_args()
    result = validate()
    print(json.dumps(result, indent=2) if args.json else result["status"])
