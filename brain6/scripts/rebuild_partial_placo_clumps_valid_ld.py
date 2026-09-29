"""Rebuild partial PLACO candidate clumps after excluding invalid raw LD edges.

This is a new, noncanonical diagnostic output family. It applies the frozen
candidate threshold and clumping rules to exact-position reference matches,
but an edge whose raw LAVA R is outside [-1, 1] is recorded and excluded from
clumping without clipping. Existing partial candidate files are untouched.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import clump_placo_candidates as clump

RULE_PATH = ROOT / "brain6/config/shared_locus_rule_v1.json"
FAMILY_PATH = ROOT / "brain6/config/placo_family_v3/family_lock.json"
RUN_PATH = ROOT / "brain6/manifests/placo_family_v3_run.json"
MASTER_PATH = ROOT / "brain6/results/placo/placo_master.tsv"
REFERENCE_PREFIX = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
OUTPUT_DIR = ROOT / "brain6/results/loci"
VARIANT_PATH = OUTPUT_DIR / "placo_candidate_variants_valid_ld_partial_v2.tsv"
LOCUS_PATH = OUTPUT_DIR / "placo_candidate_loci_valid_ld_partial_v2.tsv"
EXCEPTION_PATH = OUTPUT_DIR / "placo_candidate_ld_range_exceptions_partial_v2.tsv"
PROVENANCE_PATH = OUTPUT_DIR / "placo_candidate_valid_ld_partial_v2.provenance.json"
EXCEPTION_FIELDS = ["pair_id", "CHR", "SNP1", "BP1", "SNP2", "BP2", "lava_r_raw", "disposition"]
SUMMARY_PATH = ROOT / "brain6/results/supplement/table_S22_partial_placo_valid_ld_rebuild.tsv"
SUMMARY_FIELDS = ["pair_id", "candidate_rows", "exact_reference_matches", "unmatched_candidates",
                  "original_lead_variants", "rebuild_valid_ld_lead_variants", "original_intervals",
                  "rebuild_intervals", "out_of_range_pair_edges", "interpretation"]
PRIOR_PARTIAL_OUTPUTS = {
    "candidate_variants": OUTPUT_DIR / "placo_candidate_variants_partial.tsv",
    "candidate_loci": OUTPUT_DIR / "placo_candidate_loci_partial.tsv",
    "candidate_loci_provenance": OUTPUT_DIR / "placo_candidate_loci_partial.provenance.json",
    "signed_ld": OUTPUT_DIR / "placo_candidate_signed_ld_partial.tsv",
    "signed_ld_provenance": OUTPUT_DIR / "placo_candidate_signed_ld_partial.provenance.json",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.name}.tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temp.replace(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def extract_windowed_ld(candidates: list[dict], reference_prefix: Path,
                        rscript: Path, window_bp: int) -> tuple[dict, dict, int]:
    with tempfile.TemporaryDirectory(prefix="brain6-placo-valid-ld-") as temp:
        temp_path = Path(temp)
        candidate_path, output = temp_path / "candidates.tsv", temp_path / "ld"
        unique: dict[tuple[int, str], dict] = {}
        for row in candidates:
            key = int(row["CHR"]), str(row["SNP"]).lower()
            previous = unique.get(key)
            if previous and int(previous["BP"]) != int(row["BP"]):
                raise ValueError(f"Conflicting candidate positions for {key}")
            unique[key] = row
        with candidate_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["SNP", "CHR", "BP"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows({"SNP": snp, "CHR": chr_id, "BP": row["BP"]}
                             for (chr_id, snp), row in sorted(unique.items()))
        helper = ROOT / "brain6/scripts/extract_placo_windowed_ld.R"
        proc = subprocess.run([str(rscript), str(helper), str(reference_prefix), str(candidate_path),
                               str(output), str(window_bp)],
                              capture_output=True, text=True, check=False)
        if proc.returncode:
            raise RuntimeError("Windowed pinned-LAVA LD extraction failed: " + (proc.stdout + proc.stderr)[-4000:])
        ref_info: dict[tuple[int, str], dict] = {}
        with (output / "reference_matches.tsv").open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                ref_info[(int(row["CHR"]), row["SNP"].lower())] = {
                    "BP": int(row["POS"]), "A1": row["A1"], "A2": row["A2"]}
        ld: dict[tuple[str, str], float] = {}
        with (output / "ld_pairs_within_window.tsv").open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                a, b = row["SNP1"].lower(), row["SNP2"].lower()
                key = tuple(sorted((a, b)))
                raw_r = float(row["R"])
                if not math.isfinite(raw_r):
                    raise ValueError(f"Nonfinite LAVA LD value for {key}")
                if a == b and abs(raw_r - 1) > 1e-8:
                    raise ValueError(f"Invalid LAVA LD diagonal for {a}: {raw_r}")
                if key in ld and ld[key] != raw_r:
                    raise ValueError(f"Duplicate inconsistent LAVA LD edge: {key}")
                ld[key] = raw_r
        ref_n = int((output / "reference_sample_size.txt").read_text().strip())
        return ref_info, ld, ref_n


def split_ld_edges(ld: dict[tuple[str, str], float]) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], float]]:
    valid, invalid = {}, {}
    for edge, raw_r in ld.items():
        (valid if abs(raw_r) <= 1 else invalid)[edge] = raw_r
    return valid, invalid


def exception_rows(candidates: dict[str, list[dict]], ref_info: dict,
                   invalid: dict[tuple[str, str], float], window_bp: int) -> list[dict[str, Any]]:
    output = []
    for pair, rows in sorted(candidates.items()):
        eligible = {row["SNP"].lower(): row for row in rows
                    if (int(row["CHR"]), row["SNP"].lower()) in ref_info
                    and ref_info[(int(row["CHR"]), row["SNP"].lower())]["BP"] == int(row["BP"])}
        for (a, b), raw_r in invalid.items():
            if a == b or a not in eligible or b not in eligible:
                continue
            ra, rb = eligible[a], eligible[b]
            if int(ra["CHR"]) != int(rb["CHR"]) or abs(int(ra["BP"]) - int(rb["BP"])) > window_bp:
                continue
            output.append({"pair_id": pair, "CHR": int(ra["CHR"]), "SNP1": a,
                           "BP1": int(ra["BP"]), "SNP2": b, "BP2": int(rb["BP"]),
                           "lava_r_raw": f"{raw_r:.17g}",
                           "disposition": "EXCLUDED_FROM_R2_CLUMPING_RAW_R_OUTSIDE_CORRELATION_RANGE"})
    return sorted(output, key=lambda x: (x["pair_id"], x["CHR"], x["BP1"], x["SNP1"], x["BP2"], x["SNP2"]))


def atomic_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    import os
    import tempfile as tempfile_module
    fd, name = tempfile_module.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        with tmp.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def build(reference_prefix: Path = REFERENCE_PREFIX, rscript: Path = RSCRIPT,
          output_dir: Path = OUTPUT_DIR) -> dict[str, Any]:
    rule = json.loads(RULE_PATH.read_text(encoding="utf-8"))
    family = json.loads(FAMILY_PATH.read_text(encoding="utf-8"))
    run = json.loads(RUN_PATH.read_text(encoding="utf-8"))
    candidates, input_hashes = clump.read_candidates(MASTER_PATH, float(family["headline_threshold"]))
    all_rows = [row for pair in sorted(candidates) for row in candidates[pair]]
    if not all_rows:
        raise ValueError("No admitted PLACO headline candidates")
    window_bp = int(rule["clumping"]["window_kb"]) * 1000
    ref_info, raw_ld, ref_n = extract_windowed_ld(all_rows, reference_prefix, rscript, window_bp)
    if ref_n != int(rule["clumping"]["reference_sample_size"]):
        raise ValueError(f"Reference N mismatch: {ref_n}")
    valid_ld, invalid_ld = split_ld_edges(raw_ld)
    exceptions = exception_rows(candidates, ref_info, invalid_ld, window_bp)
    variants, loci = clump.build_outputs(candidates, ref_info, valid_ld, rule)
    variant_path = output_dir / VARIANT_PATH.name
    locus_path = output_dir / LOCUS_PATH.name
    exception_path = output_dir / EXCEPTION_PATH.name
    provenance_path = output_dir / PROVENANCE_PATH.name
    summary_path = SUMMARY_PATH
    original_variants = read_tsv(OUTPUT_DIR / "placo_candidate_variants_partial.tsv")
    original_loci = read_tsv(OUTPUT_DIR / "placo_candidate_loci_partial.tsv")
    summary = []
    for pair in sorted(candidates):
        pair_variants = [row for row in variants if row["pair_id"] == pair]
        pair_original_variants = [row for row in original_variants if row["pair_id"] == pair]
        pair_loci = [row for row in loci if row["pair_id"] == pair]
        pair_original_loci = [row for row in original_loci if row["pair_id"] == pair]
        pair_exceptions = [row for row in exceptions if row["pair_id"] == pair]
        unique_edges = {(int(row["CHR"]), *sorted((row["SNP1"], row["SNP2"])))
                        for row in pair_exceptions}
        exact = sum(row["reference_status"] == "EXACT_MATCH" for row in pair_variants)
        summary.append({"pair_id": pair, "candidate_rows": len(pair_variants),
                        "exact_reference_matches": exact, "unmatched_candidates": len(pair_variants) - exact,
                        "original_lead_variants": sum(row["candidate_status"] == "LEAD" for row in pair_original_variants),
                        "rebuild_valid_ld_lead_variants": sum(row["candidate_status"] == "LEAD" for row in pair_variants),
                        "original_intervals": len(pair_original_loci), "rebuild_intervals": len(pair_loci),
                        "out_of_range_pair_edges": len(unique_edges),
                        "interpretation": "Partial PLACO candidate grouping only; invalid raw LD edges excluded without clipping"})
    clump.atomic_tsv(variant_path, clump.VARIANT_COLUMNS, variants)
    clump.atomic_tsv(locus_path, clump.LOCUS_COLUMNS, loci)
    atomic_tsv(exception_path, EXCEPTION_FIELDS, exceptions)
    atomic_tsv(summary_path, SUMMARY_FIELDS, summary)
    reference_provenance = reference_prefix.parent / "reference.provenance.json"
    reference_files = {}
    for chrom in sorted({int(row["CHR"]) for row in all_rows}):
        for suffix in ("info", "bcor"):
            path = reference_prefix.parent / f"{reference_prefix.name}_chr{chrom}.{suffix}"
            reference_files[str(path)] = sha256(path)
    exact_candidates = sum((int(row["CHR"]), row["SNP"]) in ref_info
                           and ref_info[(int(row["CHR"]), row["SNP"])]["BP"] == int(row["BP"])
                           for row in variants)
    family_complete = all(run["pairs"][p]["status"] == "COMPLETE" for p in run["pairs"]) and "insomnia__adhd" in candidates
    outputs = []
    for path in (variant_path, locus_path, exception_path, summary_path):
        outputs.append({"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size})
    result = {
        "schema_version": 1,
        "analysis_id": "brain6-placo-partial-candidate-valid-ld-v2",
        "status": "PASS_PARTIAL_WITH_LD_RANGE_EXCEPTIONS" if exceptions else "PASS_PARTIAL_FAMILY",
        "interpretation": "Recomputed partial PLACO candidate grouping from the four admitted pair outputs. Raw LAVA LD outside [-1,1] was preserved in a separate edge-QC table and excluded from r-squared grouping without clipping. This is not a complete five-track family, does not establish independent/causal loci, and lacks trait-effect-allele alignment.",
        "scope": {"pairs": sorted(candidates), "pair_count": len(candidates), "candidate_variants": len(variants),
                  "candidate_intervals": len(loci), "reference_exact_matches": exact_candidates,
                  "raw_ld_edges_within_window": len(raw_ld), "valid_ld_edges": len(valid_ld),
                  "out_of_range_ld_edges": len(invalid_ld), "pair_candidate_exception_rows": len(exceptions),
                  "protected_track_b_present": "insomnia__adhd" in candidates, "complete_five_track_family": family_complete},
        "method": {"reference": rule["clumping"]["reference"], "sample_size": ref_n,
                   "clumping_rule": {"algorithm": rule["clumping"]["algorithm"], "r2_threshold": rule["clumping"]["r2_threshold"],
                                     "window_kb": rule["clumping"]["window_kb"], "candidate_threshold": family["headline_threshold"]},
                   "raw_ld_source": "Pinned LAVA 0.1.5 read.ld; official UKB v1.1 .bcor/.info",
                   "correlation_range": [-1, 1], "out_of_range_policy": "Preserve and exclude; never clamp",
                   "effect_allele_alignment": "NOT_AVAILABLE_IN_PUBLISHED_PLACO_VARIANT_OUTPUT"},
        "sources": {"master": {"path": str(MASTER_PATH.relative_to(ROOT)), "sha256": sha256(MASTER_PATH)},
                    "family_lock": {"path": str(FAMILY_PATH.relative_to(ROOT)), "sha256": sha256(FAMILY_PATH)},
                    "run_manifest": {"path": str(RUN_PATH.relative_to(ROOT)), "sha256": sha256(RUN_PATH)},
                    "locus_rule": {"path": str(RULE_PATH.relative_to(ROOT)), "sha256": sha256(RULE_PATH)},
                    "pair_outputs": input_hashes,
                    "reference_provenance": {"path": str(reference_provenance), "sha256": sha256(reference_provenance)},
                    "reference_files": reference_files,
                    "extractor": {"path": "brain6/scripts/extract_placo_windowed_ld.R", "sha256": sha256(ROOT / "brain6/scripts/extract_placo_windowed_ld.R")},
                    "builder": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256(Path(__file__).resolve())},
                    "rscript": {"path": str(rscript), "sha256": sha256(rscript)},
                    "signed_ld_audit": {"path": "brain6/results/loci/placo_candidate_signed_ld_partial.provenance.json", "sha256": sha256(ROOT / "brain6/results/loci/placo_candidate_signed_ld_partial.provenance.json")}},
        "outputs": outputs,
    }
    result["sources"]["preserved_partial_outputs"] = {
        name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for name, path in PRIOR_PARTIAL_OUTPUTS.items()
    }
    atomic_json(provenance_path, result)
    result["provenance_sha256"] = sha256(provenance_path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-prefix", type=Path, default=REFERENCE_PREFIX)
    parser.add_argument("--rscript", type=Path, default=RSCRIPT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    print(json.dumps(build(args.reference_prefix, args.rscript, args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
