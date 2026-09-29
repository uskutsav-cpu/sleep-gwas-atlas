#!/usr/bin/env python3
"""Describe pair-effect directions against valid signed LD for partial PLACO clumps.

This is a post-selection, four-pair diagnostic. It preserves out-of-range LD,
does not infer causality, and does not complete the protected five-pair family.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
from brain6.artifacts import verify_artifact

PAIRS = ("insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz")
CANDIDATE_PATH = ROOT / "brain6/results/loci/placo_candidate_variants_valid_ld_partial_v2.tsv"
CLUMP_PROVENANCE = ROOT / "brain6/results/loci/placo_candidate_valid_ld_partial_v2.provenance.json"
PAIR_JOIN_ROOT = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1")
REFERENCE_PREFIX = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
REFERENCE_PROVENANCE = REFERENCE_PREFIX.parent / "reference.provenance.json"
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
R_EXTRACTOR = ROOT / "brain6/scripts/extract_placo_candidate_lead_ld.R"
OUTPUT = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
PROVENANCE = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.provenance.json"
SUMMARY = ROOT / "brain6/results/supplement/table_S23_partial_placo_signed_effect_ld.tsv"
FIELDS = [
    "pair_id", "SNP", "CHR", "BP", "candidate_status", "lead_SNP", "locus_id",
    "beta_sleep_reference_A1", "beta_disorder_reference_A1", "direction_relation",
    "reference_A1", "reference_A2", "lead_reference_A1", "lead_reference_A2",
    "effect_orientation_status", "lava_r_raw", "signed_r_reference_A1", "r2_reference",
    "ld_qc_status", "sleep_effect_vs_lead_ld", "disorder_effect_vs_lead_ld",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def allele_sign(ref_a1: str, ref_a2: str, effect_a1: str, effect_a2: str) -> int | None:
    """Return the sign converting an effect to the reference A1 allele."""
    alleles = {"A", "C", "G", "T"}
    ref_a1, ref_a2, effect_a1, effect_a2 = [x.upper() for x in (ref_a1, ref_a2, effect_a1, effect_a2)]
    if any(x not in alleles for x in (ref_a1, ref_a2, effect_a1, effect_a2)):
        return None
    if ref_a1 == ref_a2 or effect_a1 == effect_a2:
        return None
    complement = str.maketrans("ACGT", "TGCA")
    for a1, a2, sign in ((effect_a1, effect_a2, 1), (effect_a2, effect_a1, -1),
                         (effect_a1.translate(complement), effect_a2.translate(complement), 1),
                         (effect_a2.translate(complement), effect_a1.translate(complement), -1)):
        if (a1, a2) == (ref_a1, ref_a2):
            return sign
    return None


def ld_direction_consistency(beta_variant: float, beta_lead: float, r: float) -> str:
    """Check whether marginal effect signs track signed LD to the selected lead."""
    if not all(math.isfinite(x) for x in (beta_variant, beta_lead, r)) or abs(r) > 1:
        return "NOT_ASSESSED_INVALID_LD_OR_EFFECT"
    if beta_variant == 0 or beta_lead == 0 or r == 0:
        return "NOT_ASSESSED_ZERO_EFFECT_OR_LD"
    return "CONSISTENT_WITH_LD_SIGN" if (beta_variant * beta_lead > 0) == (r > 0) else "OPPOSITE_TO_LD_SIGN"


def atomic_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        with temp.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def pair_effects(pair_id: str, wanted: set[str]) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    root = PAIR_JOIN_ROOT / f"pair_{pair_id}"
    receipt_path = root / "receipt.json"
    verified = verify_artifact(root)
    if verified.get("stage") != "join_pair" or verified.get("scientific_status") != "PASS":
        raise ValueError(f"Pair join is not receipt-verified and QC-passed: {pair_id}")
    path = root / "pair.tsv.gz"
    declared = {Path(x["path"]).name: x["sha256"] for x in verified.get("outputs", [])}
    if "pair.tsv.gz" not in declared:
        raise ValueError(f"Pair artifact receipt does not bind pair.tsv.gz: {pair_id}")
    before = path.stat()
    found: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"SNP", "CHR", "BP", "A1", "A2", "BETA1", "BETA2"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Unexpected joined pair schema for {pair_id}")
        for row in reader:
            snp = row["SNP"].lower()
            if snp in wanted:
                if snp in found:
                    raise ValueError(f"Duplicate SNP in verified pair join: {pair_id}/{snp}")
                found[snp] = row
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Joined pair input changed while being read: {pair_id}")
    return found, {"pair_tsv_path": str(path), "pair_tsv_sha256": declared["pair.tsv.gz"],
                   "pair_receipt_path": str(receipt_path), "pair_receipt_sha256": sha256(receipt_path)}


def build(output: Path = OUTPUT, provenance_path: Path = PROVENANCE,
          summary_path: Path = SUMMARY) -> dict[str, Any]:
    candidates = read_tsv(CANDIDATE_PATH)
    if not candidates or {r["pair_id"] for r in candidates} != set(PAIRS):
        raise ValueError("Expected all four admitted pairs in the partial valid-LD candidate table")
    clump_provenance = json.loads(CLUMP_PROVENANCE.read_text(encoding="utf-8"))
    candidate_output = next((x for x in clump_provenance["outputs"]
                             if x["path"] == str(CANDIDATE_PATH.relative_to(ROOT))), None)
    if candidate_output is None or sha256(CANDIDATE_PATH) != candidate_output["sha256"]:
        raise ValueError("Candidate table differs from its valid-LD clumping provenance")
    wanted_by_pair: dict[str, set[str]] = {pair: set() for pair in PAIRS}
    extraction_rows = []
    for row in candidates:
        pair, snp, lead = row["pair_id"], row["SNP"].lower(), row["lead_SNP"].lower()
        if row["reference_status"] == "EXACT_MATCH":
            wanted_by_pair[pair].update((snp, lead))
            extraction_rows.append({"pair_id": pair, "SNP": snp, "CHR": row["CHR"], "BP": row["BP"], "lead_SNP": lead})

    pair_data: dict[str, dict[str, dict[str, str]]] = {}
    pair_source_provenance = {}
    for pair in PAIRS:
        pair_data[pair], pair_source_provenance[pair] = pair_effects(pair, wanted_by_pair[pair])
        missing = wanted_by_pair[pair] - set(pair_data[pair])
        if missing:
            raise ValueError(f"Range-filtered candidates are absent from the verified pair join for {pair}: {sorted(missing)[:5]}")

    with tempfile.TemporaryDirectory(prefix="brain6-placo-effect-ld-") as temp:
        temp_path = Path(temp)
        candidate_path = temp_path / "candidate_leads.tsv"
        with candidate_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["pair_id", "SNP", "CHR", "BP", "lead_SNP"], delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(extraction_rows)
        ld_dir = temp_path / "ld"
        proc = subprocess.run([str(RSCRIPT), str(R_EXTRACTOR), str(REFERENCE_PREFIX), str(candidate_path), str(ld_dir)],
                              capture_output=True, text=True, check=False)
        if proc.returncode:
            raise RuntimeError("Pinned LAVA candidate-to-lead LD extraction failed: " + (proc.stdout + proc.stderr)[-4000:])
        reference = {(int(r["CHR"]), r["SNP"].lower()): r for r in read_tsv(ld_dir / "reference_alleles.tsv")}
        edges = {(r["pair_id"], r["SNP"].lower(), r["lead_SNP"].lower()): float(r["R"])
                 for r in read_tsv(ld_dir / "candidate_lead_ld.tsv")}
        output_rows = []
        for row in candidates:
            pair, snp, lead = row["pair_id"], row["SNP"].lower(), row["lead_SNP"].lower()
            chr_id, bp = int(row["CHR"]), int(row["BP"])
            status = "NO_EXACT_REFERENCE_OR_LEAD"
            raw_r = None
            beta_sleep = beta_disorder = None
            beta_sleep_lead = beta_disorder_lead = None
            cand_ref = reference.get((chr_id, snp)); lead_ref = reference.get((chr_id, lead))
            cand_effect = pair_data[pair].get(snp); lead_effect = pair_data[pair].get(lead)
            orientation = "NOT_ASSESSED"
            if row["reference_status"] == "EXACT_MATCH" and cand_ref and lead_ref and cand_effect and lead_effect:
                if int(cand_ref["POS"]) != bp or int(cand_effect["CHR"]) != chr_id or int(cand_effect["BP"]) != bp:
                    raise ValueError(f"Candidate coordinate mismatch across reference and joined inputs: {pair}/{snp}")
                lead_sign = allele_sign(lead_ref["A1"], lead_ref["A2"], lead_effect["A1"], lead_effect["A2"])
                cand_sign = allele_sign(cand_ref["A1"], cand_ref["A2"], cand_effect["A1"], cand_effect["A2"])
                if cand_sign is not None and lead_sign is not None:
                    orientation = "PASS_REFERENCE_A1_ALIGNED"
                    beta_sleep = float(cand_effect["BETA1"]) * cand_sign
                    beta_disorder = float(cand_effect["BETA2"]) * cand_sign
                    beta_sleep_lead = float(lead_effect["BETA1"]) * lead_sign
                    beta_disorder_lead = float(lead_effect["BETA2"]) * lead_sign
                    raw_r = edges.get((pair, snp, lead))
                    if raw_r is None:
                        raise ValueError(f"Candidate-to-lead LD edge missing: {pair}/{snp}/{lead}")
                    status = "OUT_OF_RANGE_RAW_LAVA_R" if abs(raw_r) > 1 else "PASS_VALID_SIGNED_R"
                else:
                    orientation = "ALLELES_NOT_UNAMBIGUOUSLY_ALIGNABLE"
            if beta_sleep is None or beta_disorder is None:
                direction = "NOT_ASSESSED"
            elif beta_sleep == 0 or beta_disorder == 0:
                direction = "ZERO_EFFECT_UNDEFINED"
            else:
                direction = "CONCORDANT" if (beta_sleep > 0) == (beta_disorder > 0) else "OPPOSING"
            valid = raw_r is not None and abs(raw_r) <= 1 and orientation == "PASS_REFERENCE_A1_ALIGNED"
            output_rows.append({
                "pair_id": pair, "SNP": row["SNP"], "CHR": chr_id, "BP": bp,
                "candidate_status": row["candidate_status"], "lead_SNP": row["lead_SNP"], "locus_id": row["locus_id"],
                "beta_sleep_reference_A1": f"{beta_sleep:.12g}" if beta_sleep is not None else "NA",
                "beta_disorder_reference_A1": f"{beta_disorder:.12g}" if beta_disorder is not None else "NA",
                "direction_relation": direction,
                "reference_A1": cand_ref["A1"] if cand_ref else "NA", "reference_A2": cand_ref["A2"] if cand_ref else "NA",
                "lead_reference_A1": lead_ref["A1"] if lead_ref else "NA", "lead_reference_A2": lead_ref["A2"] if lead_ref else "NA",
                "effect_orientation_status": orientation,
                "lava_r_raw": f"{raw_r:.12g}" if raw_r is not None else "NA",
                "signed_r_reference_A1": f"{raw_r:.12g}" if valid else "NA",
                "r2_reference": f"{raw_r * raw_r:.12g}" if valid else "NA",
                "ld_qc_status": status,
                "sleep_effect_vs_lead_ld": ld_direction_consistency(beta_sleep, beta_sleep_lead, raw_r) if valid else "NOT_ASSESSED",
                "disorder_effect_vs_lead_ld": ld_direction_consistency(beta_disorder, beta_disorder_lead, raw_r) if valid else "NOT_ASSESSED",
            })
        output_rows.sort(key=lambda r: (r["pair_id"], r["CHR"], r["BP"], r["SNP"]))
        atomic_tsv(output, FIELDS, output_rows)
        pair_counts = {}
        for pair in PAIRS:
            rows = [r for r in output_rows if r["pair_id"] == pair]
            pair_counts[pair] = {
                "candidate_rows": len(rows),
                "valid_signed_ld_rows": sum(r["ld_qc_status"] == "PASS_VALID_SIGNED_R" for r in rows),
                "out_of_range_raw_ld_rows": sum(r["ld_qc_status"] == "OUT_OF_RANGE_RAW_LAVA_R" for r in rows),
                "unmatched_or_unalignable_rows": sum(r["effect_orientation_status"] != "PASS_REFERENCE_A1_ALIGNED" for r in rows),
                "both_effects_consistent_with_ld_sign": sum(r["sleep_effect_vs_lead_ld"] == r["disorder_effect_vs_lead_ld"] == "CONSISTENT_WITH_LD_SIGN" for r in rows),
                "at_least_one_effect_opposite_to_ld_sign": sum("OPPOSITE_TO_LD_SIGN" in (r["sleep_effect_vs_lead_ld"], r["disorder_effect_vs_lead_ld"]) for r in rows),
            }
        atomic_tsv(summary_path, ["pair_id", "candidate_rows", "valid_signed_ld_rows", "out_of_range_raw_ld_rows",
                                 "unmatched_or_unalignable_rows", "both_effects_consistent_with_ld_sign",
                                 "at_least_one_effect_opposite_to_ld_sign", "interpretation"],
                   [{"pair_id": pair, **counts, "interpretation": "Post-selection descriptive LD/effect-sign comparison; partial PLACO family only; not causal or independently replicated"}
                    for pair, counts in pair_counts.items()])
        reference_files = {}
        for path, digest in clump_provenance["sources"]["reference_files"].items():
            reference_files[path] = digest
        provenance = {
            "schema_version": 1,
            "analysis_id": "brain6-placo-partial-signed-effect-ld-v1",
            "status": "PASS_PARTIAL_FAMILY_DESCRIPTIVE_WITH_LD_RANGE_EXCEPTIONS",
            "interpretation": ("Post-selection comparison of each PLACO candidate's harmonized sleep/disorder beta signs, "
                "both oriented to the official LAVA reference A1 allele, against signed reference LD to the range-filtered clump lead. "
                "Raw out-of-range LD is retained and excluded without clipping. This partial four-pair diagnostic is not an "
                "independent-locus test, fine-mapping result, replication, or causal claim; protected Track B is unavailable."),
            "scope": {"pairs": list(PAIRS), "candidate_rows": len(output_rows), "complete_five_track_family": False,
                      "pair_counts": pair_counts},
            "methods": {"effect_alignment": "Pair-joined BETA1/BETA2 aligned to the sleep-trait A1; both effects sign-flipped together as needed to orient to each reference A1",
                        "LD": "Pinned LAVA 0.1.5 read.ld from official UKB LAVA v1.1; candidate-to-range-filtered-lead signed r",
                        "invalid_ld": "Preserve raw r outside [-1,1]; omit signed correlation and direction comparison; do not clip",
                        "effect_vs_ld": "Descriptive indicator compares sign(beta_candidate * beta_lead) with sign(r); no statistical test"},
            "sources": {"candidate_table": {"path": str(CANDIDATE_PATH.relative_to(ROOT)), "sha256": sha256(CANDIDATE_PATH)},
                        "candidate_provenance": {"path": str(CLUMP_PROVENANCE.relative_to(ROOT)), "sha256": sha256(CLUMP_PROVENANCE)},
                        "joined_pair_inputs": pair_source_provenance,
                        "pair_receipt_files": {pair: str(PAIR_JOIN_ROOT / f"pair_{pair}/receipt.json") for pair in PAIRS},
                        "reference_provenance": {"path": str(REFERENCE_PROVENANCE), "sha256": sha256(REFERENCE_PROVENANCE)},
                        "reference_files": reference_files,
                        "extractor_script": {"path": str(R_EXTRACTOR.relative_to(ROOT)), "sha256": sha256(R_EXTRACTOR)},
                        "builder_script": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256(Path(__file__).resolve())},
                        "rscript_binary": {"path": str(RSCRIPT), "sha256": sha256(RSCRIPT)}},
            "outputs": [{"path": str(output.relative_to(ROOT)), "sha256": sha256(output), "rows": len(output_rows)},
                        {"path": str(summary_path.relative_to(ROOT)), "sha256": sha256(summary_path), "rows": len(pair_counts)}],
            "reference_sample_size": (ld_dir / "reference_sample_size.txt").read_text().strip(),
        }
        atomic_json(provenance_path, provenance)
        provenance["provenance_sha256"] = sha256(provenance_path)
        return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.provenance, args.summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
