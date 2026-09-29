#!/usr/bin/env python3
"""Verify candidate-specific LDSC overlap matrices and the frozen LAVA gate."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN = Path("/Volumes/Extreme SSD/brain6-work/lava-power-optimized-sensitivity-v1/overlap_v1_retry")
ARCHIVE = RUN.parent / "reference/eur_w_ld_chr.tgz"
SCREEN = ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1"
CANONICAL = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
OUT = ROOT / "brain6/results/power_optimized_sensitivity_v1/overlap_integrity_audit_v1.json"
EXPECTED_MD5 = "76c1890c8cf22d99d05c6707cc8441b4"
PAIRS = ("scz", "bipolar", "parkinson")


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def log_intercept(path: Path, label: str) -> tuple[float, float | None]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if label == "Cross trait intercept":
        old = re.search(r"(?m)^\s*Cross[- ]trait intercept:\s*([-+0-9.eE]+)(?:\s*\(([-+0-9.eE]+)\))?", text)
        if old:
            return float(old.group(1)), float(old.group(2)) if old.group(2) else None
        section = re.search(r"(?s)Genetic Covariance\s*[-]+\s*(.*?)(?:\n\s*Genetic Correlation|\Z)", text)
        match = re.search(r"(?m)^\s*Intercept:\s*([-+0-9.eE]+)(?:\s*\(([-+0-9.eE]+)\))?", section.group(1)) if section else None
    else:
        match = re.search(rf"(?m)^\s*{re.escape(label)}:\s*([-+0-9.eE]+)(?:\s*\(([-+0-9.eE]+)\))?", text)
    if not match:
        raise ValueError(f"Could not read {label} from {path}")
    return float(match.group(1)), float(match.group(2)) if match.group(2) else None


def check_matrix(path: Path, pair: dict) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 3:
        raise ValueError(f"Malformed overlap matrix: {path}")
    h = lines[0].split()
    row1 = lines[1].split()
    row2 = lines[2].split()
    if h != [pair["trait1"], pair["trait2"]] or row1[0] != pair["trait1"] or row2[0] != pair["trait2"]:
        raise ValueError(f"Trait labels differ in overlap matrix: {path}")
    h1, c1 = map(float, row1[1:])
    c2, h2 = map(float, row2[1:])
    expected = tuple(float(pair[key]) for key in
                     ("sleep_intercept", "cross_trait_intercept", "disorder_intercept"))
    observed = (h1, c1, h2)
    if len(row1) != 3 or len(row2) != 3 or any(abs(x-y) > 5.1e-9 for x, y in zip(observed, expected)) or abs(c1-c2) > 1e-12:
        raise ValueError(f"Overlap matrix values differ from provenance: {path}")
    if h1 <= 0 or h2 <= 0 or c1*c1 >= h1*h2:
        raise ValueError(f"Overlap covariance matrix is not positive definite: {path}")


def validate() -> dict:
    provenance_path = RUN / "overlap.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("status") != "PASS_LDSC_INTERCEPTS_AND_MATRICES_CREATED":
        raise ValueError("LDSC overlap provenance is not a successful completion")
    if (not ARCHIVE.is_file() or digest(ARCHIVE, "md5") != EXPECTED_MD5 or
            digest(ARCHIVE) != provenance["reference_archive"]["archive_sha256"]):
        raise ValueError("LDSC reference archive checksum mismatch")
    ref_root = Path(provenance["reference_archive"]["reference_files_root"])
    hm3 = ref_root / "w_hm3.snplist"
    if digest(hm3) != provenance["reference_archive"]["hapmap3_snplist_sha256"]:
        raise ValueError("HapMap3 list checksum mismatch")
    if not (ref_root / "1.l2.ldscore.gz").is_file() or not (ref_root / "22.l2.M_5_50").is_file():
        raise ValueError("Reference panel is incomplete")
    candidate_receipt = ROOT / "brain6/results/power_optimized_sensitivity_v1/normalized/sleep_duration_continuous_dashti_2019/receipt.json"
    if digest(candidate_receipt) != provenance["candidate_normalization_receipt_sha256"]:
        raise ValueError("Candidate normalization receipt changed")
    if len(provenance.get("source_candidates", {})) != 4:
        raise ValueError("Expected one sleep and three disorder summary-statistic inputs")
    artifact_hashes = {"overlap_provenance_sha256": digest(provenance_path),
                       "reference_archive_sha256": digest(ARCHIVE),
                       "hapmap3_snplist_sha256": digest(hm3),
                       "candidate_normalization_receipt_sha256": digest(candidate_receipt)}
    for trait, record in provenance["source_candidates"].items():
        source = Path(record["source"])
        transformed = Path(record["ldsc_input"])
        if not source.is_file() or digest(source) != record["source_sha256"]:
            raise ValueError(f"Summary-statistic source changed for {trait}")
        if not transformed.is_file() or digest(transformed) != record["ldsc_input_sha256"]:
            raise ValueError(f"LDSC input changed for {trait}")
        munged = RUN / "munge" / f"{trait}.sumstats.gz"
        munge_log = RUN / "munge" / f"{trait}.log"
        if not munged.is_file() or munged.stat().st_size == 0 or not munge_log.is_file():
            raise ValueError(f"LDSC munging output incomplete for {trait}")
        artifact_hashes[f"munge_{trait}_sha256"] = digest(munged)
        artifact_hashes[f"munge_{trait}_log_sha256"] = digest(munge_log)
        h2_log = RUN / "ldsc" / f"h2_{trait}.log"
        if not h2_log.is_file():
            raise ValueError(f"Missing univariate LDSC intercept log for {trait}")
        artifact_hashes[f"h2_{trait}_log_sha256"] = digest(h2_log)
    if len(provenance.get("pairs", [])) != 3:
        raise ValueError("Expected three new duration-by-disorder overlap matrices")
    seen = set()
    for pair in provenance["pairs"]:
        trait = pair["trait2"]
        if trait not in PAIRS or trait in seen or pair["status"] != "PASS_LDSC_INTERCEPT_DERIVED":
            raise ValueError(f"Unexpected or duplicate overlap pair: {trait}")
        seen.add(trait)
        matrix = Path(pair["matrix_path"])
        if not matrix.is_file() or digest(matrix) != pair["matrix_sha256"]:
            raise ValueError(f"Matrix checksum mismatch for {trait}")
        check_matrix(matrix, pair)
        cross, cross_se = log_intercept(RUN / "ldsc" / f"rg_sleep_duration_continuous_dashti_2019__{trait}.log",
                                       "Cross trait intercept")
        h1, _ = log_intercept(RUN / "ldsc" / "h2_sleep_duration_continuous_dashti_2019.log", "Intercept")
        h2, _ = log_intercept(RUN / "ldsc" / f"h2_{trait}.log", "Intercept")
        if (abs(cross-float(pair["cross_trait_intercept"])) > 5.1e-9 or
                abs(float(cross_se)-float(pair["cross_trait_intercept_se"])) > 5.1e-9 or
                abs(h1-float(pair["sleep_intercept"])) > 5.1e-9 or
                abs(h2-float(pair["disorder_intercept"])) > 5.1e-9):
            raise ValueError(f"LDSC intercepts differ from recorded matrix for {trait}")
        artifact_hashes[f"rg_{trait}_log_sha256"] = digest(RUN / "ldsc" / f"rg_sleep_duration_continuous_dashti_2019__{trait}.log")
        artifact_hashes[f"matrix_{trait}_sha256"] = digest(matrix)
    if seen != set(PAIRS):
        raise ValueError("Duration-by-disorder overlap matrices are incomplete")

    screen_summary = json.loads((SCREEN / "summary.json").read_text(encoding="utf-8"))
    canonical_decision = json.loads((CANONICAL / "canonical_family_decision.json").read_text(encoding="utf-8"))
    canonical_longsleep = 1291
    canonical_total_not_run = int(canonical_decision["status_counts"]["NOT_RUN"])
    candidate_not_run = int(screen_summary["not_run"])
    family_not_run = canonical_total_not_run - canonical_longsleep + candidate_not_run
    family_cells = 17_465
    family_limit = math.floor(0.05 * family_cells)
    if (candidate_not_run != 876 or canonical_total_not_run != 3720 or
            family_not_run != 3305 or family_not_run <= family_limit):
        raise ValueError("Substitution-family missingness calculation changed")
    if (CANONICAL / "canonical_family_decision.json").is_file() is False:
        raise ValueError("Canonical decision receipt is missing")
    candidate_pairs = {"scz": ["266", "1719"], "bipolar": [], "parkinson": []}
    decision = {"family_cells": family_cells, "family_not_run": family_not_run,
                "family_not_run_fraction": family_not_run / family_cells,
                "maximum_not_run_cells": family_limit,
                "passes_frozen_5_percent_qc": family_not_run <= family_limit,
                "candidate_pair_gate_locus_intersections": candidate_pairs,
                "bivariate_lava_started": False,
                "reason_bivariate_lava_not_started": "SENSITIVITY_FAMILY_FAILS_FROZEN_UNIVARIATE_MISSINGNESS_GATE"}
    result = {"schema_version": 1,
              "analysis_id": "brain6_power_optimized_sensitivity_v1_overlap_integrity_audit",
              "status": "PASS_OVERLAP_AUDIT_FAMILY_FAILS_FROZEN_QC",
              "overlap_run_root": str(RUN),
              "overlap_provenance_sha256": digest(provenance_path),
              "artifact_sha256": artifact_hashes,
              "family_decision": decision,
              "canonical_run_id": canonical_decision["run_id"],
              "canonical_status": canonical_decision["overall_status"],
              "canonical_decision_sha256": digest(CANONICAL / "canonical_family_decision.json"),
              "interpretation": "Candidate-specific LDSC covariance inputs are valid and source-bound. The requested full local-h2 screen is complete, but the seven-trait substitution family exceeds the frozen untested-cell ceiling; no bivariate LAVA or downstream association result is authorized."}
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if OUT.exists():
        if OUT.read_text(encoding="utf-8") != serialized:
            raise FileExistsError(f"Existing audit differs; refusing to overwrite {OUT}")
    else:
        OUT.write_text(serialized, encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
