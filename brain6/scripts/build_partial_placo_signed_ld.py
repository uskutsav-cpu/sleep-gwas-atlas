"""Export signed LAVA-reference LD for the existing partial PLACO clumps.

The coefficient is oriented to each LAVA reference A1 allele. Summary-statistic
effect alleles are absent from the published PLACO variant output, so this
artifact reports LD only and does not claim aligned trait-effect signs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from clump_placo_candidates import extract_reference_ld

CANDIDATE_PATH = ROOT / "brain6/results/loci/placo_candidate_variants_partial.tsv"
LOCUS_PROVENANCE = ROOT / "brain6/results/loci/placo_candidate_loci_partial.provenance.json"
MASTER_PATH = ROOT / "brain6/results/placo/placo_master.tsv"
RULE_PATH = ROOT / "brain6/config/shared_locus_rule_v1.json"
FAMILY_PATH = ROOT / "brain6/config/placo_family_v3/family_lock.json"
RUN_PATH = ROOT / "brain6/manifests/placo_family_v3_run.json"
REFERENCE_PREFIX = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1")
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
OUTPUT_PATH = ROOT / "brain6/results/loci/placo_candidate_signed_ld_partial.tsv"
PROVENANCE_PATH = ROOT / "brain6/results/loci/placo_candidate_signed_ld_partial.provenance.json"
FIELDS = ["pair_id", "SNP", "CHR", "BP", "candidate_status", "lead_SNP",
          "reference_status", "A1_REF", "A2_REF", "LEAD_A1_REF", "LEAD_A2_REF",
          "lava_r_raw", "signed_r_ref_A1", "r2_ref", "ld_qc_status",
          "r2_matches_existing_clump", "locus_id"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def atomic_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        with temp.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t",
                                    lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def make_signed_rows(candidates: list[dict[str, str]], reference_info: dict,
                     ld: dict[tuple[str, str], float]) -> list[dict[str, Any]]:
    """Join candidate-to-lead r and cross-check existing frozen clump r² values."""
    unique_ids: dict[str, set[int]] = {}
    for row in candidates:
        unique_ids.setdefault(row["SNP"].lower(), set()).add(int(row["CHR"]))
    if any(len(chromosomes) > 1 for chromosomes in unique_ids.values()):
        raise ValueError("An rsID maps to multiple chromosomes; LD output lacks chromosome labels")

    result = []
    for row in candidates:
        chrom, snp = int(row["CHR"]), row["SNP"].lower()
        lead = row["lead_SNP"].lower()
        ref = reference_info.get((chrom, snp))
        lead_ref = reference_info.get((chrom, lead)) if lead not in {"", "na"} else None
        raw_r: float | None = None
        signed_r: float | None = None
        qc_status = "NO_EXACT_REFERENCE_OR_LEAD"
        r2_check = "NA"
        if ref and lead_ref:
            raw_r = float(ld.get(tuple(sorted((snp, lead))), 1.0 if snp == lead else math.nan))
            if not math.isfinite(raw_r):
                raise ValueError(f"Signed LD is absent for exact-reference pair {snp}/{lead}")
            if snp == lead and abs(raw_r - 1.0) > 1e-8:
                raise ValueError(f"Lead self-correlation is not one: {snp}")
            old_r2 = float(row["r2_to_lead"])
            if abs(raw_r * raw_r - old_r2) > 1e-6:
                raise ValueError(f"Signed-LD square differs from frozen clump r² for {snp}/{lead}")
            if abs(raw_r) > 1:
                qc_status = "OUT_OF_RANGE_RAW_LAVA_R"
                r2_check = "RAW_R_SQUARED_MATCHES_CLUMP_ONLY"
            else:
                signed_r = raw_r
                qc_status = "PASS_VALID_SIGNED_R"
                r2_check = "TRUE"
        result.append({
            "pair_id": row["pair_id"], "SNP": row["SNP"], "CHR": chrom, "BP": row["BP"],
            "candidate_status": row["candidate_status"], "lead_SNP": row["lead_SNP"],
            "reference_status": row["reference_status"],
            "A1_REF": ref["A1"] if ref else "NA", "A2_REF": ref["A2"] if ref else "NA",
            "LEAD_A1_REF": lead_ref["A1"] if lead_ref else "NA",
            "LEAD_A2_REF": lead_ref["A2"] if lead_ref else "NA",
            "lava_r_raw": f"{raw_r:.12g}" if raw_r is not None else "NA",
            "signed_r_ref_A1": f"{signed_r:.12g}" if signed_r is not None else "NA",
            "r2_ref": f"{signed_r * signed_r:.12g}" if signed_r is not None else "NA",
            "ld_qc_status": qc_status,
            "r2_matches_existing_clump": r2_check,
            "locus_id": row["locus_id"],
        })
    return result


def build(output: Path = OUTPUT_PATH, provenance_path: Path = PROVENANCE_PATH,
          reference_prefix: Path = REFERENCE_PREFIX, rscript: Path = RSCRIPT) -> dict:
    candidates = read_tsv(CANDIDATE_PATH)
    if len(candidates) != 2297:
        raise ValueError(f"Expected the preserved 2,297-row partial candidate set; got {len(candidates)}")
    locus_provenance = json.loads(LOCUS_PROVENANCE.read_text(encoding="utf-8"))
    candidate_record = next(x for x in locus_provenance["outputs"]
                            if x["path"] == str(CANDIDATE_PATH.relative_to(ROOT)))
    if sha256(CANDIDATE_PATH) != candidate_record["sha256"]:
        raise ValueError("Partial candidate source differs from its clumping provenance")
    master_rows = read_tsv(MASTER_PATH)
    master = {r["pair_id"]: r for r in master_rows}
    if len(master) != len(master_rows) or {r["pair_id"] for r in candidates} != set(master):
        raise ValueError("Partial candidates and PLACO master do not cover the same admitted pairs")
    for pair, record in master.items():
        source = ROOT / record["output_path"]
        if sha256(source) != record["output_sha256"]:
            raise ValueError(f"Pair variant output changed since PLACO publication: {pair}")

    unique: dict[tuple[int, str], dict[str, Any]] = {}
    for row in candidates:
        key = (int(row["CHR"]), row["SNP"].lower())
        unique[key] = {"CHR": key[0], "SNP": key[1], "BP": int(row["BP"])}
        lead = row["lead_SNP"].lower()
        if lead not in {"", "na"}:
            unique.setdefault((key[0], lead), {"CHR": key[0], "SNP": lead, "BP": None})
    reference_info, ld, ref_n = extract_reference_ld(list(unique.values()), reference_prefix, rscript)
    rows = make_signed_rows(candidates, reference_info, ld)
    atomic_tsv(output, rows)
    reference_provenance = reference_prefix.parent / "reference.provenance.json"
    chromosomes = sorted({int(row["CHR"]) for row in candidates})
    reference_files = {}
    for chrom in chromosomes:
        for suffix in ("info", "bcor"):
            path = reference_prefix.parent / f"{reference_prefix.name}_chr{chrom}.{suffix}"
            reference_files[str(path)] = sha256(path)
    signed_count = sum(row["ld_qc_status"] == "PASS_VALID_SIGNED_R" for row in rows)
    out_of_range = sum(row["ld_qc_status"] == "OUT_OF_RANGE_RAW_LAVA_R" for row in rows)
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6-placo-partial-candidate-signed-ld-v1",
        "status": "PASS_SIGNED_LD_WITH_REFERENCE_RANGE_EXCEPTIONS",
        "interpretation": (
            "UKB/LAVA-reference signed LD between each exact-matched candidate and its frozen clump lead. "
            "r sign is oriented to the reference A1 alleles shown here. Four raw LAVA values are outside "
            "the mathematical correlation range [-1,1]; they are preserved in lava_r_raw and withheld from "
            "signed_r_ref_A1. PLACO output lacks summary-statistic effect alleles, so this does not align "
            "trait Z signs to LD alleles, establish independent signals, complete the five-track family, "
            "or support causality."
        ),
        "scope": {"pairs": sorted(master), "pair_count": len(master), "candidate_rows": len(rows),
                  "signed_ld_rows": signed_count, "out_of_range_rows": out_of_range,
                  "unmatched_rows": sum(row["ld_qc_status"] == "NO_EXACT_REFERENCE_OR_LEAD" for row in rows),
                  "complete_five_track_family": False},
        "method": {"operation": "LAVA 0.1.5 read.ld on official UKB v1.1 .bcor files",
                   "unit": "candidate variant to its frozen within-pair greedy-clump lead",
                   "r2_crosscheck_tolerance": 1e-6,
                   "correlation_range": [-1, 1],
                   "out_of_range_policy": "Preserve raw LAVA output; do not label it a correlation or clamp it",
                   "effect_allele_alignment": "NOT_AVAILABLE_IN_PUBLISHED_PLACO_VARIANT_OUTPUT"},
        "sources": {
            "candidate_table": {"path": str(CANDIDATE_PATH.relative_to(ROOT)), "sha256": sha256(CANDIDATE_PATH)},
            "candidate_locus_provenance": {"path": str(LOCUS_PROVENANCE.relative_to(ROOT)), "sha256": sha256(LOCUS_PROVENANCE)},
            "placo_master": {"path": str(MASTER_PATH.relative_to(ROOT)), "sha256": sha256(MASTER_PATH)},
            "placo_family_lock": {"path": str(FAMILY_PATH.relative_to(ROOT)), "sha256": sha256(FAMILY_PATH)},
            "clumping_rule": {"path": str(RULE_PATH.relative_to(ROOT)), "sha256": sha256(RULE_PATH)},
            "placo_run_manifest": {"path": str(RUN_PATH.relative_to(ROOT)), "sha256": sha256(RUN_PATH)},
            "reference_provenance": {"path": str(reference_provenance), "sha256": sha256(reference_provenance)},
            "extractor_script": {"path": "brain6/scripts/extract_placo_ld.R",
                                 "sha256": sha256(ROOT / "brain6/scripts/extract_placo_ld.R")},
            "builder_script": {"path": str(Path(__file__).resolve().relative_to(ROOT)),
                               "sha256": sha256(Path(__file__).resolve())},
            "rscript_binary": {"path": str(rscript), "sha256": sha256(rscript)},
            "pair_variant_outputs": {pair: row["output_sha256"] for pair, row in sorted(master.items())},
            "reference_files": reference_files,
        },
        "output": {"path": str(output.relative_to(ROOT)), "sha256": sha256(output), "rows": len(rows)},
        "reference_sample_size": ref_n,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance["provenance_sha256"] = sha256(provenance_path)
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-prefix", type=Path, default=REFERENCE_PREFIX)
    parser.add_argument("--rscript", type=Path, default=RSCRIPT)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE_PATH)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.provenance, args.reference_prefix, args.rscript),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
