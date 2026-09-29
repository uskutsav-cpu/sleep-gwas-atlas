#!/usr/bin/env python3
"""Reproducible, read-only audit of frozen LAVA v3 and bounded rescue QC pilot."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
RUN = ROOT / "work/lava-canonical-v3-production" / RUN_ID
AGG = RUN / "results/canonical_family_results.tsv"
TRAITS = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
PAIRS = {
    "insomnia__adhd": ("insomnia", "adhd"),
    "insomnia__mdd": ("insomnia", "mdd"),
    "longsleep__scz": ("longsleep", "scz"),
    "longsleep__bipolar": ("longsleep", "bipolar"),
    "longsleep__parkinson": ("longsleep", "parkinson"),
}
PILOT_LOCI = ("1", "950", "1416", "2207")
REASONS = {
    "LOW_LOCAL_H2_UNDERPOWERED": "GENUINELY_UNINFORMATIVE_LOCAL_H2",
    "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS": "REFERENCE_MATCH_OR_SPARSE_LOCUS",
    "FEWER_THAN_MIN_K": "TOO_FEW_LAVA_COMPONENTS",
}
FIELDS = (
    "unit", "evidence_origin", "pair_id", "trait_id", "locus_id", "chromosome",
    "status", "reason", "cause_category", "n_snps", "n_components", "h2_obs",
    "h2_latent", "p_univariate", "strict_univariate_gate", "pair_gate_status",
    "bivariate_status", "pair_trait_1_status", "pair_trait_2_status",
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_new(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Immutable rescue artifact differs: {path}")
        return
    path.write_bytes(payload)


def json_bytes(x: dict) -> bytes:
    return (json.dumps(x, sort_keys=True, indent=2) + "\n").encode()


def inspect_sumstats(path: Path, expected_sha: str) -> dict:
    observed_sha = sha(path)
    if observed_sha != expected_sha:
        raise ValueError(f"Sumstats hash mismatch: {path}")
    count = 0
    seen = set()
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != ["SNP", "A1", "A2", "Z", "N"]:
            raise ValueError(f"Malformed sumstats schema: {path}")
        for r in reader:
            count += 1
            snp = r["SNP"].lower()
            if not snp or snp in seen:
                raise ValueError(f"Missing/duplicate SNP: {path}:{count}")
            seen.add(snp)
            a1, a2 = r["A1"].upper(), r["A2"].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or a1 == a2:
                raise ValueError(f"Malformed allele: {path}:{count}")
            if not math.isfinite(float(r["Z"])) or not math.isfinite(float(r["N"])) or float(r["N"]) <= 0:
                raise ValueError(f"Malformed Z/N: {path}:{count}")
    return {"path": str(path), "sha256": observed_sha, "rows": count, "schema": "SNP,A1,A2,Z,N", "field_qc": "PASS"}


def main() -> None:
    protocol = HERE / "lava_rescue_v1_protocol.md"
    decision_path = RUN / "canonical_family_decision.json"
    decision = load(decision_path)
    if decision["overall_status"] != "FAILED_QC_NOT_PROMOTED" or decision["run_id"] != RUN_ID:
        raise ValueError("Frozen canonical run identity/status changed")
    if sha(AGG) != decision["aggregate"]["sha256"]:
        raise ValueError("Frozen canonical aggregate hash changed")
    if len(rows(AGG)) != 17465:
        raise ValueError("Frozen canonical aggregate dimensions changed")
    provenance_path = ROOT / "brain6/results/power_optimized_sensitivity_v1/canonical_family_trait_power_audit_v1.provenance.json"
    provenance = load(provenance_path)
    input_hashes = {p: h for p, h in provenance["inputs_sha256"].items() if p.endswith(".harmonized.tsv.gz")}
    if len(input_hashes) != 7:
        raise ValueError("Missing frozen GWAS identity")
    strict = float(decision["multiple_testing"]["strict_p_threshold"])
    source = rows(AGG)
    by_cell = {}
    matrix = []
    status = Counter()
    reasons = Counter()
    trait_counts = defaultdict(Counter)
    for r in source:
        trait, locus = r["phen"], r["locus_id"]
        if trait not in TRAITS or (trait, locus) in by_cell:
            raise ValueError("Duplicate/unknown canonical cell")
        by_cell[trait, locus] = r
        st, reason = r["status"], r["reason"]
        status[st] += 1
        trait_counts[trait][st] += 1
        if reason:
            if reason not in REASONS:
                raise ValueError(f"Unclassified reason: {reason}")
            reasons[reason] += 1
            trait_counts[trait][reason] += 1
        passed = st == "TESTED" and float(r["p"]) < strict
        if passed:
            trait_counts[trait]["STRICT_GATE_PASS"] += 1
        matrix.append({
            "unit": "trait_locus", "evidence_origin": "OBSERVED_CANONICAL_V3",
            "pair_id": "", "trait_id": trait, "locus_id": locus,
            "chromosome": r["chromosome"], "status": st, "reason": reason,
            "cause_category": REASONS.get(reason, ""), "n_snps": r["n_snps"],
            "n_components": r["n_components"], "h2_obs": r["h2.obs"],
            "h2_latent": r["h2.latent"], "p_univariate": r["p"],
            "strict_univariate_gate": "PASS" if passed else "NOT_RUN" if st != "TESTED" else "FAIL",
            "pair_gate_status": "", "bivariate_status": "", "pair_trait_1_status": "", "pair_trait_2_status": "",
        })
    loci = {r["locus_id"] for r in source}
    if len(loci) != 2495 or len(by_cell) != 17465 or any(status[k] != v for k, v in decision["status_counts"].items()):
        raise ValueError("Canonical counts disagree with frozen decision")
    pair_counts = defaultdict(Counter)
    for pair, (first, second) in PAIRS.items():
        for locus in sorted(loci, key=int):
            a, b = by_cell[first, locus], by_cell[second, locus]
            a_pass = a["status"] == "TESTED" and float(a["p"]) < strict
            b_pass = b["status"] == "TESTED" and float(b["p"]) < strict
            gate = "BOTH_STRICT_UNIV_PASS" if a_pass and b_pass else "AT_LEAST_ONE_NOT_RUN" if (a["status"] != "TESTED" or b["status"] != "TESTED") else "STRICT_UNIV_GATE_MISS"
            pair_counts[pair][gate] += 1
            matrix.append({
                "unit": "pair_locus", "evidence_origin": "DERIVED_UNIVARIATE_ELIGIBILITY_ONLY",
                "pair_id": pair, "trait_id": "", "locus_id": locus,
                "chromosome": a["chromosome"], "status": "BIVARIATE_NOT_EXECUTED",
                "reason": "CANONICAL_FAMILY_FAILED_UNIVARIATE_QC", "cause_category": "FAMILY_GATE",
                "n_snps": "", "n_components": "", "h2_obs": "", "h2_latent": "",
                "p_univariate": "", "strict_univariate_gate": "", "pair_gate_status": gate,
                "bivariate_status": "NOT_EXECUTED", "pair_trait_1_status": a["status"],
                "pair_trait_2_status": b["status"],
            })
    if len(matrix) != 17465 + 12475:
        raise ValueError("Matrix incomplete")
    matrix_path = HERE / "lava_v3_failure_matrix.tsv"
    from io import StringIO
    out = StringIO(newline="")
    writer = csv.DictWriter(out, FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(matrix)
    write_new(matrix_path, out.getvalue().encode())

    manifest = {
        "analysis_id": "brain6-lava-rescue-v1", "schema_version": 1,
        "protocol_path": str(protocol.relative_to(ROOT)), "protocol_sha256": sha(protocol),
        "canonical_run_id": RUN_ID, "canonical_decision_sha256": sha(decision_path),
        "canonical_aggregate_sha256": sha(AGG), "canonical_family_status": decision["overall_status"],
        "harmonized_gwas_sha256": input_hashes,
        "reference_provenance_sha256": provenance["inputs_sha256"]["ref/lava/ukb_v1.1/reference.provenance.json"],
        "locus_definition_sha256": load(RUN / "family_lock.json")["locus_definition"]["sha256"],
        "overlap_manifest_sha256": load(RUN / "family_lock.json")["overlap_gate"]["pairwise_sample_overlap_file_manifest"]["manifest_sha256"],
        "pair_ids": list(PAIRS), "trait_ids": list(TRAITS), "pilot_loci": list(PILOT_LOCI),
        "matrix_path": str(matrix_path.relative_to(ROOT)), "matrix_sha256": sha(matrix_path),
        "matrix_rows": len(matrix), "trait_locus_rows": 17465, "pair_locus_rows": 12475,
        "status": "PROTOCOL_FROZEN_DIAGNOSTIC_PILOT_ONLY",
    }
    write_new(HERE / "lava_rescue_v1_manifest.json", json_bytes(manifest))

    pilot = []
    for locus in PILOT_LOCI:
        receipt_path = RUN / "receipts" / f"locus_{locus}.json"
        receipt = load(receipt_path)
        checks = {"output_sha256": sha(ROOT / receipt["output"]["path"]) == receipt["output"]["sha256"],
                  "input_info_sha256": sha(Path(receipt["input_info"]["path"])) == receipt["input_info"]["sha256"],
                  "reference_provenance_sha256": sha(Path(receipt["reference_provenance"]["path"])) == receipt["reference_provenance"]["sha256"],
                  "worker_config_sha256": sha(ROOT / receipt["worker_config_path"]) == receipt["config_sha256"]}
        source_rows = rows(ROOT / receipt["output"]["path"])
        if len(source_rows) != 7 or {x["phen"] for x in source_rows} != set(TRAITS):
            raise ValueError(f"Pilot locus output incomplete: {locus}")
        input_info = rows(Path(receipt["input_info"]["path"]))
        checks["input_info_schema"] = len(input_info) == 7 and {x["phenotype"] for x in input_info} == set(TRAITS)
        sumstats = {}
        for trait in TRAITS:
            item = receipt["sumstats"][trait]
            sumstats[trait] = inspect_sumstats(Path(item["path"]), item["sha256"])
            cell = next(x for x in source_rows if x["phen"] == trait)
            normalized = {k: cell.get(k, "") for k in ("phen", "locus_id", "status", "n_snps", "n_components", "h2.obs", "h2.latent", "p", "reason")}
            cell_hash = hashlib.sha256(json.dumps(normalized, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()
            checks[f"{trait}_cell_identity"] = receipt["cells"][trait] == cell_hash
            canonical_cell = by_cell[trait, locus]
            checks[f"{trait}_status_concordance"] = cell["status"] == canonical_cell["status"] and cell["reason"] == canonical_cell["reason"]
        if not all(checks.values()):
            raise ValueError(f"Pilot receipt/input check failed: {locus}: {checks}")
        pilot.append({
            "analysis_id": "brain6-lava-rescue-v1-diagnostic-pilot", "locus_id": locus,
            "protocol_sha256": sha(protocol), "canonical_receipt_path": str(receipt_path.relative_to(ROOT)),
            "canonical_receipt_sha256": sha(receipt_path), "command": "python3 brain6/results/lava_rescue_v1/build_rescue_audit.py",
            "environment": {"python": sys.version.split()[0], "platform": platform.platform()},
            "checks": checks, "sumstats": sumstats, "canonical_statuses": {x["phen"]: {"status": x["status"], "reason": x["reason"]} for x in source_rows},
            "exit_status": 0, "qc": "PASS_INPUT_RECEIPT_INTEGRITY_NO_CORRECTABLE_DEFECT_SHOWN",
        })
    for item in pilot:
        write_new(HERE / f"pilot_locus_{item['locus_id']}.receipt.json", json_bytes(item))
    result = {
        "analysis_id": "brain6-lava-rescue-v1", "schema_version": 1,
        "protocol_sha256": sha(protocol), "manifest_sha256": sha(HERE / "lava_rescue_v1_manifest.json"),
        "matrix_sha256": sha(matrix_path), "canonical_decision_sha256": sha(decision_path),
        "canonical_observed_status_counts": dict(status), "canonical_not_run_reasons": dict(reasons),
        "per_trait": {t: dict(c) for t, c in trait_counts.items()},
        "per_pair_derived_gate_status": {p: dict(c) for p, c in pair_counts.items()},
        "pilot_loci": list(PILOT_LOCI), "pilot_receipts": {x["locus_id"]: sha(HERE / f"pilot_locus_{x['locus_id']}.receipt.json") for x in pilot},
        "pilot_exit_status": 0, "correctable_large_scale_technical_cause_demonstrated": False,
        "full_rescue_executed": False, "bivariate_slots_executed": 0,
        "minimum_untested_if_all_min_k_fixed": reasons["LOW_LOCAL_H2_UNDERPOWERED"],
        "minimum_untested_if_both_sleep_traits_perfect": sum(trait_counts[t]["NOT_RUN"] for t in ("adhd", "bipolar", "mdd", "parkinson", "scz")),
        "frozen_maximum_untested": decision["maximum_allowed_untested_cells"],
        "family_complete": False, "bh_12475_complete": False,
        "classification": "FAIL_QC", "promotion_permitted": False,
        "reason": "No demonstrated technical defect can bridge 3564 low-local-h2 cells to the fixed 873-cell maximum; the diagnostic pilot verified input and receipt integrity at four prospectively selected loci.",
    }
    write_new(HERE / "lava_rescue_v1_qc.json", json_bytes(result))
    print(json.dumps({"matrix_rows": len(matrix), "status": dict(status), "reasons": dict(reasons), "pairs": {p: dict(c) for p, c in pair_counts.items()}, "classification": result["classification"]}, indent=2))


if __name__ == "__main__":
    main()
