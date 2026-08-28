#!/usr/bin/env python3
"""Apply the prespecified A/B/C prioritization rules to the extension family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def optional_by_pair(path: Path | None, status_field: str) -> dict[str, str]:
    if path is None:
        return {}
    rows = read_tsv(path)
    if not rows or "pair_id" not in rows[0] or status_field not in rows[0]:
        raise SystemExit(f"ERROR: follow-up table lacks pair_id/{status_field}: {path}")
    if len({row["pair_id"] for row in rows}) != len(rows):
        raise SystemExit(f"ERROR: duplicate pair_id in {path}")
    return {row["pair_id"]: row[status_field] for row in rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rg", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
    )
    parser.add_argument(
        "--extension-h2", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_trait_readiness.tsv"),
    )
    parser.add_argument("--core-h2", type=Path, default=Path("results/tables/h2_summary.tsv"))
    parser.add_argument(
        "--novelty-audit", type=Path,
        default=Path("discovery_extension/results/novelty/extension_novelty_audit.tsv"),
    )
    parser.add_argument("--replication", type=Path)
    parser.add_argument("--local-evidence", type=Path)
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/results/prioritization/novel_hit_priority.tsv"),
    )
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/prioritization.json"),
    )
    args = parser.parse_args()

    audit_validation = subprocess.run(
        [
            "python3", "discovery_extension/scripts/18_validate_pair_novelty_audit.py",
            "--rg", str(args.rg), "--audit", str(args.novelty_audit),
        ], check=False, text=True, capture_output=True,
    )
    if audit_validation.returncode:
        raise SystemExit(
            f"ERROR: pair-level novelty audit is not complete:\n{audit_validation.stdout}{audit_validation.stderr}"
        )

    rg_rows = read_tsv(args.rg)
    extension_h2 = {row["extension_trait_id"]: row for row in read_tsv(args.extension_h2)}
    core_h2 = {row["trait"]: row for row in read_tsv(args.core_h2)}
    audits = {row["pair_id"]: row for row in read_tsv(args.novelty_audit)}
    replication = optional_by_pair(args.replication, "replication_class")
    local = optional_by_pair(args.local_evidence, "local_support_status")
    output: list[dict[str, object]] = []
    for row in rg_rows:
        pair_id = f"{row['sleep_trait']}__{row['extension_trait_id']}"
        if row["extension_trait_id"] not in extension_h2:
            raise SystemExit(f"ERROR: extension h2 verdict missing for {pair_id}")
        if row["sleep_trait"] not in core_h2:
            raise SystemExit(f"ERROR: core sleep h2 verdict missing for {pair_id}")
        ext_pass = extension_h2[row["extension_trait_id"]]["primary_rg_eligibility"] == "PRIMARY_PASS"
        sleep_pass = core_h2[row["sleep_trait"]]["verdict"] == "PASS"
        analysis_clean = row["analysis_status"] == "PRIMARY_EXTENSION_RG_COMPLETE"
        fdr_pass = float(row["extension_fdr"]) < 0.05
        effect_pass = abs(float(row["rg"])) >= 0.15
        audit = audits.get(pair_id)
        novelty_strength = audit["novelty_strength"] if audit else "NOT_AUDITED_NOT_FDR_SIGNIFICANT"
        novelty_class = audit["novelty_class"] if audit else "NOT_AUDITED_NOT_FDR_SIGNIFICANT"
        novelty_pass = (
            novelty_strength in {"STRONG", "MODERATE"}
            and novelty_class in {"APPARENTLY_NOVEL", "NO_DIRECT_RG_FOUND"}
        )
        overlap = int(row["snp_overlap_valid_alleles"])
        overlap_strength = "HIGH_GE_900000" if overlap >= 900000 else ("MODERATE_500000_TO_899999" if overlap >= 500000 else "LOW_LT_500000")
        overlap_pass = overlap >= 500000
        replication_class = replication.get(
            pair_id, audit["independent_replication_class"] if audit else "NOT_YET_ATTEMPTED"
        )
        local_status = local.get(pair_id, "NOT_YET_ATTEMPTED")
        baseline = fdr_pass and effect_pass and ext_pass and sleep_pass and analysis_clean and novelty_pass and overlap_pass
        reinforced = replication_class == "REPLICATED" or local_status == "STRONG_LOCAL_SUPPORT"
        if baseline and reinforced:
            tier = "B"
            rationale = "Tier A discovery criteria plus independent replication or strong local support"
        elif baseline:
            tier = "A"
            rationale = "extension FDR<0.05, abs(rg)>=0.15, both h2 pass, clean QC, >=500000 valid-overlap SNPs, and Strong/Moderate APPARENTLY_NOVEL/NO_DIRECT_RG_FOUND evidence"
        else:
            tier = "C"
            failed = []
            if not fdr_pass: failed.append("extension_FDR_not_below_0.05")
            if not effect_pass: failed.append("abs_rg_below_0.15")
            if not ext_pass: failed.append("extension_h2_failed")
            if not sleep_pass: failed.append("sleep_h2_failed")
            if not analysis_clean: failed.append("analysis_status_not_clean")
            if not novelty_pass: failed.append("pair_novelty_class_or_strength_failed")
            if not overlap_pass: failed.append("valid_allele_SNP_overlap_below_500000")
            rationale = ";".join(failed)
        output.append({
            "pair_id": pair_id, "sleep_trait": row["sleep_trait"],
            "extension_trait_id": row["extension_trait_id"], "phenotype_name": row["phenotype_name"],
            "phenotype_domain": row["phenotype_domain"], "rg": row["rg"], "se": row["se"],
            "p": row["p"], "extension_fdr": row["extension_fdr"],
            "analysis_status": row["analysis_status"],
            "cross_trait_LDSC_intercept": row["cross_trait_LDSC_intercept"],
            "cross_trait_LDSC_intercept_se": row["cross_trait_LDSC_intercept_se"],
            "snp_overlap_valid_alleles": overlap, "snp_overlap_strength": overlap_strength,
            "extension_h2": extension_h2[row["extension_trait_id"]]["h2"],
            "extension_h2_z": extension_h2[row["extension_trait_id"]]["h2_z"],
            "extension_LDSC_intercept": extension_h2[row["extension_trait_id"]]["LDSC_intercept"],
            "sleep_h2": core_h2[row["sleep_trait"]]["h2"],
            "sleep_h2_z": core_h2[row["sleep_trait"]]["z"],
            "sleep_LDSC_intercept": core_h2[row["sleep_trait"]]["intercept"],
            "extension_h2_pass": str(ext_pass), "sleep_h2_pass": str(sleep_pass),
            "pair_novelty_class": novelty_class, "pair_novelty_strength": novelty_strength,
            "biological_plausibility": audit["biological_plausibility"] if audit else "NOT_ASSESSED",
            "connection_obviousness": audit["connection_obviousness"] if audit else "NOT_ASSESSED",
            "independent_replication_dataset_availability": audit["independent_replication_dataset_availability"] if audit else "NOT_ASSESSED",
            "dense_summary_statistics_available": audit["dense_summary_statistics_available"] if audit else "NOT_ASSESSED",
            "molecular_qtl_data_available": audit["molecular_qtl_data_available"] if audit else "NOT_ASSESSED",
            "replication_class": replication_class, "local_support_status": local_status,
            "priority_tier": tier, "priority_rationale": rationale,
            "claim_status": "REPLICATED_OR_LOCALLY_REINFORCED" if tier == "B" else ("DISCOVERY_ONLY" if tier == "A" else "NOT_PRIORITY_DISCOVERY"),
        })
    output.sort(key=lambda row: ({"B": 0, "A": 1, "C": 2}[str(row["priority_tier"])], float(row["extension_fdr"]), -abs(float(row["rg"])), str(row["pair_id"])))
    fields = list(output[0]) if output else []
    write_tsv(args.out, fields, output)
    counts = Counter(str(row["priority_tier"]) for row in output)
    provenance = {
        "schema_version": "1.0.0", "source_rg_sha256": sha256(args.rg),
        "extension_h2_sha256": sha256(args.extension_h2), "core_h2_sha256": sha256(args.core_h2),
        "novelty_audit_sha256": sha256(args.novelty_audit),
        "tier_rules": {
            "A": "FDR<0.05; abs(rg)>=0.15; both h2 pass; clean analysis status; >=500000 valid-overlap SNPs; Strong/Moderate APPARENTLY_NOVEL or NO_DIRECT_RG_FOUND pair",
            "B": "Tier A plus independent replication or strong local support",
            "C": "all other tested pairs, including nominal/suggestive, QC-sensitive, or novelty-ambiguous pairs",
        },
        "counts": dict(sorted(counts.items())), "output": str(args.out),
        "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"EXTENSION_PRIORITIZATION_OK pairs={len(output)} tiers={dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
