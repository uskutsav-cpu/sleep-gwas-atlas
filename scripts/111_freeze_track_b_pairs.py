#!/usr/bin/env python3
"""Build the Track B brain-pair candidate table and freeze all three roles."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from datetime import datetime, timezone
from pathlib import Path


BRAIN_DOMAINS = {"neuro", "psychiatric"}
EXPECTED_BRAIN_TRAITS = {"alz", "parkinson", "mdd", "scz", "bipolar", "adhd"}
EXPECTED_PAIR_B = ("insomnia", "adhd")
CHECKPOINT = Path("results/track_b/00_repository_checkpoint.json")
SELECTION = Path("results/track_b/01_brain_pair_selection.tsv")
RATIONALE = Path("results/track_b/PAIR_B_SELECTION_RATIONALE.md")
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
LOCK = Path("results/track_b/pair_manifest.lock.json")

INPUTS = [
    Path("results/track_b/00_repository_checkpoint.json"),
    Path("config/analysis_panel.tsv"),
    Path("results/analysis/phase1_master_analysis.tsv"),
    Path("results/analysis/literature_novelty_audit.tsv"),
    Path("results/analysis/published_rg_replication.tsv"),
    Path("results/tables/h2_summary.tsv"),
    Path("results/tables/pleiotropy_input_readiness.tsv"),
    Path("results/tables/interpretation_source_readiness.tsv"),
    Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
]

SELECTION_FIELDS = [
    "sleep_trait", "brain_trait", "rg", "SE", "P", "FDR", "absolute_rg",
    "discovery_QC", "replication_status", "replication_rg", "replication_SE",
    "replication_P", "phenotype_match", "heterogeneity", "direct_prior_rg_status",
    "novelty_status", "dense_sleep_GWAS_available", "dense_brain_GWAS_available",
    "independent_GWAS_available", "human_brain_eQTL_available",
    "single_cell_reference_available", "scATAC_available", "brain_spatial_relevance",
    "major_caveats",
]

MANIFEST_FIELDS = [
    "pair_id", "sleep_trait", "external_trait", "exact_GWAS", "source", "publication",
    "ancestry", "build", "sample_size", "case_control_counts", "discovery_rg", "SE",
    "P", "FDR", "cross_trait_intercept", "SNP_overlap", "h2",
    "replication_availability", "dense_data_readiness", "primary_scientific_role",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def text_sha256(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unique(rows: list[dict[str, str]], key) -> dict[object, dict[str, str]]:
    result: dict[object, dict[str, str]] = {}
    for row in rows:
        identity = key(row)
        if identity in result:
            fail(f"duplicate canonical row: {identity}")
        result[identity] = row
    return result


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row.get(field, "NA") for field in fields} for row in rows)
    return buffer.getvalue()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def as_float(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def evidence_rank(row: dict[str, str]) -> tuple[int, int, int, float]:
    phenotype = row.get("external_phenotype_match", "")
    independence = row.get("comparison_independence", "")
    exact = int("EXACT" in phenotype)
    independence_rank = {
        "INDEPENDENT": 0,
        "PARTIALLY_INDEPENDENT_OR_OVERLAPPING": 1,
        "NOT_FULLY_INDEPENDENT": 2,
    }.get(independence, 3)
    concordant = int(row.get("direction_concordance") == "CONCORDANT")
    p = as_float(row.get("published_rg_p", ""))
    return (independence_rank, -exact, -concordant, p if p is not None else 1.0)


def replication_class(evidence: dict[str, str] | None) -> str:
    if evidence is None:
        return "NO_DIRECT_RG_REPLICATION_IDENTIFIED"
    independence = evidence["comparison_independence"]
    direction = evidence["direction_concordance"]
    if independence == "INDEPENDENT" and direction == "CONCORDANT":
        return "STRONG_REPLICATION"
    if independence == "PARTIALLY_INDEPENDENT_OR_OVERLAPPING" and direction == "CONCORDANT":
        return "DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE"
    if direction == "CONCORDANT":
        return "CONCORDANT_PRIOR_RG_NOT_INDEPENDENT"
    return "DISCORDANT_PRIOR_RG_NOT_INDEPENDENT"


def heterogeneity(discovery: dict[str, str], evidence: dict[str, str] | None) -> str:
    if evidence is None:
        return "NOT_COMPUTABLE_NO_DIRECT_RG"
    values = [
        as_float(discovery["rg"]), as_float(discovery["se"]),
        as_float(evidence["published_rg"]), as_float(evidence["published_rg_se"]),
    ]
    if any(value is None for value in values):
        return "NOT_COMPUTABLE_MISSING_SE"
    atlas_rg, atlas_se, prior_rg, prior_se = values
    assert atlas_rg is not None and atlas_se is not None and prior_rg is not None and prior_se is not None
    z = (atlas_rg - prior_rg) / math.sqrt(atlas_se**2 + prior_se**2)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    overlap = evidence["comparison_independence"]
    return f"z={z:.6g};p={p:.6g};assumes_independence_but_dataset_relation={overlap}"


def validate_checkpoint_inputs() -> None:
    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    for path_text, record in checkpoint["input_hashes"].items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != record["sha256"]:
            fail(f"canonical input drifted after repository checkpoint: {path}")
    if checkpoint["hard_constraints"]["track_b_local_results_accessed_before_checkpoint"] is not False:
        fail("repository checkpoint does not certify result-free Track B selection")


def make_selection() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    panel_rows = read_tsv(Path("config/analysis_panel.tsv"))
    panel = unique(panel_rows, lambda row: row["trait_id"])
    brain_traits = {row["trait_id"] for row in panel_rows if row["domain"] in BRAIN_DOMAINS}
    if brain_traits != EXPECTED_BRAIN_TRAITS:
        fail(f"eligible locked-core brain traits drifted: {sorted(brain_traits)}")
    sleep_traits = [row["trait_id"] for row in panel_rows if row["domain"] == "sleep"]
    master = unique(
        read_tsv(Path("results/analysis/phase1_master_analysis.tsv")),
        lambda row: (row["sleep_trait"], row["external_trait"]),
    )
    audit = unique(
        read_tsv(Path("results/analysis/literature_novelty_audit.tsv")),
        lambda row: (row["sleep_trait"], row["external_trait"]),
    )
    h2 = unique(read_tsv(Path("results/tables/h2_summary.tsv")), lambda row: row["trait"])
    dense = unique(read_tsv(Path("results/tables/pleiotropy_input_readiness.tsv")), lambda row: row["trait_id"])
    evidence_by_pair: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in read_tsv(Path("results/analysis/published_rg_replication.tsv")):
        evidence_by_pair.setdefault((row["sleep_trait"], row["external_trait"]), []).append(row)

    selection_rows: list[dict[str, object]] = []
    for sleep_trait in sleep_traits:
        for brain_trait in sorted(brain_traits):
            identity = (sleep_trait, brain_trait)
            if identity not in master:
                fail(f"candidate missing from canonical Phase-1 analysis: {identity}")
            discovery = master[identity]
            prior = audit.get(identity)
            evidence_rows = sorted(evidence_by_pair.get(identity, []), key=evidence_rank)
            evidence = evidence_rows[0] if evidence_rows else None
            discovery_qc = (
                "PASS" if discovery["primary_or_sensitivity"] == "PRIMARY_PHASE1"
                and discovery["interpretation_status"] == "PRIMARY"
                and h2[sleep_trait]["verdict"] == "PASS"
                and h2[brain_trait]["verdict"] == "PASS"
                else "FAIL"
            )
            sleep_dense = dense[sleep_trait]["input_status"] == "READY_FULL_SUMSTATS"
            brain_dense = dense[brain_trait]["input_status"] == "READY_FULL_SUMSTATS"
            repl_class = replication_class(evidence)
            caveats = []
            if float(discovery["fdr"]) >= 0.05:
                caveats.append("GLOBAL_FDR_NOT_SIGNIFICANT")
            if not brain_dense:
                caveats.append("BRAIN_GWAS_HAPMAP3_ONLY_FOR_DOWNSTREAM")
            if repl_class != "STRONG_REPLICATION":
                caveats.append("NO_FULLY_INDEPENDENT_EXACT_REPLICATION_AT_SELECTION")
            if prior is not None and prior["direct_rg_found"] == "True":
                caveats.append("DIRECT_RG_ALREADY_REPORTED")
            if prior is None:
                caveats.append("NO_PAIR_SPECIFIC_LITERATURE_AUDIT_NONDISCOVERY")
            if evidence and evidence["comparison_independence"] != "INDEPENDENT":
                caveats.append("PRIOR_RG_DATASET_OVERLAP_OR_PARTIAL_INDEPENDENCE")
            caveats.append("CELL_AND_QTL_RESOURCES_ARE_GENERIC_NOT_LOCUS_VALIDATION")
            selection_rows.append({
                "sleep_trait": sleep_trait,
                "brain_trait": brain_trait,
                "rg": discovery["rg"],
                "SE": discovery["se"],
                "P": discovery["p"],
                "FDR": discovery["fdr"],
                "absolute_rg": discovery["abs_rg"],
                "discovery_QC": discovery_qc,
                "replication_status": repl_class,
                "replication_rg": evidence["published_rg"] if evidence else "NA",
                "replication_SE": evidence["published_rg_se"] if evidence else "NA",
                "replication_P": evidence["published_rg_p"] if evidence else "NA",
                "phenotype_match": (
                    evidence["external_phenotype_match"] if evidence
                    else prior["phenotype_definition_match"] if prior is not None
                    else "NOT_AUDITED_NONDISCOVERY"
                ),
                "heterogeneity": heterogeneity(discovery, evidence),
                "direct_prior_rg_status": (
                    "DIRECT_PRIOR_RG_FOUND" if prior is not None and prior["direct_rg_found"] == "True"
                    else "NO_DIRECT_PRIOR_RG_FOUND" if prior is not None
                    else "NOT_AUDITED_NONDISCOVERY"
                ),
                "novelty_status": (
                    prior["novelty_classification"] if prior is not None
                    else "NOT_APPLICABLE_GLOBAL_FDR_NOT_SIGNIFICANT"
                ),
                "dense_sleep_GWAS_available": "YES" if sleep_dense else "NO",
                "dense_brain_GWAS_available": "YES" if brain_dense else "NO",
                "independent_GWAS_available": (
                    prior["replication_feasibility"] if prior is not None
                    else "NOT_ASSESSED_NONDISCOVERY"
                ),
                "human_brain_eQTL_available": "SEARCH_INFRASTRUCTURE_ONLY;NO_PAIR_SPECIFIC_QTL_LOCK",
                "single_cell_reference_available": "YES_GENERIC_FUMA_SCRNA;PAIR_SPECIFIC_ANALYSIS_NOT_RUN",
                "scATAC_available": "YES_GENERIC_CATLAS_ADULT;PAIR_SPECIFIC_ANALYSIS_NOT_RUN",
                "brain_spatial_relevance": "BIOLOGICALLY_RELEVANT;NO_AUTHORIZED_SPATIAL_DATASET_FOUND",
                "major_caveats": ";".join(caveats),
            })

    if len(selection_rows) != 72:
        fail(f"expected complete 12 x 6 brain-pair table, observed {len(selection_rows)}")
    shortlist = [
        row for row in selection_rows
        if row["discovery_QC"] == "PASS"
        and float(row["FDR"]) < 0.05
        and row["dense_sleep_GWAS_available"] == "YES"
        and row["dense_brain_GWAS_available"] == "YES"
        and "EXACT" in str(row["phenotype_match"])
        and row["replication_status"] in {"STRONG_REPLICATION", "DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE"}
    ]
    shortlist.sort(key=lambda row: (-float(row["absolute_rg"]), float(row["P"]), row["sleep_trait"], row["brain_trait"]))
    if not shortlist or (shortlist[0]["sleep_trait"], shortlist[0]["brain_trait"]) != EXPECTED_PAIR_B:
        observed = (shortlist[0]["sleep_trait"], shortlist[0]["brain_trait"]) if shortlist else None
        fail(f"frozen lexicographic selection no longer yields Pair B {EXPECTED_PAIR_B}: {observed}")
    return selection_rows, shortlist


def make_pair_manifest(selection_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    panel = unique(read_tsv(Path("config/analysis_panel.tsv")), lambda row: row["trait_id"])
    dense = unique(read_tsv(Path("results/tables/pleiotropy_input_readiness.tsv")), lambda row: row["trait_id"])
    master = unique(
        read_tsv(Path("results/analysis/phase1_master_analysis.tsv")),
        lambda row: (row["sleep_trait"], row["external_trait"]),
    )
    selection = {(str(row["sleep_trait"]), str(row["brain_trait"])): row for row in selection_rows}
    roles = [
        ("A", "snoring", "parental_lifespan", "PRIMARY_DISCOVERY_PAIR_A_FIXED_A_PRIORI"),
        ("B", EXPECTED_PAIR_B[0], EXPECTED_PAIR_B[1], "PRIMARY_DISCOVERY_PAIR_B_FROZEN_BY_SELECTION_RULE"),
        ("CONTROL", "insomnia", "frailty", "POSITIVE_CONTROL_KNOWN_INSOMNIA_FRAILTY_ARCHITECTURE"),
    ]
    output: list[dict[str, object]] = []
    for pair_id, sleep_trait, external_trait, role in roles:
        row = master[(sleep_trait, external_trait)]
        sleep = panel[sleep_trait]
        external = panel[external_trait]
        if pair_id == "A":
            replication = "NO_VALID_EXACT_NONOVERLAPPING_PUBLIC_DATASET_IN_PREEXISTING_AUDIT;FRESH_REAUDIT_REQUIRED"
        elif pair_id == "B":
            replication = str(selection[(sleep_trait, external_trait)]["replication_status"]) + ";FRESH_INDEPENDENCE_AUDIT_REQUIRED"
        else:
            replication = "DIRECT_PRIOR_RG_NOT_FULLY_INDEPENDENT;POSITIVE_CONTROL_NOT_NOVELTY_CLAIM"
        dense_ready = (
            "READY_FULL_SUMSTATS_BOTH"
            if dense[sleep_trait]["input_status"] == "READY_FULL_SUMSTATS"
            and dense[external_trait]["input_status"] == "READY_FULL_SUMSTATS"
            else "BLOCKED_FULL_SUMSTATS"
        )
        output.append({
            "pair_id": pair_id,
            "sleep_trait": sleep_trait,
            "external_trait": external_trait,
            "exact_GWAS": f"sleep={sleep['dataset_version']};external={external['dataset_version']}",
            "source": f"sleep={sleep['source_id']};external={external['source_id']}",
            "publication": f"sleep=PMID:{sleep['pmid']},DOI:{sleep['doi']};external=PMID:{external['pmid']},DOI:{external['doi']}",
            "ancestry": f"sleep={sleep['ancestry']};external={external['ancestry']}",
            "build": f"sleep={sleep['build']};external={external['build']}",
            "sample_size": f"sleep={sleep['n_total']};external={external['n_total']}",
            "case_control_counts": (
                f"sleep={sleep['ncase']}/{sleep['ncontrol']};external={external['ncase']}/{external['ncontrol']}"
            ),
            "discovery_rg": row["rg"], "SE": row["se"], "P": row["p"], "FDR": row["fdr"],
            "cross_trait_intercept": f"{row['gcov_intercept']} (SE {row['gcov_intercept_se']})",
            "SNP_overlap": row["SNP_overlap"],
            "h2": f"sleep={row['sleep_h2']} (SE {row['sleep_h2_se']});external={row['external_h2']} (SE {row['external_h2_se']})",
            "replication_availability": replication,
            "dense_data_readiness": dense_ready,
            "primary_scientific_role": role,
        })
    return output


def rationale_text(shortlist: list[dict[str, object]], selection_hash: str, manifest_hash: str, pair_b_hash: str) -> str:
    lines = [
        "# Track B Pair B selection rationale",
        "",
        "## Frozen decision",
        "",
        "Pair B is **insomnia ↔ ADHD**. This identity was frozen before any Track B local, pleiotropy, fine-mapping, colocalization, tissue, cell, or molecular-QTL result existed.",
        "",
        f"- Pair B identity SHA-256 (`B\\tinsomnia\\tadhd\\n`): `{pair_b_hash}`",
        f"- Complete 72-row candidate table SHA-256: `{selection_hash}`",
        f"- Three-role pair manifest SHA-256: `{manifest_hash}`",
        "",
        "## Eligible family",
        "",
        "The locked core contains 12 sleep traits and six eligible diagnosed brain/psychiatric disorders (Alzheimer's disease, Parkinson's disease, MDD, schizophrenia, bipolar disorder, and ADHD), yielding 72 candidates. The existing extension contains psychosocial traits but no additional diagnosed brain-disorder phenotype meeting this definition, so it contributes no eligible Pair B row.",
        "Non-significant candidates were retained in the complete table and explicitly marked `NOT_AUDITED_NONDISCOVERY`; the pre-existing literature audit was discovery-only and cannot support pair-specific novelty conclusions for those rows.",
        "",
        "## Prospective selection procedure",
        "",
        "Candidates were filtered lexicographically, not assigned an opaque score:",
        "",
        "1. clean primary Phase-1 status and PASS h² QC for both traits;",
        "2. locked core FDR < 0.05;",
        "3. dense signed GWAS available for both traits;",
        "4. exact-concept phenotype match in the existing evidence map;",
        "5. concordant direct rg evidence from a changed/expanded dataset with at least partial independence;",
        "6. among survivors, descending |rg|, then ascending P, then stable trait identifiers.",
        "",
        "Generic human single-cell and scATAC resources were recorded as feasibility evidence only. They were not treated as pair-specific cell/QTL validation, and no downstream result was inspected.",
        "",
        "## Ordered survivors",
        "",
        "| Rank | Sleep trait | Brain trait | rg | FDR | Existing replication class |",
        "|---:|---|---|---:|---:|---|",
    ]
    for rank, row in enumerate(shortlist, 1):
        lines.append(
            f"| {rank} | {row['sleep_trait']} | {row['brain_trait']} | {float(row['rg']):.4f} | "
            f"{float(row['FDR']):.4g} | {row['replication_status']} |"
        )
    lines.extend([
        "",
        "## Why Pair B won",
        "",
        "Insomnia–ADHD was the leading survivor because it combined clean h²/QC, an extremely strong locked global signal, dense inputs for both traits, exact-concept prior evidence, and the largest absolute rg among candidates meeting all preceding gates. This is not a novelty claim: direct rg has already been reported, and the prior comparison is only partially independent/overlapping.",
        "",
        "## Required next gate",
        "",
        "The existing evidence is not a fully independent replication. Pair B cannot advance beyond a global finding until a fresh cohort-overlap audit and, if feasible, standardized LDSC against a genuinely independent external ADHD GWAS are completed. A materially incompatible, well-powered opposite-direction result is a no-go for mechanistic escalation.",
        "",
        "## Non-selection caveats",
        "",
        "MDD candidates were not eligible for the dense-input gate because the local file is HapMap3-prefiltered. Parkinson candidates likewise lack dense downstream-ready input. Literature saturation was retained as a caveat rather than used post hoc to override the prospective gates. No candidate may replace Pair B because a later analysis is easier or more positive.",
        "",
    ])
    return "\n".join(lines)


def verify() -> None:
    for path in (SELECTION, RATIONALE, PAIR_MANIFEST, LOCK):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"missing frozen Track B pair artifact: {path}")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    if lock["pair_b"] != {"sleep_trait": EXPECTED_PAIR_B[0], "brain_trait": EXPECTED_PAIR_B[1]}:
        fail("Pair B identity differs from the frozen expected pair")
    if lock["pair_b_identity_sha256"] != text_sha256("B\tinsomnia\tadhd\n"):
        fail("Pair B identity hash is invalid")
    expected_outputs = {
        str(SELECTION): sha256(SELECTION),
        str(RATIONALE): sha256(RATIONALE),
        str(PAIR_MANIFEST): sha256(PAIR_MANIFEST),
    }
    if lock["output_sha256"] != expected_outputs:
        fail("a frozen Track B pair artifact changed")
    observed_inputs = {str(path): sha256(path) for path in INPUTS}
    if lock["input_sha256"] != observed_inputs:
        fail("a Pair B selection input changed")
    selection_rows = read_tsv(SELECTION)
    manifest_rows = read_tsv(PAIR_MANIFEST)
    if len(selection_rows) != 72 or [row["pair_id"] for row in manifest_rows] != ["A", "B", "CONTROL"]:
        fail("frozen candidate table or pair-manifest membership is invalid")
    pair_b = next(row for row in manifest_rows if row["pair_id"] == "B")
    if (pair_b["sleep_trait"], pair_b["external_trait"]) != EXPECTED_PAIR_B:
        fail("Pair B was replaced after freezing")
    if lock["downstream_results_accessed_before_pair_freeze"] is not False:
        fail("pair lock does not certify a result-free freeze")
    print(
        "TRACK_B_PAIR_MANIFEST_OK "
        f"pair_B={EXPECTED_PAIR_B[0]}__{EXPECTED_PAIR_B[1]} sha256={sha256(PAIR_MANIFEST)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify()
        return
    if any(path.exists() for path in (SELECTION, RATIONALE, PAIR_MANIFEST, LOCK)):
        fail("Track B pair artifacts already exist; use --verify instead of regenerating")
    validate_checkpoint_inputs()
    selection_rows, shortlist = make_selection()
    atomic_text(SELECTION, table_text(SELECTION_FIELDS, selection_rows))
    manifest_rows = make_pair_manifest(selection_rows)
    atomic_text(PAIR_MANIFEST, table_text(MANIFEST_FIELDS, manifest_rows))
    pair_b_identity_hash = text_sha256("B\tinsomnia\tadhd\n")
    rationale = rationale_text(shortlist, sha256(SELECTION), sha256(PAIR_MANIFEST), pair_b_identity_hash)
    atomic_text(RATIONALE, rationale)
    lock = {
        "schema_version": "1.0.0",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pair_b": {"sleep_trait": EXPECTED_PAIR_B[0], "brain_trait": EXPECTED_PAIR_B[1]},
        "pair_b_identity_sha256": pair_b_identity_hash,
        "candidate_family": "12 locked sleep traits x 6 locked core brain/psychiatric disorders",
        "candidate_count": len(selection_rows),
        "selection_rule": "PASS_QC -> core_FDR<0.05 -> dense_both -> exact_concept -> concordant_partial_or_full_independence -> descending_abs_rg -> ascending_P -> stable_ID",
        "thresholds_frozen_before_downstream": {"core_FDR": 0.05},
        "downstream_results_accessed_before_pair_freeze": False,
        "pair_replacement_policy": "FORBIDDEN; failures remain attached to the frozen pair",
        "input_sha256": {str(path): sha256(path) for path in INPUTS},
        "output_sha256": {
            str(SELECTION): sha256(SELECTION),
            str(RATIONALE): sha256(RATIONALE),
            str(PAIR_MANIFEST): sha256(PAIR_MANIFEST),
        },
    }
    atomic_text(LOCK, json.dumps(lock, indent=2, sort_keys=True) + "\n")
    print(
        "TRACK_B_PAIRS_FROZEN "
        f"pair_B={EXPECTED_PAIR_B[0]}__{EXPECTED_PAIR_B[1]} candidates={len(selection_rows)} "
        f"pair_b_sha256={pair_b_identity_hash}"
    )


if __name__ == "__main__":
    main()
