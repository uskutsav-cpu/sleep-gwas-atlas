#!/usr/bin/env python3
"""Build a source-hashed interim sensitivity ledger for current Brain6 claims.

Untested sensitivities remain explicit; this is not a final paper-readiness pass.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/sensitivity/headline_claim_sensitivity.tsv"
FIELDS = [
    "claim_id", "manuscript_claim", "evidence_status", "claim_scope",
    "evidence_sources", "alternative_gwas", "independent_replication",
    "sample_overlap_sensitivity", "phenotype_definition_sensitivity",
    "original_396_family_significance", "LAVA_UKB_v1_1",
    "profile_leave_one_sleep_trait_out",
    "HDL_L_local_rg", "local_h2_thresholds", "alternate_LD_reference",
    "MHC_exclusion", "fine_mapping_stability", "coloc_prior_sensitivity",
    "alternate_gene_mapping", "objective_vs_subjective_sleep",
    "ancestry_sensitivity", "overall_sensitivity_status", "interpretation_limit",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def evidence(*paths: Path) -> str:
    return ";".join(f"{p.relative_to(ROOT)}@sha256:{digest(p)}" for p in paths)


def status_defaults(**overrides: str) -> dict[str, str]:
    row = {field: "NOT_RUN" for field in FIELDS}
    row.update({
        "alternative_gwas": "NOT_RUN",
        "independent_replication": "NOT_ASSESSED",
        "sample_overlap_sensitivity": "NOT_RUN",
        "phenotype_definition_sensitivity": "NOT_RUN",
        "original_396_family_significance": "NOT_APPLICABLE",
        "LAVA_UKB_v1_1": "IN_PROGRESS",
        "HDL_L_local_rg": "NOT_AVAILABLE",
        "local_h2_thresholds": "IN_PROGRESS_WITH_LAVA",
        "alternate_LD_reference": "NOT_RUN",
        "MHC_exclusion": "NOT_RUN",
        "fine_mapping_stability": "NOT_RUN",
        "coloc_prior_sensitivity": "NOT_RUN",
        "alternate_gene_mapping": "NOT_RUN",
        "objective_vs_subjective_sleep": "NOT_RUN",
        "profile_leave_one_sleep_trait_out": "NOT_APPLICABLE_TO_CLAIM",
        "ancestry_sensitivity": "EUR_DISCOVERY_ONLY; NOT_TESTED",
        "overall_sensitivity_status": "INCOMPLETE_NOT_A_SURVIVAL_CLAIM",
    })
    row.update(overrides)
    return row


def main() -> None:
    global_path = ROOT / "brain6/results/global/brain6_72_locked.tsv"
    profile_path = ROOT / "brain6/results/global/disorder_profile_similarity.tsv"
    replication_path = ROOT / "brain6/results/replication/replication_master.tsv"
    placo_path = ROOT / "brain6/results/placo/placo_master.tsv"
    locus_path = ROOT / "brain6/results/loci/placo_candidate_loci_partial.provenance.json"
    context_source_path = ROOT / "brain6/manifests/published_context_literature.json"
    context_table_path = ROOT / "brain6/results/supplement/table_S18_published_rg_context.tsv"
    profile_loo_path = ROOT / "brain6/results/global/disorder_profile_leave_one_sleep_trait_out.tsv"
    profile_loo_provenance_path = ROOT / "brain6/results/global/disorder_profile_leave_one_sleep_trait_out.provenance.json"
    profile_loo_plan_path = ROOT / "brain6/config/brain6_profile_leave_one_out_v1.json"
    mhc_plan_path = ROOT / "brain6/config/placo_candidate_mhc_exclusion_v1.json"
    mhc_path = ROOT / "brain6/results/sensitivity/placo_candidate_mhc_exclusion.tsv"
    mhc_provenance_path = ROOT / "brain6/results/sensitivity/placo_candidate_mhc_exclusion.provenance.json"
    mhc_script_path = ROOT / "brain6/scripts/analyze_placo_candidate_mhc_exclusion.py"
    lava_path = ROOT / "brain6/manifests/lava_family_v2_run.json"
    for path in (global_path, profile_path, replication_path, placo_path, locus_path,
                 context_source_path, context_table_path, lava_path, profile_loo_path,
                 profile_loo_provenance_path, profile_loo_plan_path,
                 mhc_plan_path, mhc_path, mhc_provenance_path, mhc_script_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    global_rows = read_tsv(global_path)
    profiles = read_tsv(profile_path)
    reps = read_tsv(replication_path)
    placo = read_tsv(placo_path)
    overlap_rows = read_tsv(ROOT / "brain6/manifests/cohort_overlap.tsv")
    overlap_by_pair = {r["pair"]: r for r in overlap_rows}
    context_rows = read_tsv(context_table_path)
    loci = json.loads(locus_path.read_text(encoding="utf-8"))
    profile_loo = read_tsv(profile_loo_path)
    profile_loo_provenance = json.loads(profile_loo_provenance_path.read_text(encoding="utf-8"))
    profile_loo_plan = json.loads(profile_loo_plan_path.read_text(encoding="utf-8"))
    mhc_plan = json.loads(mhc_plan_path.read_text(encoding="utf-8"))
    mhc_provenance = json.loads(mhc_provenance_path.read_text(encoding="utf-8"))
    profile_loo_plan_sha = profile_loo_plan_path.with_suffix(profile_loo_plan_path.suffix + ".sha256")
    if profile_loo_provenance.get("status") != "PASS":
        raise ValueError("Leave-one-sleep-trait-out profile sensitivity has not passed provenance checks")
    if (digest(profile_loo_plan_path) != profile_loo_plan_sha.read_text(encoding="utf-8").split()[0] or
            profile_loo_provenance.get("plan_sha256") != digest(profile_loo_plan_path) or
            profile_loo_provenance.get("global_map_sha256") != digest(global_path) or
            profile_loo_provenance.get("locked_family_config_sha256") != digest(ROOT / "brain6/config/brain6_locked_family.yaml") or
            profile_loo_provenance.get("script_sha256") != digest(ROOT / "brain6/scripts/build_profile_leave_one_out.py") or
            profile_loo_provenance.get("output_sha256") != digest(profile_loo_path)):
        raise ValueError("Leave-one-sleep-trait-out profile sensitivity provenance is stale")
    if (mhc_provenance.get("status") != "PASS_PARTIAL_FAMILY_DESCRIPTIVE_COUNTS" or
            mhc_provenance.get("analysis_id") != mhc_plan.get("analysis_id") or
            mhc_provenance.get("plan_sha256") != digest(mhc_plan_path) or
            mhc_provenance.get("script_path") != str(mhc_script_path.relative_to(ROOT)) or
            mhc_provenance.get("script_sha256") != digest(mhc_script_path) or
            mhc_provenance.get("output_sha256") != digest(mhc_path) or
            mhc_provenance.get("family_complete") is not False or
            mhc_provenance.get("protected_track_b_present") is not False):
        raise ValueError("Partial PLACO candidate MHC sensitivity provenance is missing or stale")
    mhc_rows = read_tsv(mhc_path)
    mhc_summary = mhc_rows[-1] if mhc_rows else {}
    mhc_pair_rows = mhc_rows[:-1]
    count_fields = ("candidate_variants_all", "candidate_variants_position_validated",
                    "candidate_variants_inside_mhc", "candidate_variants_outside_mhc",
                    "candidate_variants_unpositioned_or_mismatched", "candidate_intervals_all",
                    "candidate_intervals_overlapping_mhc", "candidate_intervals_wholly_outside_mhc")
    if (len(mhc_pair_rows) != len(mhc_plan["expected_pairs"]) or
            {row["pair_id"] for row in mhc_pair_rows} != set(mhc_plan["expected_pairs"]) or
            mhc_summary.get("pair_id") != "ALL_PARTIAL_PAIRS" or any(
                int(mhc_summary.get(field, -1)) != sum(int(row[field]) for row in mhc_pair_rows)
                for field in count_fields)):
        raise ValueError("Partial PLACO candidate MHC sensitivity table is incomplete or inconsistent")
    configured_profile_claims = {row["claim_id"] for row in profile_loo_plan["profiles"]}
    observed_profile_claims = {row["claim_id"] for row in profile_loo}
    expected_sleep_n = int(profile_loo_plan["method"]["expected_sleep_trait_count"])
    profile_loo_summary: dict[str, dict[str, float]] = {}
    for claim_id in configured_profile_claims:
        rows_for_claim = [row for row in profile_loo if row["claim_id"] == claim_id]
        if (len(rows_for_claim) != expected_sleep_n or
                len({row["omitted_sleep_trait"] for row in rows_for_claim}) != expected_sleep_n):
            raise ValueError(f"Incomplete leave-one-out profile rows for {claim_id}")
        minimum = min(float(row["leave_one_out_r"]) for row in rows_for_claim)
        maximum = max(float(row["leave_one_out_r"]) for row in rows_for_claim)
        if (abs(minimum - float(rows_for_claim[0]["leave_one_out_min_r"])) > 1e-10 or
                abs(maximum - float(rows_for_claim[0]["leave_one_out_max_r"])) > 1e-10):
            raise ValueError(f"Leave-one-out profile summary disagrees with rows for {claim_id}")
        profile_loo_summary[claim_id] = {"minimum": minimum, "maximum": maximum}
    if observed_profile_claims != configured_profile_claims:
        raise ValueError("Leave-one-out profile claims differ from the frozen plan")
    for record in placo:
        result_path = ROOT / record["output_path"]
        if not result_path.is_file() or digest(result_path) != record["output_sha256"]:
            raise ValueError(f"Published PLACO result does not match its master checksum: {result_path}")
    if len(global_rows) != 72:
        raise ValueError(f"Expected the locked 72-row global map, got {len(global_rows)}")
    sig = [r for r in global_rows if r["significance_under_original_396_family"].lower() == "true"]
    if len(sig) != 35:
        raise ValueError(f"Expected 35 inherited significant rows, got {len(sig)}")
    if len(context_rows) != 1 or context_rows[0]["admission_status"] != "CONTEXT_ONLY_NOT_PHASE5_REPLICATION":
        raise ValueError("Published long-sleep context must remain a single non-replication comparison")

    evidence_global = evidence(global_path, lava_path)
    evidence_profiles = evidence(global_path, profile_path, lava_path, profile_loo_path,
                                 profile_loo_provenance_path, profile_loo_plan_path)
    evidence_replication = evidence(replication_path, lava_path)
    evidence_loci = evidence(placo_path, locus_path, lava_path, mhc_plan_path, mhc_script_path,
                             mhc_path, mhc_provenance_path)
    evidence_context = evidence(context_source_path, context_table_path)
    rows: list[dict[str, str]] = []
    mhc_positioned = int(mhc_summary["candidate_variants_position_validated"])
    mhc_inside = int(mhc_summary["candidate_variants_inside_mhc"])
    mhc_unknown = int(mhc_summary["candidate_variants_unpositioned_or_mismatched"])
    mhc_intervals = int(mhc_summary["candidate_intervals_all"])
    mhc_overlap_intervals = int(mhc_summary["candidate_intervals_overlapping_mhc"])
    mhc_filter_status = ("NO_CANDIDATES_IN_INTERVAL" if mhc_inside == 0 and mhc_overlap_intervals == 0
                         else "CANDIDATES_OVERLAP_INTERVAL")
    mhc_filter_summary = (
        f"{mhc_filter_status}; {mhc_inside}/{mhc_positioned} position-validated candidates and "
        f"{mhc_overlap_intervals}/{mhc_intervals} candidate intervals overlap GRCh37 chr6:25-34 Mb; "
        f"{mhc_unknown} candidate positions remain unmatched; four pairs only"
    )
    mhc_by_pair = {row["pair_id"]: row for row in mhc_pair_rows}
    evidence_placo_mhc = evidence(placo_path, lava_path, mhc_plan_path, mhc_script_path,
                                  mhc_path, mhc_provenance_path)

    rows.append(status_defaults(
        claim_id="global_map_35_of_72",
        manuscript_claim="The locked 12-by-6 map contains 72 relationships, of which 35 retain significance under the original 396-test FDR family.",
        evidence_status="SUPPORTED_LOCKED_ATLAS_EXTRACTION",
        claim_scope="Aggregate descriptive map; no new discovery-family correction",
        evidence_sources=evidence_global,
        original_396_family_significance="PASS_INHERITED_FROM_LOCKED_ATLAS",
        LAVA_UKB_v1_1="IN_PROGRESS_FOR_FIVE_PRIORITIZED_PAIRS",
        overall_sensitivity_status="LOCKED_EXTRACTION_VERIFIED; DOWNSTREAM_ROBUSTNESS_INCOMPLETE",
        interpretation_limit="Does not test robustness to alternate GWAS or establish variant sharing/causality.",
    ))

    alz = [r for r in global_rows if r["brain_disorder"].lower() == "alz"]
    alz_sig = sum(r["significance_under_original_396_family"].lower() == "true" for r in alz)
    rows.append(status_defaults(
        claim_id="alzheimer_no_inherited_signals",
        manuscript_claim=f"No Alzheimer's pair in the locked 12-trait subset passes the original 396-family FDR ({alz_sig}/12 significant).",
        evidence_status="SUPPORTED_LOCKED_ATLAS_NULL_WITHIN_PANEL",
        claim_scope="Locked atlas subset only",
        evidence_sources=evidence_global,
        original_396_family_significance="PASS_INHERITED_NULL_WITHIN_THIS_PANEL",
        overall_sensitivity_status="LOCKED_NULL_VERIFIED; GENERALIZABILITY_UNTESTED",
        interpretation_limit="Not evidence that all sleep phenotypes or Alzheimer's GWAS have no genetic sharing.",
    ))

    ins_mdd = next(r for r in global_rows if r["sleep_trait"] == "insomnia" and r["brain_disorder"] == "mdd")
    rows.append(status_defaults(
        claim_id="insomnia_mdd_largest_global_rg",
        manuscript_claim=f"Insomnia-MDD has the largest absolute global estimate in the map (rg={float(ins_mdd['rg']):.4f}).",
        evidence_status="SUPPORTED_LOCKED_ATLAS_ESTIMATE",
        claim_scope="Global rg ranking among the 72 locked rows",
        evidence_sources=evidence_global,
        independent_replication="NO_RESULT_IN_CURRENT_PACKAGE",
        sample_overlap_sensitivity="UKB_COMPONENTS_IN_BOTH; PARTICIPANT_OVERLAP_UNVERIFIED; SENSITIVITY_NOT_RUN",
        original_396_family_significance="PASS_INHERITED_SIGNIFICANT",
        LAVA_UKB_v1_1="IN_PROGRESS; NO_FAMILY_RESULT",
        overall_sensitivity_status="GLOBAL_ESTIMATE_VERIFIED; ROBUSTNESS_INCOMPLETE",
        interpretation_limit="Selected post-atlas for deep follow-up; overlap and selection limit interpretation.",
    ))

    # Preserve the manuscript's two cited comparisons and their exact rows.
    profile_matrix = {r["brain_disorder"]: r for r in profiles}
    profile_claims = [("profile_adhd_mdd", "adhd", "mdd"),
                      ("profile_scz_bipolar", "scz", "bipolar")]
    for claim_id, a, b in profile_claims:
        value = float(profile_matrix[a][b])
        loo = profile_loo_summary[claim_id]
        claim = f"{a.upper()}-{b.upper()} disorder-profile similarity is descriptive (r={value:.3f})."
        rows.append(status_defaults(
            claim_id=claim_id, manuscript_claim=claim,
            evidence_status="SUPPORTED_DESCRIPTIVE_PROFILE_STATISTIC",
            claim_scope="Pearson similarity across shared 12 sleep traits",
            evidence_sources=evidence_profiles,
            profile_leave_one_sleep_trait_out=(
                f"12 omissions; Pearson r range {loo['minimum']:.3f}–{loo['maximum']:.3f}; descriptive only"
            ),
            original_396_family_significance="NOT_APPLICABLE_TO_PROFILE_CORRELATION",
            overall_sensitivity_status="DESCRIPTIVE_ONLY; LOO_RANGES_REPORTED; NO_INDEPENDENT_VALIDATION",
            interpretation_limit=(f"Reuses the same sleep-trait profiles; leave-one-out is not independent validation "
                                  f"and does not estimate uncertainty for {a}-{b}."),
        ))

    adhd_rep = next(r for r in reps if r["sleep_trait"] == "insomnia" and r["brain_disorder"] == "adhd")
    rows.append(status_defaults(
        claim_id="insomnia_adhd_directional_replication",
        manuscript_claim=f"Insomnia-ADHD is directionally concordant in FinnGen R13 (discovery rg={adhd_rep['discovery_rg']}; replication rg={adhd_rep['replication_rg']}).",
        evidence_status="SUPPORTED_DIRECTIONAL_REPLICATION_LEVEL_B",
        claim_scope="Same sleep GWAS retained; external FinnGen ADHD GWAS",
        evidence_sources=evidence_replication,
        independent_replication="PARTIAL_LEVEL_B_DIRECTIONAL_ONLY",
        sample_overlap_sensitivity="EXTERNAL_FINNGEN_COHORT; PARTICIPANT_LEVEL_CHECK_NOT_PROVIDED",
        phenotype_definition_sensitivity="COMPARABLE_REGISTER_ADHD; NOT_EXACT_PHENOTYPE",
        original_396_family_significance="PASS_INHERITED_SIGNIFICANT",
        overall_sensitivity_status="DIRECTION_REPLICATED; FULL_INDEPENDENCE_AND_HETEROGENEITY_UNTESTED",
        interpretation_limit="Do not call fully independent replication: discovery insomnia GWAS is reused; Finnish phenotype/population differ.",
    ))

    context = context_rows[0]
    rows.append(status_defaults(
        claim_id="longsleep_scz_published_context",
        manuscript_claim=(f"Published long-sleep-SCZ context is directionally concordant (rg={context['published_rg']}, "
                          f"SE={context['published_se']}, P={context['published_p']}); it is not independent replication."),
        evidence_status="SUPPORTED_PUBLISHED_CONTEXT_ONLY_LEVEL_D",
        claim_scope="Published comparable-phenotype estimate; not a Brain6 reanalysis",
        evidence_sources=evidence_context,
        independent_replication="LEVEL_D_CONTEXT_ONLY; UKB_SLEEP_OVERLAP_POSSIBLE; EXACT_PGC_CROSSWALK_UNRESOLVED",
        sample_overlap_sensitivity="UKB_IN_PUBLISHED_SLEEP_SAMPLE; PARTICIPANT_OVERLAP_POSSIBLE; NOT_TESTED",
        phenotype_definition_sensitivity="PUBLISHED_LONG_SLEEP_GE10H; BRAIN6_LONG_SLEEP_GE9H; COMPARABLE_NOT_IDENTICAL",
        original_396_family_significance="NOT_APPLICABLE_TO_EXTERNAL_PUBLISHED_ESTIMATE",
        overall_sensitivity_status="CONTEXTUAL_DIRECTIONAL_AGREEMENT; NOT_REPLICATION; OVERLAP_LIMITED",
        interpretation_limit="Published rg values are not combined; the covariance between estimates is unknown.",
    ))

    for pair in sorted(placo, key=lambda r: r["pair_id"]):
        if pair["stage_status"].startswith("COMPLETE_QC_PASS"):
            overlap_record = overlap_by_pair[pair["pair_id"]]
            overlap = (f"{overlap_record['overlap_classification']}; "
                       f"{overlap_record['evidence_and_limit']}; SENSITIVITY_NOT_RUN")
            pair_mhc = mhc_by_pair[pair["pair_id"]]
            pair_mhc_status = ("NO_CANDIDATES_IN_INTERVAL" if
                               int(pair_mhc["candidate_variants_inside_mhc"]) == 0 and
                               int(pair_mhc["candidate_intervals_overlapping_mhc"]) == 0
                               else "CANDIDATES_OVERLAP_INTERVAL")
            pair_mhc_sensitivity = (
                f"{pair_mhc_status}; {pair_mhc['candidate_variants_inside_mhc']}/"
                f"{pair_mhc['candidate_variants_position_validated']} position-validated candidates and "
                f"{pair_mhc['candidate_intervals_overlapping_mhc']}/{pair_mhc['candidate_intervals_all']} "
                "candidate intervals overlap GRCh37 chr6:25-34 Mb; "
                f"{pair_mhc['candidate_variants_unpositioned_or_mismatched']} positions unmatched; "
                "candidate set only"
            )
            rows.append(status_defaults(
                claim_id=f"placo_{pair['pair_id']}_headline_variants",
                manuscript_claim=f"PLACO flagged {pair['headline_variants']} variants below the frozen {pair['across_track_p_threshold']} threshold for {pair['sleep_trait']}-{pair['brain_disorder']}.",
                evidence_status="PAIR_QC_PASS_PARTIAL_FAMILY",
                claim_scope=("Pair-level QC passed; full five-track correction family remains incomplete"
                             if not loci["family_complete"] else "Pair-level QC passed within the complete frozen family"),
                evidence_sources=evidence_placo_mhc,
                independent_replication="NOT_ATTEMPTED",
                sample_overlap_sensitivity=overlap,
                original_396_family_significance="PASS_INHERITED_FOR_GLOBAL_PAIR",
                LAVA_UKB_v1_1="IN_PROGRESS; NO_LOCAL_RESULT",
                MHC_exclusion=pair_mhc_sensitivity,
                overall_sensitivity_status=f"PAIR_RESULT_VERIFIED; MHC_{pair_mhc_status}_FOR_PARTIAL_CANDIDATE_SET; FULL_FAMILY_AND_LOCUS_SENSITIVITY_INCOMPLETE",
                interpretation_limit="PLACO association does not establish independent loci, causal variants, or biological pleiotropy.",
            ))

    rows.append(status_defaults(
        claim_id="placo_partial_candidate_regions",
        manuscript_claim=(f"The partial-family candidate grouping across {len(loci['pairs_included'])} published "
                          f"QC-passed PLACO pairs ({', '.join(loci['pairs_included'])}) contains "
                          f"{loci['n_candidate_variants']} variants in {loci['n_candidate_loci']} "
                          "UKB-reference intervals."),
        evidence_status="PASS_PARTIAL_FAMILY_PLACO_ONLY",
        claim_scope=("All required PLACO tracks and protected Track B are included"
                     if loci["family_complete"] else
                     "Four v3 PLACO pairs only; protected Track B is absent and full-family deduplication remains pending"),
        evidence_sources=evidence_loci,
        independent_replication="NOT_TESTED_FOR_CANDIDATE_REGIONS",
        sample_overlap_sensitivity="PAIRWISE_OVERLAP_CAVEATS_RETAINED; NO_SENSITIVITY",
        original_396_family_significance="LEAD PAIRS_INHERITED_SIGNIFICANT",
        LAVA_UKB_v1_1="LD_REFERENCE_USED_FOR_CLUMPING; LOCAL_RG_FAMILY_IN_PROGRESS",
        alternate_LD_reference="NOT_RUN",
        MHC_exclusion=mhc_filter_summary,
        overall_sensitivity_status=f"PRELIMINARY_CANDIDATE_GROUPING_ONLY; MHC_{mhc_filter_status}_FOR_OBSERVED_PARTIAL_SET; NOT_FINAL_LOCI",
        interpretation_limit="MHC filter is descriptive for four PLACO-only candidate outputs and does not establish full-family, local-rg, fine-mapping, colocalization, or molecular robustness.",
    ))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    temp = OUT.with_suffix(OUT.suffix + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(OUT)
    print(json.dumps({"path": str(OUT), "claims": len(rows), "evidence_status": "INTERIM"}, indent=2))


if __name__ == "__main__":
    main()
