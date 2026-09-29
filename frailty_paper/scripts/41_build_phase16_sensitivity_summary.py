#!/usr/bin/env python3
"""Build a source-linked conclusion-by-sensitivity status matrix.

The matrix separates sensitivities that were actually evaluated from those
that are unavailable, unresolved, blocked, or not justified by upstream data.
It does not rerun a GWAS or create new statistical thresholds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


INPUTS = {
    "fi_frozen": "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv",
    "fi_rerun": "frailty_paper/results/frailty_v1/global_rg_sleep_frailty.tsv",
    "fi_correction": "frailty_paper/analysis/frozen_fi_correction_sensitivity.tsv",
    "latent_correction": "frailty_paper/analysis/latent_factor_correction_sensitivity.tsv",
    "h2_master": "frailty_paper/results/frailty_v1/ldsc/h2_master.tsv",
    "overlap": "frailty_paper/analysis/sample_overlap_assessment.tsv",
    "replication": "frailty_paper/analysis/replication_resource_audit.tsv",
    "physical_qc": "frailty_paper/manifests/physical_component_qc.tsv",
    "analysis_plan": "frailty_paper/config/analysis_plan_v1.yaml",
}

OUTPUT = "frailty_paper/analysis/phase16_sensitivity_conclusion_matrix.tsv"
PROVENANCE = "frailty_paper/analysis/phase16_sensitivity_conclusion_matrix.provenance.json"

COLUMNS = [
    "conclusion_id", "headline_conclusion", "sensitivity_domain",
    "sensitivity_check", "assessment", "observed_result",
    "conclusion_survival", "interpretation", "evidence_files",
]

CLAIMS = [
    ("FI_INSOMNIA", "Insomnia has positive global genetic correlation with the registered FI.", "insomnia"),
    ("FI_SHORTSLEEP", "Short sleep (<7 h) has positive global genetic correlation with the registered FI.", "shortsleep"),
    ("FI_SLEEP_APNEA", "Sleep apnea has positive global genetic correlation with the registered FI.", "sleep_apnea"),
    ("FI_SLEEP_DURATION", "Continuous sleep duration has negative global genetic correlation with the registered FI.", "sleepdur"),
    ("FI_SLEEP_EFFICIENCY", "Actigraphy sleep efficiency has negative global genetic correlation with the registered FI.", "sleep_efficiency"),
    ("FI_PANEL", "Nine of 12 registered sleep/FI rg estimates pass the inherited all-396 BH correction.", "__panel__"),
    ("LATENT_PANEL", "Forty-eight of 84 secondary sleep/latent-factor rg estimates pass the fixed 84-pair BH correction.", "__latent__"),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def truth(value: str) -> bool:
    return value.strip().upper() in {"TRUE", "1", "YES", "PASS"}


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(stream.getvalue(), encoding="utf-8")


def build(repo: Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    paths = {key: repo / rel for key, rel in INPUTS.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing sensitivity-matrix source(s): {missing}")
    src = {key: read_tsv(path) for key, path in paths.items() if path.suffix == ".tsv"}

    fi = src["fi_frozen"]
    rerun = src["fi_rerun"]
    fi_corr = src["fi_correction"]
    latent = src["latent_correction"]
    h2 = src["h2_master"]
    overlap = [r for r in src["overlap"] if r["frailty_endpoint"] == "frailty"]
    replication = src["replication"]
    physical_qc = src["physical_qc"]
    if (len(fi), len(rerun), len(fi_corr), len(latent), len(h2), len(overlap)) != (12, 12, 12, 84, 20, 12):
        raise ValueError("Source families changed: expected FI=12, rerun=12, FI correction=12, latent=84, h2=20, overlap=12")
    if len({r["sleep_trait"] for r in fi}) != 12 or len({r["sleep_trait"] for r in overlap}) != 12:
        raise ValueError("FI/overlap panel is not one-to-one")
    fi_corr_by_sleep = {r["sleep_trait"]: r for r in fi_corr}
    fi_by_sleep = {r["sleep_trait"]: r for r in fi}
    rerun_by_sleep = {r["sleep_trait"]: r for r in rerun}
    latent_general = {r["sleep_trait"]: r for r in latent if r["frailty_factor"] == "frailty_general"}
    h2_by_trait = {r["trait_id"]: r for r in h2}
    if len(latent_general) != 12 or len(h2_by_trait) != 20:
        raise ValueError("Missing one-to-one latent general-factor or heritability rows")

    fi_bh = sum(truth(r["bh_all_396_significant_at_0.05"]) for r in fi_corr)
    fi_bonf = sum(truth(r["bonferroni_all_396_significant_at_0.05"]) for r in fi_corr)
    latent_bh = sum(truth(r["bh_significant_at_0.05"]) for r in latent)
    latent_bonf = sum(truth(r["bonferroni_significant_at_0.05"]) for r in latent)
    if (fi_bh, fi_bonf, latent_bh, latent_bonf) != (9, 9, 48, 37):
        raise ValueError("Observed correction counts no longer match the locked sensitivity outputs")

    exact_rg = sum(float(fi_by_sleep[name]["global_rg"]) == float(rerun_by_sleep[name]["rg"]) for name in fi_by_sleep)
    exact_p = sum(float(fi_by_sleep[name]["global_rg_p"]) == float(rerun_by_sleep[name]["p"]) for name in fi_by_sleep)
    max_abs_delta = max(abs(float(fi_by_sleep[name]["global_rg"]) - float(rerun_by_sleep[name]["rg"])) for name in fi_by_sleep)
    max_relative_p_delta = max(abs(float(fi_by_sleep[name]["global_rg_p"]) - float(rerun_by_sleep[name]["p"])) / max(abs(float(fi_by_sleep[name]["global_rg_p"])), abs(float(rerun_by_sleep[name]["p"]))) for name in fi_by_sleep)
    if exact_rg != 12 or exact_p != 11 or max_abs_delta > 1e-12 or max_relative_p_delta > 1e-14:
        raise ValueError("Frozen-versus-rerun FI rg concordance no longer matches the documented audit")

    relevant_h2 = [r for r in h2 if r["trait_id"] in {"frailty", *fi_by_sleep.keys()}]
    relevant_factor_h2 = [r for r in h2 if r["trait_id"] == "frailty_general" or r["trait_id"].startswith("frailty_factor_")]
    fi_h2_pass = len(relevant_h2) == 13 and all(r["verdict"] == "PASS" and truth(r["pass_z"]) and truth(r["pass_intercept"]) for r in relevant_h2)
    latent_h2_pass = len(relevant_factor_h2) == 7 and all(r["verdict"] == "PASS" and truth(r["pass_z"]) and truth(r["pass_intercept"]) for r in relevant_factor_h2)

    actigraphy = [fi_by_sleep[name] for name in ("accel_sleep_duration", "sleep_efficiency", "sleep_timing")]
    actigraphy_sig = sum(float(r["global_rg_fdr_all_396"]) <= 0.05 for r in actigraphy)
    actigraphy_sig_names = ", ".join(r["sleep_trait"] for r in actigraphy if float(r["global_rg_fdr_all_396"]) <= 0.05) or "none"
    ukb_expected = sum("EXPECTED" in r["cohort_overlap_status"].upper() for r in overlap)
    overlap_exact_unknown = all(r["exact_participant_overlap"].upper() == "UNKNOWN" for r in overlap)
    no_pairwise_replication = not any("PAIRWISE_REPLICATION" in r["pairwise_replication_eligibility"] and "NOT" not in r["pairwise_replication_eligibility"] for r in replication)
    if not physical_qc:
        raise ValueError("Physical-component QC register is empty")

    ev = "; ".join(INPUTS[key] for key in INPUTS)
    rows: list[dict[str, str]] = []

    def add(claim_id: str, headline: str, domain: str, check: str, assessment: str,
            result: str, survival: str, interpretation: str, evidence: str = ev) -> None:
        rows.append({"conclusion_id": claim_id, "headline_conclusion": headline,
                     "sensitivity_domain": domain, "sensitivity_check": check,
                     "assessment": assessment, "observed_result": result,
                     "conclusion_survival": survival, "interpretation": interpretation,
                     "evidence_files": evidence})

    for claim_id, headline, sleep in CLAIMS:
        is_fi_pair = sleep not in {"__panel__", "__latent__"}
        if sleep == "__latent__":
            add(claim_id, headline, "alternate_multiple_testing", "BH versus same-family Bonferroni",
                "COMPLETED", f"BH: {latent_bh}/84; Bonferroni: {latent_bonf}/84",
                "PARTIAL_FAMILY_SURVIVAL", "The latent family retains 37 associations under its locked 84-pair Bonferroni sensitivity; this remains secondary and overlap-limited.")
        else:
            relevant = [fi_corr_by_sleep[sleep]] if is_fi_pair else fi_corr
            bh_n = sum(truth(r["bh_all_396_significant_at_0.05"]) for r in relevant)
            bonf_n = sum(truth(r["bonferroni_all_396_significant_at_0.05"]) for r in relevant)
            denom = len(relevant)
            add(claim_id, headline, "alternate_multiple_testing", "Locked BH q-values versus all-396 Bonferroni",
                "COMPLETED", f"BH-significant {bh_n}/{denom}; all-396 Bonferroni-significant {bonf_n}/{denom}",
                "SURVIVES" if bonf_n == denom else "PARTIAL_OR_NOT_SURVIVING",
                "The point estimates are unchanged; the stricter correction changes only the significance decision.")

        if is_fi_pair or sleep == "__panel__":
            add(claim_id, headline, "workflow_reproducibility", "Frozen atlas versus separate frailty_v1 LDSC rerun",
                "COMPLETED", f"12/12 rows compared; rg values identical at saved precision; P values numerically identical for {exact_p}/12 (maximum relative difference {max_relative_p_delta:.3g})",
                "WORKFLOW_CONCORDANT", "Same source GWAS and analysis family; this is reproducibility evidence, not independent replication.")
        else:
            add(claim_id, headline, "workflow_reproducibility", "Independent rerun of the latent-factor family",
                "NOT_RUN", "No separate latent-family rerun registered", "NOT_TESTED",
                "Do not imply independent workflow reproduction for the latent-factor family.")

        gate_pass = latent_h2_pass if sleep == "__latent__" else fi_h2_pass
        n_gate = 7 if sleep == "__latent__" else 13
        add(claim_id, headline, "weak_h2_exclusion", "Apply frozen h² Z/intercept input gate",
            "COMPLETED_INPUT_GATE", f"{n_gate}/{n_gate} relevant traits pass; no relevant trait excluded by this gate" if gate_pass else "One or more relevant traits fail the frozen h² gate",
            "INPUTS_PASS" if gate_pass else "INPUT_GATE_FAILS",
            "This records the pre-analysis eligibility gate; it is not an additional effect-estimate sensitivity analysis.")

        add(claim_id, headline, "alternate_frailty_definitions", "Compare FI with Fried score and HFRS",
            "BLOCKED", "No eligible paired sleep–Fried or sleep–HFRS rg estimate", "NOT_TESTED",
            "The Fried file is not fully acquired and the exact HFRS endpoint/file is unverified; no cross-definition robustness claim is supported.")
        add(claim_id, headline, "physical_frailty_components", "Compare with five Fried components",
            "BLOCKED_SOURCE_PROVENANCE", "No eligible component sleep-pair estimates", "NOT_TESTED",
            "Component files are not eligible for harmonization because sufficient build/effect/model/generating-study provenance is missing.")

        if sleep == "__latent__":
            add(claim_id, headline, "latent_frailty_factors", "Secondary 12×7 latent-factor family",
                "COMPLETED_SENSITIVITY_FAMILY", f"BH {latent_bh}/84; Bonferroni {latent_bonf}/84; exact participant overlap UNKNOWN",
                "PARTIAL_WITH_OVERLAP_LIMIT", "This is a separate secondary family and does not identify a factor driving FI overlap.")
        elif sleep == "__panel__":
            add(claim_id, headline, "latent_frailty_factors", "Secondary 12×7 latent-factor family",
                "COMPLETED_SENSITIVITY_FAMILY", f"BH {latent_bh}/84; Bonferroni {latent_bonf}/84; exact participant overlap UNKNOWN",
                "PARTIAL_CONVERGENCE", "Related frailty-factor evidence exists, but it is not a direct alternate-FI test and does not resolve cohort overlap.")
        else:
            lr = latent_general[sleep]
            sign_match = (float(fi_by_sleep[sleep]["global_rg"]) * float(lr["rg"])) > 0
            lat_sig = truth(lr["bonferroni_significant_at_0.05"])
            add(claim_id, headline, "latent_frailty_factors", "Compare with the general latent frailty factor",
                "COMPLETED_DIFFERENT_ENDPOINT", f"General factor rg={lr['rg']}; 84-pair Bonferroni P={lr['bonferroni_p_secondary_84']}; sign concordance={sign_match}",
                "PARTIAL_CONVERGENCE" if sign_match and lat_sig else "NO_CLEAR_CONVERGENCE",
                "A related latent phenotype is not the same FI endpoint, and estimates are not a statistical test of heterogeneity; exact overlap remains unknown.")

        add(claim_id, headline, "independent_replication", "Require paired independent sleep and frailty estimate",
            "NOT_ESTABLISHED" if no_pairwise_replication else "CANDIDATE_REQUIRES_REVIEW",
            "No eligible independent pairwise sleep–frailty rg replication established" if no_pairwise_replication else "At least one candidate row needs eligibility review",
            "NOT_ESTABLISHED" if no_pairwise_replication else "UNRESOLVED",
            "Selected loci/PRS, outcome-only sources, and overlapping meta-analyses do not replicate a paired global rg.")
        add(claim_id, headline, "exact_vs_comparable_replication", "Classify eligible pairwise results as exact or comparable phenotype replication",
            "CANDIDATE_MATCHES_AUDITED_NO_PAIRWISE_ESTIMATE",
            "Phenotype-match distinctions are recorded for candidate sources, but no eligible paired estimate has been classified as exact or comparable replication",
            "NOT_ESTABLISHED",
            "Do not promote candidate phenotype similarity, lead-locus validation, or outcome-side validation to a pairwise replication label.")
        add(claim_id, headline, "UKB_overlap", "Assess exact participant overlap and cohort-level UKB sharing",
            "UNRESOLVED", f"UKB overlap expected in {ukb_expected}/12 primary pairs; exact overlap unknown for all 12" if overlap_exact_unknown else "See pairwise overlap crosswalk",
            "UNRESOLVED", "No overlap-adjusted or UKB-excluded paired rg estimate is available; unknown is not zero.")

        if sleep == "__latent__":
            objective_result = f"Primary FI panel includes 3 actigraphy traits; {actigraphy_sig}/3 retain inherited q≤0.05: {actigraphy_sig_names}"
            objective_assess, objective_survival = "DESCRIPTIVE_ONLY", "NOT_DIRECTLY_COMPARABLE"
        elif sleep == "__panel__":
            objective_result = f"Actigraphy FI pairs: {actigraphy_sig}/3 retain inherited q≤0.05 ({actigraphy_sig_names})"
            objective_assess, objective_survival = "DESCRIPTIVE_ONLY", "PARTIAL_MEASUREMENT_CONTRAST"
        else:
            objective_result = f"No direct objective counterpart for `{sleep}`; actigraphy FI pairs overall: {actigraphy_sig}/3 pass inherited q≤0.05 ({actigraphy_sig_names})"
            objective_assess, objective_survival = "NOT_DIRECTLY_MATCHED", "NOT_TESTED"
        add(claim_id, headline, "self_report_vs_objective_sleep", "Compare comparable self-report and actigraphy sleep measures",
            objective_assess, objective_result, objective_survival,
            "Actigraphy phenotypes are distinct measures, not replications of insomnia, sleep apnea, or tail-defined sleep duration.")

        add(claim_id, headline, "ancestry", "Evaluate ancestry-matched replication/LD support",
            "NOT_AVAILABLE", "Current sleep–frailty genetic-correlation families are European ancestry", "NOT_TESTED",
            "No ancestry-diverse paired frailty estimate with matched ancestry-specific LD was identified in the eligible inputs.")
        add(claim_id, headline, "independent_cohort_sensitivity", "Repeat the paired analysis in an independent cohort with compatible phenotypes",
            "NOT_AVAILABLE", "No complete independent cohort pair with compatible sleep and frailty summary statistics is eligible",
            "NOT_TESTED", "HRS/ELSA and CHARGE candidates are outcome-only or selected-locus/score sources; UKB–MVP public meta-analysis overlaps UKB.")
        for domain, check in [
            ("LD_reference", "Re-estimate eligible local results across appropriate LD references"),
            ("fine_mapping_stability", "Compare locked fine-mapping settings/LD sensitivity"),
            ("coloc_prior_sensitivity", "Vary locked trait-trait coloc priors"),
            ("MHC_exclusion", "Re-evaluate eligible locus analyses excluding MHC"),
        ]:
            add(claim_id, headline, domain, check, "NOT_JUSTIFIED_UPSTREAM_GATES_FAIL",
                "No frailty-specific LAVA-prioritized loci or downstream locus analysis", "NOT_TESTED",
                "Do not create a post hoc locus-based sensitivity analysis without eligible shared-locus inputs.")
        add(claim_id, headline, "leave_one_dataset_out", "Repeat paired rg while removing one contributing cohort/dataset",
            "NOT_AVAILABLE", "No cohort-specific FI summary-statistic decomposition available; exact intersections unknown", "NOT_TESTED",
            "A UKB-only sleep result and a UKB+TwinGene FI meta-analysis do not provide sufficient cohort-specific inputs for this sensitivity.")

    if len(rows) != len(CLAIMS) * 17:
        raise ValueError(f"Unexpected sensitivity matrix size: {len(rows)}")
    return rows, {key: {"path": rel, "sha256": sha256(paths[key])} for key, rel in INPUTS.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--replace-outputs", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    rows, inputs = build(repo)
    out = repo / OUTPUT
    prov = repo / PROVENANCE
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    data = stream.getvalue().encode("utf-8")
    if out.exists() and out.read_bytes() != data and not args.replace_outputs:
        raise SystemExit(f"Refusing to overwrite non-identical output: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    provenance = {
        "schema_version": "frailty_phase16_sensitivity_matrix.v1",
        "script": Path(__file__).resolve().relative_to(repo).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "command": "python frailty_paper/scripts/41_build_phase16_sensitivity_summary.py --repo <repo>",
        "inputs": inputs,
        "output": {"path": OUTPUT, "sha256": sha256(out), "rows": len(rows)},
        "claims": len(CLAIMS),
        "sensitivity_domains_per_claim": 17,
        "assessment_legend": {
            "SURVIVES": "The stated conclusion retains its locked significance decision under a completed sensitivity.",
            "PARTIAL*": "Related evidence is available but is not a direct test of the same estimand or endpoint.",
            "NOT_TESTED": "The test is unavailable, unresolved, blocked, or not justified; do not interpret as a null result.",
            "WORKFLOW_CONCORDANT": "Same-source reproducibility evidence only; not independent replication.",
        },
        "limitations": [
            "This status matrix does not create new effect estimates or relax frozen thresholds.",
            "All blocked/not-tested cells remain open sensitivity domains, not evidence of robustness or absence.",
            "Latent factors and actigraphy are different endpoints; their results are descriptive cross-phenotype checks, not equivalence tests.",
        ],
    }
    prov_data = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if prov.exists() and prov.read_bytes() != prov_data and not args.replace_outputs:
        raise SystemExit(f"Refusing to overwrite non-identical provenance: {prov}")
    prov.write_bytes(prov_data)
    print(json.dumps({"output": OUTPUT, "rows": len(rows), "sha256": sha256(out), "provenance": PROVENANCE}, indent=2))


if __name__ == "__main__":
    main()
