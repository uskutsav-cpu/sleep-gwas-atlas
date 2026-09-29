#!/usr/bin/env python3
"""Build source-linked supplementary tables from currently available outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


INPUTS = {
    "resource_manifest": "frailty_paper/manifests/all_acquired_resources.tsv",
    "h2_master": "frailty_paper/results/frailty_v1/ldsc/h2_master.tsv",
    "fi_rg": "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv",
    "fi_overlap": "frailty_paper/analysis/sample_overlap_assessment.tsv",
    "aging_rg": "frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv",
    "aging_overlap": "frailty_paper/analysis/frozen_atlas_aging_context_overlap.tsv",
    "latent_rg": "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv",
    "replication_audit": "frailty_paper/analysis/replication_resource_audit.tsv",
    "fi_correction": "frailty_paper/analysis/frozen_fi_correction_sensitivity.tsv",
    "latent_correction": "frailty_paper/analysis/latent_factor_correction_sensitivity.tsv",
    "sensitivity_summary": "frailty_paper/analysis/phase16_sensitivity_conclusion_matrix.tsv",
}
TOOL_PINS = "environment/tool_versions.tsv"
VERSION_EVIDENCE = {
    "fi_h2_reproduction": "frailty_paper/analysis/frailty_index_h2_reproduction_manifest_2026-09-23.json",
    "figure1": "frailty_paper/analysis/figure1_evidence_framework.provenance.json",
    "figure2": "frailty_paper/analysis/figure2_global_rg_matrix.provenance.json",
    "figure4_current": "frailty_paper/analysis/figures4_5_downstream_evidence_status.provenance.json",
    "fi_forest_preview": "frailty_paper/analysis/preliminary_frozen_fi_global_rg_forest.provenance.json",
    "latent_heatmap_preview": "frailty_paper/analysis/sleep_latent_frailty_rg_heatmap.provenance.json",
    "dimension_forest_preview": "frailty_paper/analysis/insomnia_sleep_apnea_frailty_dimensions_forest.provenance.json",
    "workflow_validation": "frailty_paper/analysis/pinned_workflow_environment_validation_2026-09-23.json",
    "clean_venv_validation": "frailty_paper/analysis/clean_venv_setup_validation_2026-09-23.json",
    "ldsc_fi_log": "frailty_paper/results/frailty_v1/logs_fi_verification/h2_frailty_reproduction_2026-09-23.log",
    "integrated_preflight": "frailty_paper/analysis/preflight_reverification_2026-09-24_350_94.log",
    "latest_integrated_preflight": "frailty_paper/analysis/preflight_current_head_py311_2026-09-24.log",
    "post_s16_tests": "frailty_paper/analysis/current_head_py311_unittest_2026-09-24.log",
    "latest_suite": "frailty_paper/analysis/test_suite_frailty_py311_2026-09-27_0320.log",
    "current_canonical_suite": "frailty_paper/analysis/make_test_unittest_2026-09-28_2112.log",
    "bibliography_audit_latest": "frailty_paper/analysis/manuscript_bibliography_audit_2026-09-28_2052.json",
    "quantitative_claims_audit_latest": "frailty_paper/analysis/manuscript_quantitative_claims_audit_2026-09-28_2052.json",
    "reporting_locator_audit_latest": "frailty_paper/analysis/reporting_locator_audit_2026-09-28_2052.json",
}

TABLE_FILES = {
    "S1": "table_s1_gwas_metadata.tsv",
    "S4": "table_s4_snp_heritability.tsv",
    "S5": "table_s5_global_rg.tsv",
    "S7": "table_s7_latent_factor_rg.tsv",
    "S8": "table_s8_replication_audit.tsv",
    "S15": "table_s15_correction_sensitivity.tsv",
    "S15_matrix": "table_s15_sensitivity_conclusion_matrix.tsv",
    "S16": "table_s16_software_resources_versions.tsv",
}

S5_COLUMNS = [
    "analysis_family", "analysis_class", "pair_id", "sleep_trait", "outcome_trait",
    "rg", "se", "z", "p_value", "q_value", "q_family", "family_denominator",
    "direction", "cohort_overlap_status", "exact_participant_overlap",
    "replication_status", "interpretation_status", "claim_limit", "source_table",
]

S7_COLUMNS = [
    "sleep_trait", "disease_trait", "rg", "se", "z", "p_value", "h2_obs",
    "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se",
    "input_log", "q_value", "analysis_family", "family_denominator",
    "cohort_overlap_status", "interpretation_status", "claim_limit",
    "indicator_overlap_status", "indicator_overlap_note",
]

S15_COLUMNS = [
    "analysis_family", "pair_id", "sleep_trait", "outcome_trait", "rg", "se", "raw_p",
    "bh_q", "bh_family_denominator", "bonferroni_family_size", "bonferroni_p",
    "bh_significant_at_0.05", "bonferroni_significant_at_0.05", "cohort_overlap_status",
    "exact_participant_overlap", "comparison_interpretation", "source_table",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def tsv_bytes(rows: list[dict[str, str]], columns: list[str] | None = None) -> bytes:
    if columns is None:
        if not rows:
            raise ValueError("Cannot infer columns from an empty table")
        columns = list(rows[0])
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def gwas_manifest_scope(path: Path) -> tuple[str, int]:
    """Hash only the GWAS-summary rows consumed by Supplementary Table S1."""
    rows = [r for r in read_tsv(path) if r["resource_type"] == "GWAS summary statistics"]
    if not rows:
        raise ValueError(f"No GWAS summary-statistic rows found in {path}")
    return hashlib.sha256(tsv_bytes(rows)).hexdigest(), len(rows)


def software_version_table(repo: Path) -> list[dict[str, str]]:
    pin_path = repo / TOOL_PINS
    version_paths = {key: repo / value for key, value in VERSION_EVIDENCE.items()}
    missing = [str(p) for p in [pin_path, *version_paths.values()] if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing software-version evidence: {missing}")
    pins = read_tsv(pin_path)
    manifests = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in version_paths.items()
        if key not in {
            "ldsc_fi_log", "integrated_preflight", "latest_integrated_preflight", "post_s16_tests",
            "latest_suite", "current_canonical_suite", "bibliography_audit_latest", "quantitative_claims_audit_latest",
            "reporting_locator_audit_latest",
        }
    }
    ldsc_log = version_paths["ldsc_fi_log"].read_text(encoding="utf-8", errors="replace")
    if "3.0.1" not in ldsc_log:
        raise ValueError("LDSC run log does not record the expected 3.0.1 release")
    preflight_log = version_paths["integrated_preflight"].read_text(encoding="utf-8", errors="replace")
    required_preflight_evidence = (
        "2026-09-24 02:18–02:30 UTC",
        "PYTHON='../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python'",
        "exit_status: 0",
        "RESOURCE_MANIFEST_OK rows=350 files_verified=350",
        "unittest: Ran 94 tests; OK",
    )
    if any(marker not in preflight_log for marker in required_preflight_evidence):
        raise ValueError("Captured integrated preflight does not match the expected 94-test clean-Python snapshot")
    latest_preflight_log = version_paths["latest_integrated_preflight"].read_text(encoding="utf-8", errors="replace")
    required_latest_preflight_evidence = (
        "RESOURCE_MANIFEST_OK rows=351 files_verified=351",
        "MANIFEST_METADATA_OK rows=351 required_columns=22",
        "ANALYSIS_PLAN_LOCK_OK version=1 sleep_traits=12 multiplicity=396",
        "REVIEW_SCREENING_OK",
        "PASS empty_abstract_pmids=935 source_xml_found=935 with_abstract=0 duplicate_pmids=53",
        "PASS manuscript_claims=28 passed=28",
        "LATENT_QSNP_AUDIT_OK factors=7",
        "source_rows=32035589 significant_q_variants=0 excluded_rows=0",
        "REPORTING_LOCATORS_PASS rows=92 line_references=93 errors=0",
        "Ran 97 tests",
        "OK",
    )
    if any(marker not in latest_preflight_log for marker in required_latest_preflight_evidence):
        raise ValueError("Latest integrated preflight log does not match the 97-test, 351-resource, manuscript-claim-audited snapshot")
    post_s16_test_log = version_paths["post_s16_tests"].read_text(encoding="utf-8", errors="replace")
    if "Ran 97 tests" not in post_s16_test_log or "\nOK\n" not in post_s16_test_log:
        raise ValueError("Post-S16-update package test log does not show all 97 tests passing")
    latest_suite_log = version_paths["latest_suite"].read_text(encoding="utf-8", errors="replace")
    latest_suite_match = re.search(r"Ran (\d+) tests? in ([0-9.]+)s", latest_suite_log)
    latest_suite_capture = re.search(r"Captured UTC: (\d{4}-\d{2}-\d{2} \d{2}:\d{2})", latest_suite_log)
    if (
        "Python: 3.11.11" not in latest_suite_log
        or "Result: 169 tests, OK" not in latest_suite_log
        or latest_suite_match is None
        or latest_suite_match.group(1) != "169"
        or latest_suite_capture is None
        or "\nOK\n" not in latest_suite_log
    ):
        raise ValueError("Latest pinned Python 3.11.11 package test log does not show all 169 tests passing")
    current_suite_log = version_paths["current_canonical_suite"].read_text(encoding="utf-8", errors="replace")
    current_suite_match = re.search(r"Ran (\d+) tests? in ([0-9.]+)s", current_suite_log)
    current_suite_capture = re.search(r"Captured UTC: (\d{4}-\d{2}-\d{2} \d{2}:\d{2})", current_suite_log)
    current_suite_python = re.search(r"^Python: (.+)$", current_suite_log, re.MULTILINE)
    current_suite_head = re.search(r"^Git HEAD: ([0-9a-f]{40})$", current_suite_log, re.MULTILINE)
    current_suite_builder_sha = re.search(r"^Supplementary builder SHA-256: ([0-9a-f]{64})$", current_suite_log, re.MULTILINE)
    current_suite_tests_sha = re.search(r"^Supplementary builder tests SHA-256: ([0-9a-f]{64})$", current_suite_log, re.MULTILINE)
    hashes_match_current = (
        current_suite_builder_sha is not None
        and current_suite_tests_sha is not None
        and current_suite_builder_sha.group(1) == sha256(repo / "frailty_paper/scripts/39_build_supplementary_tables.py")
        and current_suite_tests_sha.group(1) == sha256(repo / "frailty_paper/tests/test_supplementary_tables.py")
    )
    if (
        current_suite_match is None
        or current_suite_capture is None
        or current_suite_python is None
        or current_suite_head is None
        or not current_suite_python.group(1).startswith("Python 3.13.")
        or "Exit status: 0" not in current_suite_log
        or "unittest discover -s frailty_paper/tests -v" not in current_suite_log
        or "\nOK\n" not in current_suite_log
    ):
        raise ValueError("Current canonical Python 3.13 package test log is missing interpreter, command, HEAD, or success evidence")
    bibliography_audit = json.loads(version_paths["bibliography_audit_latest"].read_text(encoding="utf-8"))
    quantitative_audit = json.loads(version_paths["quantitative_claims_audit_latest"].read_text(encoding="utf-8"))
    locator_audit = json.loads(version_paths["reporting_locator_audit_latest"].read_text(encoding="utf-8"))
    if (
        bibliography_audit.get("status") != "PASS"
        or bibliography_audit.get("missing_keys")
        or bibliography_audit.get("duplicate_keys")
        or bibliography_audit.get("uncited_keys")
        or quantitative_audit.get("status") != "PASS"
        or quantitative_audit.get("claims_checked") != 28
        or quantitative_audit.get("claims_passed") != 28
        or locator_audit.get("result") != "PASS"
        or locator_audit.get("errors")
        or locator_audit.get("checklist_rows") != 92
        or locator_audit.get("line_locators_checked") != 93
    ):
        raise ValueError("Latest manuscript citation, quantitative-claim, or reporting-locator audit did not pass")
    h2_software = manifests["fi_h2_reproduction"]["software"]
    workflow_software = manifests["workflow_validation"]["software"]
    environment_snapshot = manifests["workflow_validation"]["environment_snapshot"]
    for rel_path, expected_hash in environment_snapshot["lock_files"].items():
        lock_path = repo / rel_path
        if not lock_path.is_file() or sha256(lock_path) != expected_hash:
            raise ValueError(f"Workflow environment lock is missing or stale: {rel_path}")
    python_plot = sorted({manifests[k]["python"] for k in ("figure1", "figure2")})
    plot_keys = ("fi_forest_preview", "latent_heatmap_preview", "dimension_forest_preview")
    plot_python = sorted({manifests[k]["python"] for k in plot_keys})
    plot_mpl = sorted({manifests[k]["matplotlib"] for k in plot_keys})
    clean_venv = manifests["clean_venv_validation"]
    clean_venv_run = clean_venv.get("current_revalidation", clean_venv["test_run"])
    workflow_python_versions = set(python_plot) | set(plot_python) | {workflow_software["python"], clean_venv["python"]}
    workflow_python_pin = next((row["version_or_commit"] for row in pins if row["component"] == "python_workflow"), None)
    workflow_python_status = (
        "EXECUTED_MATCHES_PIN"
        if workflow_python_pin is not None and workflow_python_versions == {workflow_python_pin}
        else "RUN_VERSIONS_VARY_FROM_PIN"
    )
    actual = {
        "python_workflow": ("; ".join(f"{v} (Figure 1–2)" for v in python_plot) + "; " + "; ".join(f"{v} (FI/latent/dimension plots)" for v in plot_python) + f"; {workflow_software['python']} ({manifests['workflow_validation']['tests']['total']}-test workflow validation)", workflow_python_status),
        "python_ldsc": (h2_software["python"], "EXECUTED_MATCHES_PIN"),
        "numpy": (f"{h2_software['numpy']} (FI h2 runtime); {workflow_software['numpy']} (workflow validation)", "H2_RUNTIME_DIFFERS_TEST_RUNTIME_MATCHES_PIN"),
        "pandas": (f"{h2_software['pandas']} (FI h2 runtime); {workflow_software['pandas']} (workflow validation)", "H2_RUNTIME_DIFFERS_TEST_RUNTIME_MATCHES_PIN"),
        "matplotlib": ("; ".join(plot_mpl) + " (FI/latent/dimension plots); " + workflow_software["matplotlib"] + " (full-suite validation)", "EXECUTED_MATCHES_PIN" if plot_mpl == ["3.9.4"] and workflow_software["matplotlib"] == "3.9.4" else "EXECUTED_VERSION_DIFFERS_FROM_PIN"),
        "ldsc_CBIIT": ("LDSC 3.0.1 reported in run logs; source commit not verified from run manifest", "EXECUTED_RELEASE_RECORDED_SOURCE_IDENTITY_PARTIAL"),
    }
    for component in ("snakemake", "scipy", "lxml", "requests", "biopython", "beautifulsoup4", "openpyxl", "pyyaml", "tqdm"):
        pin = next((row["version_or_commit"] for row in pins if row["component"] == component), None)
        observed = workflow_software[component]
        actual[component] = (observed, "EXECUTED_MATCHES_PIN" if observed == pin else "EXECUTED_VERSION_DIFFERS_FROM_PIN")
    actual["python_tests"] = (f"{workflow_software['python']}: {manifests['workflow_validation']['tests']['passed']}/{manifests['workflow_validation']['tests']['total']} passed", "FULL_PACKAGE_TEST_SUITE_PASS")
    evidence_for = {
        "python_workflow": "; ".join(VERSION_EVIDENCE[k] for k in ("figure1", "figure2", "fi_forest_preview", "latent_heatmap_preview", "dimension_forest_preview", "workflow_validation")),
        "python_ldsc": VERSION_EVIDENCE["fi_h2_reproduction"],
        "numpy": VERSION_EVIDENCE["fi_h2_reproduction"] + "; " + VERSION_EVIDENCE["workflow_validation"],
        "pandas": VERSION_EVIDENCE["fi_h2_reproduction"] + "; " + VERSION_EVIDENCE["workflow_validation"],
        "matplotlib": "; ".join(VERSION_EVIDENCE[k] for k in ("fi_forest_preview", "latent_heatmap_preview", "dimension_forest_preview", "workflow_validation")),
        "ldsc_CBIIT": VERSION_EVIDENCE["ldsc_fi_log"] + "; " + VERSION_EVIDENCE["fi_h2_reproduction"],
    }
    for component in ("snakemake", "scipy", "lxml", "requests", "biopython", "beautifulsoup4", "openpyxl", "pyyaml", "tqdm", "python_tests"):
        evidence_for[component] = VERSION_EVIDENCE["workflow_validation"]
    rows: list[dict[str, str]] = []
    for pin in pins:
        component = pin["component"]
        observed, state = actual.get(component, ("NOT_OBSERVED_IN_CURRENT_FRAILTY_PACKAGE", "PINNED_NOT_USED_OR_NOT_EVIDENCED_HERE"))
        rows.append({
            "component": component, "pinned_version_or_commit": pin["version_or_commit"],
            "observed_version_or_commit": observed, "status": state,
            "purpose": pin["purpose"], "evidence": evidence_for.get(component, TOOL_PINS),
            "notes": "Pinned repository policy is not proof of execution; see status and cited run evidence.",
        })
    rows.extend([
        {"component": "python_tests", "pinned_version_or_commit": "Python 3.11.11; direct requirements files",
         "observed_version_or_commit": f"{workflow_software['python']}: {manifests['workflow_validation']['tests']['passed']}/{manifests['workflow_validation']['tests']['total']} passed",
         "status": "FULL_PACKAGE_TEST_SUITE_PASS", "purpose": "Full frailty package unittest suite in the isolated workflow environment",
         "evidence": VERSION_EVIDENCE["workflow_validation"], "notes": "Canonical frailty_paper/.venv remains separate; transitive dependencies are not hash-locked."},
        {"component": "python_clean_venv_tests", "pinned_version_or_commit": "Python 3.11.11; direct requirements plus resolved freeze",
         "observed_version_or_commit": f"{clean_venv['python']}: {clean_venv_run['tests']['passed']}/{clean_venv_run['tests']['total']} passed; pip check clean",
         "status": "FRESH_ISOLATED_VENV_SETUP_AND_REVALIDATION_PASS", "purpose": "Fresh setup-target build, offline setup replay, dependency check and full package suite",
         "evidence": VERSION_EVIDENCE["clean_venv_validation"] + "; " + clean_venv_run["tests"]["log"], "notes": "This setup record captures a 91/91 run. A later integrated preflight recorded 94/94 at its dated snapshot (see python_integrated_preflight_tests); separate from frailty_paper/.venv. Freeze records exact versions but not PyPI artifact hashes."},
        {"component": "python_integrated_preflight_tests", "pinned_version_or_commit": "Python 3.11.11; clean isolated environment",
         "observed_version_or_commit": "94/94 tests passed; 350/350 registered resource files verified; integrated preflight exit status 0",
         "status": "INTEGRATED_PREFLIGHT_PASS_AT_CAPTURED_SNAPSHOT", "purpose": "Combined resource, plan, overlap, reporting and package-test preflight",
         "evidence": VERSION_EVIDENCE["integrated_preflight"], "notes": "Captured 2026-09-24 02:18–02:30 UTC. The 94/94 count is the recorded run snapshot and does not attest to later worktree changes; 350/350 is the manifest size at that run."},
        {"component": "python_integrated_preflight_tests_latest", "pinned_version_or_commit": "Python 3.11.11; clean isolated environment",
         "observed_version_or_commit": "97/97 tests passed; 351/351 registered resource files verified; 935 missing-abstract PMIDs source-checked; 28 manuscript claims passed; source HEAD 6b906ba",
         "status": "INTEGRATED_PREFLIGHT_PASS_AT_CAPTURED_SNAPSHOT", "purpose": "Combined resource, plan, review-queue, bibliography, manuscript-claim, source-QC, reporting and package-test preflight",
         "evidence": VERSION_EVIDENCE["latest_integrated_preflight"], "notes": "Captured 2026-09-24 07:01–07:03 UTC from source HEAD 6b906ba in a mixed worktree. Screening remains at zero decisions and exact participant intersections remain unknown; this engineering/input-integrity run is not a scientific replay."},
        {"component": "python_post_s16_update_tests", "pinned_version_or_commit": "Python 3.11.11; clean isolated environment",
         "observed_version_or_commit": "97/97 tests passed after S16 builder and regression update; source HEAD 6b906ba",
         "status": "FULL_PACKAGE_TEST_SUITE_PASS_AT_UPDATED_SOURCE", "purpose": "Full frailty package unittest suite after S16 evidence integration",
         "evidence": VERSION_EVIDENCE["post_s16_tests"], "notes": "Standalone suite run after updating the S16 builder/test expectation; review decisions remain blank. Mixed-worktree validation, not a scientific replay."},
        {"component": "python_latest_package_suite_tests", "pinned_version_or_commit": "Python 3.11.11; isolated workflow environment",
         "observed_version_or_commit": f"Python 3.11.11: {latest_suite_match.group(1)}/{latest_suite_match.group(1)} tests passed at {latest_suite_capture.group(1)} UTC",
         "status": "FULL_PACKAGE_TEST_SUITE_PASS_AT_CAPTURED_SNAPSHOT", "purpose": "Latest full Frailty package unittest suite after current audit/figure updates",
         "evidence": VERSION_EVIDENCE["latest_suite"], "notes": "The test log does not embed a Git HEAD; this records the captured test execution and environment, not a clean-checkout scientific replay."},
        {"component": "python_current_canonical_package_suite_tests", "pinned_version_or_commit": "Python 3.13; canonical frailty_paper/.venv",
         "observed_version_or_commit": f"{current_suite_python.group(1)}: {current_suite_match.group(1)}/{current_suite_match.group(1)} tests passed at {current_suite_capture.group(1)} UTC; source HEAD {current_suite_head.group(1)[:12]}",
         "status": "FULL_PACKAGE_TEST_SUITE_PASS_AT_CAPTURED_SOURCE" if hashes_match_current else "PRIOR_SOURCE_SNAPSHOT_NOT_BOUND_TO_CURRENT_WORKTREE", "purpose": "Latest full Frailty package unittest suite in the canonical project venv",
         "evidence": VERSION_EVIDENCE["current_canonical_suite"], "notes": (f"Separate from the pinned Python 3.11.11 workflow suite. Log binds the tested supplementary-table builder/test module: {current_suite_builder_sha.group(1)[:12]} / {current_suite_tests_sha.group(1)[:12]}." if hashes_match_current else "Separate from the pinned Python 3.11.11 workflow suite. Captured test success is retained as historical evidence; its builder/test source hashes do not match the current worktree.")},
        {"component": "current_manuscript_validation_audits", "pinned_version_or_commit": "Python 3.11.11; isolated workflow environment",
         "observed_version_or_commit": (f"Bibliography PASS ({bibliography_audit['citation_uses']} citation uses/{bibliography_audit['bibliography_entries']} entries); "
                                         f"quantitative claims PASS ({quantitative_audit['claims_passed']}/{quantitative_audit['claims_checked']}); "
                                         f"reporting locators PASS ({locator_audit['checklist_rows']} rows/{locator_audit['line_locators_checked']} references)"),
         "status": "CURRENT_MANUSCRIPT_AUDITS_PASS", "purpose": "Current citation-key, quantitative-claim, and checklist-locator validation",
         "evidence": "; ".join(VERSION_EVIDENCE[k] for k in ("bibliography_audit_latest", "quantitative_claims_audit_latest", "reporting_locator_audit_latest")),
         "notes": "Range and claim checks do not establish item-level reporting compliance; the manuscript remains provisional."},
        {"component": "workflow_environment_snapshot",
         "pinned_version_or_commit": "osx-arm64 Conda explicit lock + exact PyPI package versions",
         "observed_version_or_commit": f"{environment_snapshot['conda_package_count']} Conda packages SHA256-pinned; {environment_snapshot['pypi_package_count']} PyPI packages version-pinned",
         "status": "PLATFORM_SCOPED_SNAPSHOT_TESTED",
         "purpose": "Reconstruct the tested Python 3.11.11 workflow environment on Apple silicon",
         "evidence": VERSION_EVIDENCE["workflow_validation"],
         "notes": "Conda package artifacts have SHA-256 values; PyPI artifact hashes are not recorded. R/system binaries and the canonical .venv are outside this snapshot."},
        {"component": "rsvg-convert", "pinned_version_or_commit": "NOT_IN_TOOL_PIN_FILE",
         "observed_version_or_commit": (
             f"Figures 1–2: {manifests['figure2']['rsvg_convert_version'].splitlines()[0]}; "
             f"Figure 4: {manifests['figure4_current']['renderer'].splitlines()[0]}"
         ),
         "status": "EXECUTED_FOR_FIGURES1_2_4", "purpose": "SVG to vector PDF and raster PNG rendering",
         "evidence": "; ".join((VERSION_EVIDENCE["figure1"], VERSION_EVIDENCE["figure2"], VERSION_EVIDENCE["figure4_current"])),
         "notes": "Renderer reports librsvg 2.62.3; cairo/pango/harfbuzz/fontconfig versions are retained in figure provenance."},
        {"component": "scipy_runtime_ldsc", "pinned_version_or_commit": "NOT_LISTED_IN_TOOL_PIN_FILE",
         "observed_version_or_commit": h2_software["scipy"], "status": "FI_H2_REPRODUCTION_RUNTIME",
         "purpose": "Numerical runtime for the documented FI h2 reproduction",
         "evidence": VERSION_EVIDENCE["fi_h2_reproduction"], "notes": "Exact version recorded in the reproduction manifest."},
        {"component": "bitarray_runtime_ldsc", "pinned_version_or_commit": "NOT_LISTED_IN_TOOL_PIN_FILE",
         "observed_version_or_commit": h2_software["bitarray"], "status": "FI_H2_REPRODUCTION_RUNTIME",
         "purpose": "Runtime dependency for the documented FI h2 reproduction",
         "evidence": VERSION_EVIDENCE["fi_h2_reproduction"], "notes": "Exact version recorded in the reproduction manifest."},
    ])
    return rows


def build_tables(repo: Path) -> tuple[dict[str, list[dict[str, str]]], dict[str, int]]:
    paths = {key: repo / rel for key, rel in INPUTS.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing locked table input(s): {missing}")
    src = {key: read_tsv(path) for key, path in paths.items()}

    s1 = [r for r in src["resource_manifest"] if r["resource_type"] == "GWAS summary statistics"]
    if len(s1) != 33 or len({r["resource_id"] for r in s1}) != 33:
        raise ValueError("Registered GWAS metadata family changed; expected 33 unique resources")

    s4 = src["h2_master"]
    if len(s4) != 20 or len({r["trait_id"] for r in s4}) != 20:
        raise ValueError("Available h2 master family changed; expected 20 unique traits")

    fi = src["fi_rg"]
    latent = src["latent_rg"]
    aging = src["aging_rg"]
    fi_overlap = [r for r in src["fi_overlap"] if r["frailty_endpoint"] == "frailty"]
    aging_overlap = src["aging_overlap"]
    if (len(fi), len(latent), len(aging)) != (12, 84, 72):
        raise ValueError("Global-rg source family sizes changed; expected 12 FI, 84 latent, 72 aging")
    if len(fi_overlap) != 12 or len(aging_overlap) != 72:
        raise ValueError("Overlap crosswalk no longer covers the FI and aging pairs exactly")
    fi_overlap_by_sleep = {r["sleep_trait"]: r for r in fi_overlap}
    aging_overlap_by_pair = {r["pair_id"]: r for r in aging_overlap}
    if set(fi_overlap_by_sleep) != {r["sleep_trait"] for r in fi}:
        raise ValueError("FI overlap crosswalk does not map one-to-one to the frozen FI panel")
    if set(aging_overlap_by_pair) != {r["pair_id"] for r in aging}:
        raise ValueError("Aging overlap crosswalk does not map one-to-one to the frozen aging family")

    s5: list[dict[str, str]] = []
    for r in fi:
        ov = fi_overlap_by_sleep[r["sleep_trait"]]
        s5.append({
            "analysis_family": "frozen_atlas_sleep_x_FI",
            "analysis_class": "PRIMARY_FROZEN_ATLAS_REUSE",
            "pair_id": r["pair_id"], "sleep_trait": r["sleep_trait"], "outcome_trait": r["non_sleep_trait"],
            "rg": r["global_rg"], "se": r["global_rg_se"], "z": r["global_rg_z"],
            "p_value": r["global_rg_p"], "q_value": r["global_rg_fdr_all_396"],
            "q_family": "frozen atlas all-396 BH", "family_denominator": "396",
            "direction": r["effect_direction"], "cohort_overlap_status": ov["cohort_overlap_status"],
            "exact_participant_overlap": ov["exact_participant_overlap"],
            "replication_status": "NO_INDEPENDENT_PAIRWISE_REPLICATION_ESTABLISHED",
            "interpretation_status": r["interpretation_status"],
            "claim_limit": "Frozen FI atlas reuse; global rg only; not independent replication or causality",
            "source_table": INPUTS["fi_rg"],
        })
    for r in latent:
        pair = (r["sleep_trait"], r["disease_trait"])
        if pair == ("insomnia", "frailty_general"):
            indicator_overlap_status = "DIRECT_CONSTRUCT_INDICATOR_MATCH_PART_WHOLE_DEPENDENCE"
            indicator_overlap_note = (
                "Insomnia is one of the 30 deficit GWAS used to construct the general factor. "
                "The project and factor studies use different GWAS releases, but UK Biobank "
                "overlap is expected and exact participant intersection is unknown; this is "
                "not independent validation."
            )
        elif r["sleep_trait"] == "sleepiness":
            indicator_overlap_status = "RELATED_SLEEP_INDICATOR_IN_SOURCE_MODEL_FACTOR_LOADING_UNRESOLVED"
            indicator_overlap_note = (
                "The primary factor model includes tiredness/lethargy (TIR), related but not "
                "identical to the project daytime-sleepiness phenotype. Exact factor loading "
                "was not resolved in this audit; do not treat as a direct indicator match."
            )
        else:
            indicator_overlap_status = "NO_DIRECT_MATCH_IDENTIFIED_IN_MAIN_SOURCE_INDICATOR_LIST"
            indicator_overlap_note = (
                "This label was not an exact match to an indicator listed in the primary "
                "article's 30-deficit figure. This does not resolve participant-level overlap."
            )
        s5.append({
            "analysis_family": "secondary_sleep_x_latent_frailty",
            "analysis_class": "SENSITIVITY_ONLY",
            "pair_id": f'{r["sleep_trait"]}__{r["disease_trait"]}',
            "sleep_trait": r["sleep_trait"], "outcome_trait": r["disease_trait"],
            "rg": r["rg"], "se": r["se"], "z": r["z"], "p_value": r["p_value"], "q_value": r["q_value"],
            "q_family": "fixed secondary 84-pair BH", "family_denominator": r["family_denominator"],
            "direction": "POSITIVE" if float(r["rg"]) > 0 else "NEGATIVE" if float(r["rg"]) < 0 else "ZERO",
            "cohort_overlap_status": r["cohort_overlap_status"], "exact_participant_overlap": "UNKNOWN",
            "replication_status": "NO_INDEPENDENT_PAIRWISE_REPLICATION_ESTABLISHED",
            "interpretation_status": r["interpretation_status"], "claim_limit": r["claim_limit"],
            "source_table": INPUTS["latent_rg"],
        })
    for r in aging:
        ov = aging_overlap_by_pair[r["pair_id"]]
        s5.append({
            "analysis_family": "frozen_atlas_sleep_x_aging_context",
            "analysis_class": "READ_ONLY_CONTEXT_NOT_REPLICATION",
            "pair_id": r["pair_id"], "sleep_trait": r["sleep_trait"], "outcome_trait": r["non_sleep_trait"],
            "rg": r["global_rg"], "se": r["global_rg_se"], "z": r["global_rg_z"],
            "p_value": r["global_rg_p"], "q_value": r["global_rg_fdr_all_396"],
            "q_family": "frozen atlas all-396 BH", "family_denominator": r["fdr_denominator"].split()[0],
            "direction": r["effect_direction"], "cohort_overlap_status": ov["overall_cohort_status"],
            "exact_participant_overlap": ov["exact_participant_intersection"],
            "replication_status": "READ_ONLY_ATLAS_REUSE_NOT_REPLICATION",
            "interpretation_status": r["reuse_status"], "claim_limit": r["claim_limit"],
            "source_table": INPUTS["aging_rg"],
        })
    if len(s5) != 168:
        raise ValueError(f"Combined available global-rg table has {len(s5)} rows, expected 168")

    replication = src["replication_audit"]
    correction = src["fi_correction"]
    latent_correction = src["latent_correction"]
    if not replication or len(correction) != 12 or len(latent_correction) != 84:
        raise ValueError("Replication audit or correction-sensitivity family is incomplete")
    s15 = []
    for r in correction:
        s15.append({
            "analysis_family": "primary_sleep_x_FI",
            "pair_id": r["pair_id"], "sleep_trait": r["sleep_trait"],
            "outcome_trait": r["non_sleep_trait"], "rg": r["global_rg"], "se": r["global_rg_se"],
            "raw_p": r["global_rg_p_frozen"], "bh_q": r["bh_q_all_396_frozen"],
            "bh_family_denominator": "396", "bonferroni_family_size": r["bonferroni_family_size"],
            "bonferroni_p": r["bonferroni_p_all_396"],
            "bh_significant_at_0.05": r["bh_all_396_significant_at_0.05"],
            "bonferroni_significant_at_0.05": r["bonferroni_all_396_significant_at_0.05"],
            "cohort_overlap_status": "EXPECTED_OR_POSSIBLE; exact participant overlap UNKNOWN",
            "exact_participant_overlap": "UNKNOWN",
            "comparison_interpretation": "Same frozen P values; BH remains primary; all-396 Bonferroni sensitivity.",
            "source_table": INPUTS["fi_correction"],
        })
    for r in latent_correction:
        s15.append({
            "analysis_family": r["analysis_family"],
            "pair_id": f'{r["sleep_trait"]}__{r["frailty_factor"]}',
            "sleep_trait": r["sleep_trait"], "outcome_trait": r["frailty_factor"],
            "rg": r["rg"], "se": r["se"], "raw_p": r["raw_p"],
            "bh_q": r["bh_q_secondary_84"], "bh_family_denominator": "84",
            "bonferroni_family_size": r["bonferroni_family_size"],
            "bonferroni_p": r["bonferroni_p_secondary_84"],
            "bh_significant_at_0.05": r["bh_significant_at_0.05"],
            "bonferroni_significant_at_0.05": r["bonferroni_significant_at_0.05"],
            "cohort_overlap_status": r["cohort_overlap_status"],
            "exact_participant_overlap": r["exact_participant_overlap"],
            "comparison_interpretation": r["comparison_interpretation"] + "; BH remains primary; same-family 84-pair Bonferroni sensitivity.",
            "source_table": INPUTS["latent_correction"],
        })

    s7 = []
    for r in latent:
        row = dict(r)
        pair = (r["sleep_trait"], r["disease_trait"])
        if pair == ("insomnia", "frailty_general"):
            row["indicator_overlap_status"] = "DIRECT_CONSTRUCT_INDICATOR_MATCH_PART_WHOLE_DEPENDENCE"
            row["indicator_overlap_note"] = (
                "Insomnia is one of the 30 deficit GWAS used to construct the general factor. "
                "The project and factor studies use different GWAS releases, but UK Biobank "
                "overlap is expected and exact participant intersection is unknown; this is "
                "not independent validation."
            )
        elif r["sleep_trait"] == "sleepiness":
            row["indicator_overlap_status"] = "RELATED_SLEEP_INDICATOR_IN_SOURCE_MODEL_FACTOR_LOADING_UNRESOLVED"
            row["indicator_overlap_note"] = (
                "The primary factor model includes tiredness/lethargy (TIR), related but not "
                "identical to the project daytime-sleepiness phenotype. Exact factor loading "
                "was not resolved in this audit; do not treat as a direct indicator match."
            )
        else:
            row["indicator_overlap_status"] = "NO_DIRECT_MATCH_IDENTIFIED_IN_MAIN_SOURCE_INDICATOR_LIST"
            row["indicator_overlap_note"] = (
                "This label was not an exact match to an indicator listed in the primary "
                "article's 30-deficit figure. This does not resolve participant-level overlap."
            )
        s7.append(row)

    tables = {
        TABLE_FILES["S1"]: s1,
        TABLE_FILES["S4"]: s4,
        TABLE_FILES["S5"]: s5,
        TABLE_FILES["S7"]: s7,
        TABLE_FILES["S8"]: replication,
        TABLE_FILES["S15"]: s15,
        TABLE_FILES["S15_matrix"]: src["sensitivity_summary"],
        TABLE_FILES["S16"]: software_version_table(repo),
    }
    return tables, {key: len(value) for key, value in tables.items()}


def table_status(counts: dict[str, int]) -> list[dict[str, str]]:
    available = {
        "S1": ("PARTIAL_METADATA_TABLE", "All 33 registered GWAS metadata rows; source eligibility differs by phenotype."),
        "S4": ("PARTIAL_H2_TABLE", "20 current sleep/FI/latent h2 rows; physical, Fried and HFRS rows unavailable."),
        "S5": ("PARTIAL_GLOBAL_RG_TABLE", "168 FI, latent and read-only aging-context estimates; no eligible physical/Fried/HFRS estimates."),
        "S7": ("SENSITIVITY_ONLY", "84 secondary latent-factor rg estimates; direct insomnia/general-factor part-whole dependence is flagged; exact participant overlap remains unknown."),
        "S8": ("PARTIAL_REPLICATION_FEASIBILITY", "Audited candidate resources; no eligible independent pairwise sleep-frailty rg replication established."),
        "S15": ("PARTIAL_CORRECTION_SENSITIVITY", "The numeric correction table reports primary sleep × FI all-396 Bonferroni and latent-factor same-family 84-pair Bonferroni. The companion conclusion matrix distinguishes completed, partial, unresolved, blocked, and not-justified sensitivity domains; unavailable tests are not null results."),
        "S16": ("PARTIAL_VERSION_INVENTORY", "Includes locked versions and captured analysis, validation, and current manuscript-audit executions; does not cover every future tool run or the final full scientific replay."),
    }
    descriptions = {
        "S1": "GWAS metadata", "S2": "Systematic-review studies", "S3": "Risk-of-bias results",
        "S4": "SNP heritability", "S5": "Global genetic correlation", "S6": "Physical-component results",
        "S7": "Latent-factor results", "S8": "Replication", "S9": "LAVA", "S10": "PLACO/shared loci",
        "S11": "Fine-mapping", "S12": "Trait-trait colocalization", "S13": "Molecular QTL evidence",
        "S14": "Cell-type evidence", "S15": "Sensitivity analyses", "S16": "Software/resources/versions",
    }
    missing_state = {
        "S2": ("NOT_AVAILABLE_UNSCREENED", "No included studies; 56,117 records remain unscreened and licensed exports are absent."),
        "S3": ("NOT_AVAILABLE_NO_APPRAISALS", "No included studies have been appraised."),
        "S6": ("BLOCKED_SOURCE_PROVENANCE", "Five candidate files are structurally readable but lack eligible source/build/effect provenance."),
        "S9": ("NOT_JUSTIFIED", "No eligible independent replicated frailty pair has passed upstream gates."),
        "S10": ("NOT_JUSTIFIED", "No eligible frailty-specific LAVA-prioritized loci."),
        "S11": ("NOT_JUSTIFIED", "No eligible frailty-specific shared loci."),
        "S12": ("NOT_JUSTIFIED", "No eligible frailty-specific locus inputs."),
        "S13": ("NOT_JUSTIFIED", "No eligible frailty-specific shared loci for molecular-QTL analysis."),
        "S14": ("NOT_JUSTIFIED", "No eligible candidate loci/genes for cell-type analysis."),
    }
    rows = []
    for key in [f"S{i}" for i in range(1, 17)]:
        if key in available:
            status, note = available[key]
            file = TABLE_FILES[key]
            n = str(counts[file])
            source = file
            if key == "S15":
                matrix_file = TABLE_FILES["S15_matrix"]
                n = f'{counts[file]} correction rows + {counts[matrix_file]} conclusion-by-sensitivity rows'
                file = f'{file}; {matrix_file}'
                source = file
        else:
            status, note = missing_state[key]
            file = "NOT_GENERATED"
            n = "NA"
            source = ""
        rows.append({"table": key, "title": descriptions[key], "status": status,
                     "data_rows": n, "file": file, "source": source, "limitation": note})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--replace-outputs", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    outdir = (args.outdir or repo / "frailty_paper/paper/supplementary_tables").resolve()
    if repo not in outdir.parents:
        raise SystemExit("Output directory must remain under the repository")
    tables, counts = build_tables(repo)
    tables["supplementary_table_status.tsv"] = table_status(counts)
    outdir.mkdir(parents=True, exist_ok=True)
    output_paths: dict[str, Path] = {}
    for name, rows in tables.items():
        columns = S15_COLUMNS if name == TABLE_FILES["S15"] else S7_COLUMNS if name == TABLE_FILES["S7"] else None
        data = tsv_bytes(rows, columns)
        path = outdir / name
        if path.exists() and path.read_bytes() != data and not args.replace_outputs:
            raise SystemExit(f"Refusing to overwrite non-identical output: {path}")
        path.write_bytes(data)
        output_paths[name] = path

    input_paths = {name: repo / rel for name, rel in INPUTS.items()}
    version_inputs = {"tool_version_pins": repo / TOOL_PINS}
    version_inputs.update({f"version_evidence_{key}": repo / rel for key, rel in VERSION_EVIDENCE.items()})
    all_input_paths = {name: path for name, path in input_paths.items() if name != "resource_manifest"}
    all_input_paths.update(version_inputs)
    provenance_inputs = {
        name: {"path": path.relative_to(repo).as_posix(), "sha256": sha256(path)}
        for name, path in all_input_paths.items()
    }
    manifest_path = input_paths["resource_manifest"]
    manifest_hash, manifest_rows = gwas_manifest_scope(manifest_path)
    provenance_inputs["resource_manifest_gwas_summary_statistics"] = {
        "path": manifest_path.relative_to(repo).as_posix(),
        "sha256": manifest_hash,
        "scope": "Rows with resource_type == 'GWAS summary statistics'; exact subset used to build table_s1_gwas_metadata.tsv",
        "row_count": manifest_rows,
    }
    provenance = {
        "schema_version": "frailty_supplementary_tables.v1",
        "script": Path(__file__).resolve().relative_to(repo).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": provenance_inputs,
        "outputs": {path.name: {"path": path.relative_to(repo).as_posix(), "sha256": sha256(path)}
                    for path in output_paths.values()},
        "row_counts": counts,
        "global_rg_families": {"primary_FI": 12, "secondary_latent_sensitivity": 84,
                                "read_only_aging_context": 72, "combined": 168},
        "limitations": [
            "The combined global-rg table keeps distinct correction families and analysis classes explicit.",
            "Aging-context estimates are read-only frozen-atlas reuse, not replication.",
            "FFS summary statistics are acquired and registered, with GRCh37/hg19 coordinates empirically verified; the signed-effect convention is inferred from official FUMA/BOLT documentation, while exact participant overlap remains unresolved and limits cross-trait replication claims. Physical-component provenance and full HFRS statistics remain incomplete, and frailty-specific downstream locus/molecular tables are not justified.",
        ],
    }
    prov_path = outdir / "supplementary_tables.provenance.json"
    prov_data = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode()
    if prov_path.exists() and prov_path.read_bytes() != prov_data and not args.replace_outputs:
        raise SystemExit(f"Refusing to overwrite non-identical output: {prov_path}")
    prov_path.write_bytes(prov_data)
    print(json.dumps({"row_counts": counts, "status_rows": 16,
                      "outputs": {p.name: sha256(p) for p in output_paths.values()},
                      "provenance": prov_path.relative_to(repo).as_posix()}, indent=2))


if __name__ == "__main__":
    main()
