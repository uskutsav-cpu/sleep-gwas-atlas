#!/usr/bin/env python3
"""Build a compact, path-free candidate release from frozen V4 result tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[3]
PACKAGE = Path(__file__).resolve().parents[1]
V4 = ROOT / "sleep_unified_research_v4"
SOURCE_DICT = ROOT / "sleep_integrated_discovery_methods_v1/manuscript/sleep_construct_source_dictionary_v1.tsv"

RG_FIELDS = [
    "analysis_family", "family_size", "sleep_trait", "sleep_source_id",
    "outcome_trait", "discovery_outcome_trait_id", "outcome_source_id", "rg", "se", "z", "p",
    "ci_lower_95", "ci_upper_95", "adjusted_p_or_q", "adjustment_method",
    "threshold_pass", "result_class", "reproduction_status",
    "historical_full_precision_available", "pairwise_h2_z_ge_4_both",
    "pairwise_h2_intercept_le_1_2_both", "fully_independent_two_trait_replication",
    "interpretation",
]
H2_FIELDS = [
    "analysis_family", "trait_id", "source_id", "h2_observed", "se_observed",
    "h2_reported", "se_reported", "scale", "h2_z", "intercept", "intercept_se",
    "input_snp_count", "regression_snp_count", "h2_z_ge_4", "intercept_le_1_2",
    "sample_prevalence_argument", "population_prevalence_argument", "liability_factor",
    "interpretation", "historical_full_precision_available",
]
SENS_FIELDS = [
    "job_id", "kind", "sleep_trait", "outcome_trait", "status",
    "baseline_estimate", "baseline_se", "baseline_ci95_lower", "baseline_ci95_upper",
    "sensitivity_estimate", "sensitivity_se", "sensitivity_ci95_lower", "sensitivity_ci95_upper",
    "descriptive_delta", "delta_uncertainty", "baseline_p", "sensitivity_p",
    "difference_p_value", "cross_fit_delete_alignment_certified",
]
SLEEP_FIELDS = [
    "trait_id", "source_id", "source_release", "construct", "wording_device_and_coding",
    "assessment_window", "source_effect_units", "effect_allele", "effect_column",
    "SE_column", "P_column", "positive_direction", "N", "cases", "controls",
    "cohort", "ancestry", "build", "reliability_ascertainment", "sample_overlap",
    "primary_source_urls", "missing_or_unverified_fields",
]
OUTCOME_FIELDS = [
    "analysis_family", "trait_id", "label", "domain", "source_id", "source_release_label",
    "phenotype_definition", "ancestry", "build", "sample_size", "cases", "controls",
    "phenotype_type", "PMID", "DOI", "source_note",
]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], data: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(data)


def boolean(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def build_rg(stage: str, path: Path, size: int, sleep_sources: dict[str, str], outcome_sources: dict[str, str], out: Path) -> list[dict[str, str]]:
    source = rows(path)
    if len(source) != size:
        raise ValueError(f"{path.name}: expected {size} rows, found {len(source)}")
    result = []
    for row in source:
        sleep = row["sleep_trait"]
        outcome = row["outcome_trait"]
        if sleep not in sleep_sources:
            raise ValueError(f"No sleep source mapping for {sleep}")
        if stage != "validation" and outcome not in outcome_sources:
            raise ValueError(f"No outcome source mapping for {outcome}")
        if stage == "validation":
            adjusted = ""
            method = "Bonferroni_0.05_over_217"
            passes = row["native_p_pass_original_217_bonferroni"]
            result_class = row["historical_classification"]
            independent = row["independent_two_trait_replication"]
            interpretation = row["current_interpretation"]
            outcome_source = row.get("original_replication_source_id", "")
            discovery_outcome_id = row.get("original_extension_trait_id", "")
            if not discovery_outcome_id:
                raise ValueError(f"No discovery extension trait ID for validation pair {sleep}/{outcome}")
        else:
            adjusted = row["frozen_fdr_original_complete_family"]
            method = f"frozen_BH_{size}"
            passes = row["frozen_fdr_pass_0_05"]
            result_class = "FDR_PASS" if boolean(passes) else "BELOW_FROZEN_FDR_THRESHOLD"
            independent = "False"
            interpretation = "Marginal genetic correlation; not a causal or between-phenotype contrast."
            outcome_source = outcome_sources[outcome]
            discovery_outcome_id = outcome if stage == "extension" else ""
        result.append({
            "analysis_family": stage,
            "family_size": str(217 if stage == "validation" else size),
            "sleep_trait": sleep,
            "sleep_source_id": sleep_sources[sleep],
            "outcome_trait": outcome,
            "discovery_outcome_trait_id": discovery_outcome_id,
            "outcome_source_id": outcome_source,
            "rg": row["rg"], "se": row["se"], "z": row["z"], "p": row["p"],
            "ci_lower_95": row["ci_lower_95"], "ci_upper_95": row["ci_upper_95"],
            "adjusted_p_or_q": adjusted,
            "adjustment_method": method,
            "threshold_pass": passes,
            "result_class": result_class,
            "reproduction_status": row["reproduction_status"],
            "historical_full_precision_available": row["historical_full_precision_available"],
            "pairwise_h2_z_ge_4_both": row.get("pairwise_h2_z_ge_4_both", ""),
            "pairwise_h2_intercept_le_1_2_both": row.get("pairwise_h2_intercept_le_1_2_both", ""),
            "fully_independent_two_trait_replication": independent,
            "interpretation": interpretation,
        })
    write_tsv(out, RG_FIELDS, result)
    return result


def build_h2(path: Path, family: str, source_map: dict[str, str]) -> list[dict[str, str]]:
    result = []
    for row in rows(path):
        trait = row["trait_id"]
        source = source_map.get(trait, row.get("original_replication_source_id", ""))
        if not source:
            raise ValueError(f"No source ID for heritability row {trait}")
        result.append({
            "analysis_family": family,
            "trait_id": trait,
            "source_id": source,
            "h2_observed": row["h2_observed_full_precision"],
            "se_observed": row["h2_observed_se_full_precision"],
            "h2_reported": row["reported_h2_full_precision"],
            "se_reported": row["reported_h2_se_full_precision"],
            "scale": row["scale"],
            "h2_z": row["h2_z_full_precision"],
            "intercept": row["intercept_full_precision"],
            "intercept_se": row["intercept_se_full_precision"],
            "input_snp_count": row["input_snp_count"],
            "regression_snp_count": row["regression_snp_count"],
            "h2_z_ge_4": row["h2_z_ge_4_full_precision"],
            "intercept_le_1_2": row["intercept_le_1_2_full_precision"],
            "sample_prevalence_argument": row["sample_prevalence_argument"],
            "population_prevalence_argument": row["population_prevalence_argument"],
            "liability_factor": row["liability_factor"],
            "interpretation": row["observed_fit_interpretation"],
            "historical_full_precision_available": row["historical_full_precision_available"],
        })
    return result


def build_sensitivity(path: Path) -> list[dict[str, str]]:
    result = []
    for row in rows(path):
        if row["kind"] == "rg":
            base = ("baseline_rg_ratio", "baseline_rg_se", "baseline_CI95_lower", "baseline_CI95_upper")
            current = ("sensitivity_rg_ratio", "sensitivity_rg_se", "sensitivity_CI95_lower", "sensitivity_CI95_upper")
            delta = row["delta_rg"]
            bp, sp = row["baseline_p"], row["sensitivity_p"]
        else:
            base = ("baseline_tot", "baseline_tot_se", "baseline_CI95_lower", "baseline_CI95_upper")
            current = ("sensitivity_tot", "sensitivity_tot_se", "sensitivity_CI95_lower", "sensitivity_CI95_upper")
            delta = ""
            bp, sp = "", ""
        result.append({
            "job_id": row["job_id"], "kind": row["kind"],
            "sleep_trait": row["sleep_trait"], "outcome_trait": row["outcome_trait"],
            "status": row["status"], "baseline_estimate": row[base[0]], "baseline_se": row[base[1]],
            "baseline_ci95_lower": row[base[2]], "baseline_ci95_upper": row[base[3]],
            "sensitivity_estimate": row[current[0]], "sensitivity_se": row[current[1]],
            "sensitivity_ci95_lower": row[current[2]], "sensitivity_ci95_upper": row[current[3]],
            "descriptive_delta": delta,
            "delta_uncertainty": "NOT_CALIBRATED; do not interpret as a difference test or interval",
            "baseline_p": bp, "sensitivity_p": sp, "difference_p_value": "",
            "cross_fit_delete_alignment_certified": row["cross_fit_delete_alignment_certified"],
        })
    return result


def build_outcome_catalog(core_path: Path, extension_path: Path, validation_h2_path: Path) -> list[dict[str, str]]:
    catalog = []
    for row in rows(core_path):
        catalog.append({
            "analysis_family": "core", "trait_id": row["trait_id"], "label": row["label"],
            "domain": row["domain"], "source_id": row["source_id"],
            "source_release_label": row["dataset_version"], "phenotype_definition": row["phenotype_definition"],
            "ancestry": row["ancestry"], "build": row["build"], "sample_size": row["n_total"],
            "cases": row["ncase"], "controls": row["ncontrol"], "phenotype_type": row["type"],
            "PMID": row["pmid"], "DOI": row["doi"], "source_note": row["source_note"],
        })
    for row in rows(extension_path):
        catalog.append({
            "analysis_family": "extension", "trait_id": row["extension_trait_id"],
            "label": row["phenotype_name"], "domain": row["phenotype_domain"],
            "source_id": row["study_accession"], "source_release_label": row["source"],
            "phenotype_definition": row["phenotype_definition"], "ancestry": row["ancestry"],
            "build": row["build"], "sample_size": row["sample_size"], "cases": row["cases"],
            "controls": row["controls"], "phenotype_type": row["binary_or_continuous"],
            "PMID": row["PMID"], "DOI": row["DOI"],
            "source_note": "Source release label is copied from the frozen panel; it is not independent proof of a separate regenerated release.",
        })
    validation_seen = set()
    for row in rows(validation_h2_path):
        source_id = row["original_replication_source_id"]
        if source_id in validation_seen:
            continue
        validation_seen.add(source_id)
        catalog.append({
            "analysis_family": "validation", "trait_id": source_id,
            "label": row["original_phenotype_definition"],
            "domain": "external validation outcome", "source_id": source_id,
            "source_release_label": row["original_study_accession"],
            "phenotype_definition": row["original_phenotype_definition"],
            "ancestry": row["original_ancestry"], "build": "", "sample_size": "",
            # The legacy validation field records the h2 scale (e.g. observed
            # or liability), not the phenotype's binary/continuous type.
            # Keep phenotype type unknown instead of conflating these concepts.
            "cases": "", "controls": "", "phenotype_type": "UNKNOWN",
            "PMID": "", "DOI": "",
            "source_note": "Sample counts/build/rights are not included because the frozen validation metadata did not establish them in this release crosswalk.",
        })
    keys = [(row["analysis_family"], row["trait_id"]) for row in catalog]
    if len(catalog) != 158 or len(keys) != len(set(keys)):
        raise ValueError(f"Outcome catalog must have 158 unique family/trait rows; found {len(catalog)}")
    return catalog


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=PACKAGE / "release_candidate")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sleep_meta = rows(SOURCE_DICT)
    sleep_sources = {row["trait_id"]: row["source_id"] for row in sleep_meta}
    core = rows(ROOT / "config/analysis_panel.tsv")
    core_sources = {row["trait_id"]: row["source_id"] for row in core}
    ext = rows(ROOT / "discovery_extension/config/candidate_traits.tsv")
    ext_sources = {row["extension_trait_id"]: row["study_accession"] for row in ext}
    outcome_sources = {**core_sources, **ext_sources}

    groups = [
        ("core_rg.tsv", "core", "core_native_full_precision_rg_v4.tsv", 396, outcome_sources),
        ("extension_rg.tsv", "extension", "extension_native_full_precision_rg_v4.tsv", 1200, outcome_sources),
        ("external_validation_rg.tsv", "validation", "validation_native_full_precision_rg_v4.tsv", 41, {}),
    ]
    source_files = [SOURCE_DICT, ROOT / "config/analysis_panel.tsv", ROOT / "discovery_extension/config/candidate_traits.tsv"]
    manifest_files: dict[str, dict[str, object]] = {}
    rg_counts = {}
    for name, stage, table, count, maps in groups:
        source = V4 / "tables" / table
        source_files.append(source)
        built = build_rg(stage, source, count, sleep_sources, maps, out / name)
        rg_counts[name] = len(built)
        manifest_files[name] = {"rows": len(built), "sha256": sha256(out / name)}

    h2_outputs = []
    for family, filename, mapping in [
        ("core", "core_native_full_precision_h2_v4.tsv", core_sources),
        ("extension", "extension_native_full_precision_h2_v4.tsv", ext_sources),
        ("validation", "validation_native_full_precision_h2_v4.tsv", {}),
    ]:
        src = V4 / "tables" / filename
        source_files.append(src)
        h2_outputs.extend(build_h2(src, family, mapping))
    outcome_source = V4 / "tables/validation_native_full_precision_h2_v4.tsv"
    source_files.append(outcome_source)
    outcome_catalog = build_outcome_catalog(
        ROOT / "config/analysis_panel.tsv",
        ROOT / "discovery_extension/config/candidate_traits.tsv",
        outcome_source,
    )
    write_tsv(out / "outcome_phenotypes.tsv", OUTCOME_FIELDS, outcome_catalog)
    manifest_files["outcome_phenotypes.tsv"] = {"rows": len(outcome_catalog), "sha256": sha256(out / "outcome_phenotypes.tsv")}
    if len(h2_outputs) != 158:
        raise ValueError(f"Expected 158 heritability estimates, found {len(h2_outputs)}")
    write_tsv(out / "heritability.tsv", H2_FIELDS, h2_outputs)
    manifest_files["heritability.tsv"] = {"rows": len(h2_outputs), "sha256": sha256(out / "heritability.tsv")}

    sens_source = V4 / "tables/independent_sensitivity_62_full_precision_with_CI_v4_6.tsv"
    source_files.append(sens_source)
    sensitivities = build_sensitivity(sens_source)
    if len(sensitivities) != 62:
        raise ValueError(f"Expected 62 sensitivity rows, found {len(sensitivities)}")
    write_tsv(out / "sensitivities.tsv", SENS_FIELDS, sensitivities)
    manifest_files["sensitivities.tsv"] = {"rows": len(sensitivities), "sha256": sha256(out / "sensitivities.tsv")}

    public_sleep_meta = [{k: row.get(k, "") for k in SLEEP_FIELDS} for row in sleep_meta]
    if len(public_sleep_meta) != 12:
        raise ValueError(f"Expected 12 sleep construct rows, found {len(public_sleep_meta)}")
    write_tsv(out / "sleep_constructs.tsv", SLEEP_FIELDS, public_sleep_meta)
    manifest_files["sleep_constructs.tsv"] = {"rows": len(public_sleep_meta), "sha256": sha256(out / "sleep_constructs.tsv")}
    source_files.append(V4 / "tables/independent_validation_admissibility217_v4_1.tsv")
    classes = {}
    for row in rows(source_files[-1]):
        cls = row["current_class"]
        classes[cls] = classes.get(cls, 0) + 1
    class_rows = [{"classification": key, "n_pairs": str(classes[key])} for key in sorted(classes)]
    write_tsv(out / "validation_217_class_counts.tsv", ["classification", "n_pairs"], class_rows)
    manifest_files["validation_217_class_counts.tsv"] = {"rows": len(class_rows), "sha256": sha256(out / "validation_217_class_counts.tsv")}

    for filename in ("SCHEMA.md", "schema.json"):
        source = PACKAGE / "release_candidate" / filename
        target = out / filename
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
        manifest_files[filename] = {"sha256": sha256(target), "file_role": "schema"}

    source_records = []
    for source in dict.fromkeys(source_files):
        source_records.append({"path": source.relative_to(ROOT).as_posix(), "sha256": sha256(source)})
    manifest = {
        "schema_version": "sleep-atlas-candidate-release/1.0",
        "release_status": "CANDIDATE_FOR_RIGHTS_AND_METADATA_REVIEW",
        "builder": {
            "path": Path(__file__).resolve().relative_to(ROOT).as_posix(),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "analysis_families_kept_separate": {"core_rg": 396, "extension_rg": 1200, "validation_candidates": 217},
        "counts": {"rg_core": 396, "rg_extension": 1200, "rg_external_estimates": 41, "heritability": 158, "sensitivities": 62, "sleep_constructs": 12, "outcome_phenotypes": 158},
        "files": manifest_files,
        "input_files": source_records,
        "excluded_fields": ["local paths", "source-body hashes", "raw GWAS rows", "credentials or signed URLs", "intermediate logs"],
        "rights_note": "This is a review candidate. Source owners/investigators must confirm rights and metadata before public deposition.",
    }
    (out / "release_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"CANDIDATE_RELEASE_OK files={len(manifest_files)} rg={sum(rg_counts.values())} h2={len(h2_outputs)} sensitivity={len(sensitivities)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
