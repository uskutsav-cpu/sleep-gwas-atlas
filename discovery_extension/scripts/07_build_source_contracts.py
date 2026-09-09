#!/usr/bin/env python3
"""Build fail-closed source, schema, mapping, and h2-QC contracts for the lock."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


CONTINUOUS_HEADER = (
    "chr\tpos\tref\talt\taf_meta_hq\tbeta_meta_hq\tse_meta_hq\t"
    "neglog10_pval_meta_hq\tneglog10_pval_heterogeneity_hq\taf_meta\t"
    "beta_meta\tse_meta\tneglog10_pval_meta\tneglog10_pval_heterogeneity\t"
    "af_AFR\taf_AMR\taf_CSA\taf_EAS\taf_EUR\taf_MID\tbeta_AFR\tbeta_AMR\t"
    "beta_CSA\tbeta_EAS\tbeta_EUR\tbeta_MID\tse_AFR\tse_AMR\tse_CSA\t"
    "se_EAS\tse_EUR\tse_MID\tneglog10_pval_AFR\tneglog10_pval_AMR\t"
    "neglog10_pval_CSA\tneglog10_pval_EAS\tneglog10_pval_EUR\t"
    "neglog10_pval_MID\tlow_confidence_AFR\tlow_confidence_AMR\t"
    "low_confidence_CSA\tlow_confidence_EAS\tlow_confidence_EUR\t"
    "low_confidence_MID"
)
BINARY_HEADER = (
    "chr\tpos\tref\talt\taf_cases_meta_hq\taf_controls_meta_hq\tbeta_meta_hq\t"
    "se_meta_hq\tneglog10_pval_meta_hq\tneglog10_pval_heterogeneity_hq\t"
    "af_cases_meta\taf_controls_meta\tbeta_meta\tse_meta\tneglog10_pval_meta\t"
    "neglog10_pval_heterogeneity\taf_cases_AFR\taf_cases_CSA\taf_cases_EUR\t"
    "af_cases_MID\taf_controls_AFR\taf_controls_CSA\taf_controls_EUR\t"
    "af_controls_MID\tbeta_AFR\tbeta_CSA\tbeta_EUR\tbeta_MID\tse_AFR\tse_CSA\t"
    "se_EUR\tse_MID\tneglog10_pval_AFR\tneglog10_pval_CSA\t"
    "neglog10_pval_EUR\tneglog10_pval_MID\tlow_confidence_AFR\t"
    "low_confidence_CSA\tlow_confidence_EUR\tlow_confidence_MID"
)
VARIANT_HEADER = (
    "chrom\tpos\tref\talt\trsid\tvarid\tpass_gnomad_genomes\t"
    "n_passing_populations\thigh_quality\tnearest_genes\tinfo\tac_AFR\taf_AFR\t"
    "an_AFR\tgnomad_genomes_ac_AFR\tgnomad_genomes_af_AFR\t"
    "gnomad_genomes_an_AFR\tac_AMR\taf_AMR\tan_AMR\tgnomad_genomes_ac_AMR\t"
    "gnomad_genomes_af_AMR\tgnomad_genomes_an_AMR\tac_CSA\taf_CSA\tan_CSA\t"
    "ac_EAS\taf_EAS\tan_EAS\tgnomad_genomes_ac_EAS\tgnomad_genomes_af_EAS\t"
    "gnomad_genomes_an_EAS\tac_EUR\taf_EUR\tan_EUR\tgnomad_genomes_ac_EUR\t"
    "gnomad_genomes_af_EUR\tgnomad_genomes_an_EUR\tac_MID\taf_MID\tan_MID"
)
VARIANT_MANIFEST_URL = (
    "https://pan-ukb-us-east-1.s3.amazonaws.com/"
    "sumstats_release/full_variant_qc_metrics.txt.bgz"
)
VARIANT_TABIX_URL = VARIANT_MANIFEST_URL + ".tbi"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel", type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--lock", type=Path,
        default=Path("discovery_extension/config/extension_panel.lock.json"),
    )
    parser.add_argument(
        "--preflight", type=Path,
        default=Path("discovery_extension/provenance/acquisition_preflight.json"),
    )
    parser.add_argument(
        "--sources-out", type=Path,
        default=Path("discovery_extension/config/extension_gwas_sources.tsv"),
    )
    parser.add_argument(
        "--schemas-out", type=Path,
        default=Path("discovery_extension/config/extension_gwas_schemas.tsv"),
    )
    parser.add_argument(
        "--mapping-out", type=Path,
        default=Path("discovery_extension/config/extension_variant_mapping_plan.tsv"),
    )
    parser.add_argument(
        "--policy-out", type=Path,
        default=Path("discovery_extension/config/extension_harmonization_policy.json"),
    )
    parser.add_argument(
        "--pilots-out", type=Path,
        default=Path("discovery_extension/provenance/panukbb/schema_pilots.json"),
    )
    parser.add_argument(
        "--h2-out", type=Path,
        default=Path("discovery_extension/results/extension_source_h2_qc.tsv"),
    )
    args = parser.parse_args()

    validator = subprocess.run(
        ["python3", "discovery_extension/scripts/05_validate_extension_panel.py"],
        check=False, text=True, capture_output=True,
    )
    if validator.returncode:
        raise SystemExit(f"ERROR: extension panel validation failed:\n{validator.stdout}{validator.stderr}")

    panel = read_tsv(args.panel)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    acquisition_status = preflight["status"]

    source_fields = [
        "extension_trait_id", "phenotype_name", "phenotype_domain", "source_release",
        "study_accession", "PMID", "DOI", "ancestry", "build", "schema_id",
        "source_url", "source_checksum", "source_filename", "source_file_size_bytes",
        "tabix_url", "tabix_checksum", "tabix_filename", "tabix_size_bytes",
        "acquisition_status", "analysis_readiness", "retention_policy",
    ]
    source_rows: list[dict[str, object]] = []
    for row in panel:
        source_rows.append(
            {
                "extension_trait_id": row["extension_trait_id"],
                "phenotype_name": row["phenotype_name"],
                "phenotype_domain": row["phenotype_domain"],
                "source_release": row["source"],
                "study_accession": row["study_accession"],
                "PMID": row["PMID"],
                "DOI": row["DOI"],
                "ancestry": row["ancestry"],
                "build": row["build"],
                "schema_id": f"PANUKBB_EUR_{row['binary_or_continuous'].upper()}_V1",
                "source_url": row["source_url"],
                "source_checksum": row["checksum"],
                "source_filename": row["source_filename"],
                "source_file_size_bytes": row["source_file_size_bytes"],
                "tabix_url": row["source_tabix_url"],
                "tabix_checksum": row["source_tabix_checksum"],
                "tabix_filename": row["source_tabix_filename"],
                "tabix_size_bytes": row["source_tabix_size_bytes"],
                "acquisition_status": acquisition_status,
                "analysis_readiness": "BLOCKED_RAW_SOURCE_NOT_LOCAL" if acquisition_status != "READY" else "ACQUISITION_AUTHORIZATION_PENDING",
                "retention_policy": "retain_exact_full_resolution_bgzip_and_tabix_with_checksum",
            }
        )
    write_tsv(args.sources_out, source_fields, source_rows)

    schema_fields = [
        "schema_id", "status", "container", "delimiter", "chromosome", "position",
        "effect_allele", "other_allele", "beta", "standard_error", "p_value_source",
        "p_value_transform", "effect_allele_frequency", "imputation_info",
        "low_confidence", "sample_size", "effect_scale", "notes",
    ]
    schemas = [
        {
            "schema_id": "PANUKBB_EUR_CONTINUOUS_V1", "status": "HEADER_RANGE_VERIFIED",
            "container": "bgzip_tsv", "delimiter": "TAB", "chromosome": "chr",
            "position": "pos", "effect_allele": "alt", "other_allele": "ref",
            "beta": "beta_EUR", "standard_error": "se_EUR",
            "p_value_source": "neglog10_pval_EUR",
            "p_value_transform": "P=max(10^(-neglog10_pval_EUR),1e-300)",
            "effect_allele_frequency": "af_EUR", "imputation_info": "JOIN_VARIANT_MANIFEST.info",
            "low_confidence": "low_confidence_EUR", "sample_size": "CONSTANT_N_FROM_LOCKED_MANIFEST",
            "effect_scale": "beta_per_alt_effect_allele",
            "notes": "Full header verified from official 1 MiB range probe; coordinates are GRCh37.",
        },
        {
            "schema_id": "PANUKBB_EUR_BINARY_V1", "status": "HEADER_RANGE_VERIFIED",
            "container": "bgzip_tsv", "delimiter": "TAB", "chromosome": "chr",
            "position": "pos", "effect_allele": "alt", "other_allele": "ref",
            "beta": "beta_EUR", "standard_error": "se_EUR",
            "p_value_source": "neglog10_pval_EUR",
            "p_value_transform": "P=max(10^(-neglog10_pval_EUR),1e-300)",
            "effect_allele_frequency": "weighted(af_cases_EUR,af_controls_EUR;cases,controls)",
            "imputation_info": "JOIN_VARIANT_MANIFEST.info", "low_confidence": "low_confidence_EUR",
            "sample_size": "CONSTANT_N_EFF=4/(1/cases+1/controls)",
            "effect_scale": "log_odds_or_linear_mixed_model_beta_per_alt_effect_allele",
            "notes": "Full header verified from official 1 MiB range probe; phenotype-specific model interpretation remains in Pan-UKB metadata.",
        },
    ]
    write_tsv(args.schemas_out, schema_fields, schemas)

    mapping_fields = [
        "mapping_id", "status", "source_url", "source_checksum", "source_size_bytes",
        "tabix_url", "tabix_checksum", "tabix_size_bytes", "build", "join_key",
        "rsid_field", "info_field", "ambiguity_policy", "ldsc_reference",
        "ldsc_reference_sha256", "notes",
    ]
    hm3 = Path("ref/hm3_grch37_variant_map.tsv.gz")
    mapping_rows = [{
        "mapping_id": "PANUKBB_FULL_VARIANT_QC_V1",
        "status": "HEADER_AND_REMOTE_OBJECT_METADATA_VERIFIED__FULL_FILE_NOT_LOCAL",
        "source_url": VARIANT_MANIFEST_URL,
        "source_checksum": "etag:e70ebc8289f762dd8d5086f54e766654",
        "source_size_bytes": 2701503051,
        "tabix_url": VARIANT_TABIX_URL,
        "tabix_checksum": "etag:555b93fcb299e59c6e75a5a41ec3de40",
        "tabix_size_bytes": 2261519,
        "build": "GRCh37",
        "join_key": "exact(chrom,pos,ref,alt)",
        "rsid_field": "rsid",
        "info_field": "info",
        "ambiguity_policy": "fail_or_drop_on_missing_non_rsid_duplicate_or_conflicting_exact_key;never_guess",
        "ldsc_reference": str(hm3),
        "ldsc_reference_sha256": sha256(hm3),
        "notes": "Dense path joins the full Pan-UKB variant manifest; LDSC path then restricts to the pinned EUR HapMap3 identity map.",
    }]
    write_tsv(args.mapping_out, mapping_fields, mapping_rows)

    policy = {
        "schema_version": "1.0.0",
        "panel_name": lock["panel_name"],
        "panel_sha256": sha256(args.panel),
        "selection_timing": lock["selection_timing"],
        "input_build": "GRCh37",
        "output_build": "GRCh37",
        "canonical_columns": ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N", "INFO"],
        "allele_semantics": {"A1": "source alt and beta effect allele", "A2": "source ref"},
        "binary_frequency": "(af_cases_EUR*cases + af_controls_EUR*controls)/(cases+controls)",
        "binary_effective_n": "4/(1/cases+1/controls)",
        "p_value": {"source": "neglog10_pval_EUR", "transform": "10**(-x)", "numeric_floor": 1e-300},
        "filters_in_order": [
            "source low_confidence_EUR must equal false",
            "exact variant-manifest join on chrom,pos,ref,alt",
            "valid rsID required; no guessed or ambiguous identity",
            "autosomal chromosomes 1-22 only",
            "single-base A/C/G/T alleles only; A1 must differ from A2",
            "drop all strand-ambiguous A/T and C/G SNPs",
            "finite beta; finite positive SE; 0<P<=1; 0<=FRQ<=1",
            "variant-manifest imputation INFO >= 0.9 and <= 1",
            "common variant MAF >= 0.01",
            "exclude extended MHC chr6:25000000-34000000 inclusive",
            "one row per rsID; duplicate/conflicting rows fail closed",
        ],
        "ldsc_restriction": {
            "reference": str(hm3),
            "reference_sha256": sha256(hm3),
            "policy": "retain only exact GRCh37 HapMap3 coordinate-and-allele matches",
        },
        "primary_h2_gate": {
            "rerun_required_after_harmonization": True,
            "h2_z_minimum_inclusive": 4.0,
            "LDSC_intercept_maximum_inclusive": 1.2,
            "failure_disposition": "exclude trait from primary extension rg/FDR family; retain in sensitivity ledger without replacement",
        },
        "primary_rg_family": {
            "sleep_trait_count": 12,
            "locked_extension_trait_count": len(panel),
            "planned_raw_tests_before_rerun_h2_exclusion": 12 * len(panel),
            "multiple_testing": "Benjamini-Hochberg across all primary extension pairs only",
            "fdr_threshold": 0.05,
            "never_merge_with_core_family": True,
        },
        "analysis_status": "BLOCKED_RAW_SOURCE_NOT_LOCAL" if acquisition_status != "READY" else "SOURCE_ACQUISITION_PENDING",
    }
    write_json(args.policy_out, policy)

    pilots = {
        "schema_version": "1.0.0",
        "verified_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "method": "HTTP Range bytes=0-1048575 from official Pan-UKB S3 objects; full releases were not downloaded",
        "probes": [
            {
                "schema_id": "PANUKBB_EUR_CONTINUOUS_V1",
                "url": "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30520-both_sexes-irnt.tsv.bgz",
                "full_object_size_bytes": 2287347717,
                "full_object_md5_from_official_manifest": "c6316476ff81d15cb11aaf82ea3ab4c9",
                "range": "bytes=0-1048575",
                "range_sha256": "daf6dcfc7e94688b4ab5453ea61534b120a769525127662de208e176d3d3afe2",
                "header": CONTINUOUS_HEADER,
            },
            {
                "schema_id": "PANUKBB_EUR_BINARY_V1",
                "url": "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/phecode-594-both_sexes.tsv.bgz",
                "full_object_size_bytes": 2089852693,
                "full_object_md5_from_official_manifest": "2fedaf06312a142a31e404616e338a91",
                "range": "bytes=0-1048575",
                "range_sha256": "20ffa66760469cc1fc90a9f781378cf1a4c68bf380a0cb830b4d178b3450177c",
                "header": BINARY_HEADER,
            },
            {
                "schema_id": "PANUKBB_FULL_VARIANT_QC_V1",
                "url": VARIANT_MANIFEST_URL,
                "full_object_size_bytes": 2701503051,
                "full_object_etag": "e70ebc8289f762dd8d5086f54e766654",
                "last_modified": "2020-08-28T23:21:14Z",
                "range": "bytes=0-1048575",
                "range_sha256": "785bed25b735e9f2324eeeb38c17fbd0baa449fd5dcfd113de490c0ac4f954a0",
                "header": VARIANT_HEADER,
            },
        ],
        "limitations": "Range probes verify headers and starting rows, not full-file integrity. Full checksums and complete traversal remain mandatory after acquisition.",
    }
    write_json(args.pilots_out, pilots)

    h2_fields = [
        "extension_trait_id", "phenotype_name", "phenotype_domain", "source_h2",
        "source_h2_se", "source_h2_z", "source_LDSC_intercept", "source_precheck",
        "extension_rerun_h2_status", "primary_rg_eligibility", "failure_disposition",
    ]
    h2_rows: list[dict[str, object]] = []
    for row in panel:
        source_pass = float(row["source_h2_z"]) >= 4 and float(row["source_LDSC_intercept"]) <= 1.2
        h2_rows.append({
            "extension_trait_id": row["extension_trait_id"],
            "phenotype_name": row["phenotype_name"],
            "phenotype_domain": row["phenotype_domain"],
            "source_h2": row["source_h2"],
            "source_h2_se": row["source_h2_se"],
            "source_h2_z": row["source_h2_z"],
            "source_LDSC_intercept": row["source_LDSC_intercept"],
            "source_precheck": "PASS" if source_pass else "FAIL",
            "extension_rerun_h2_status": "BLOCKED_RAW_SOURCE_NOT_LOCAL" if acquisition_status != "READY" else "PENDING_HARMONIZATION",
            "primary_rg_eligibility": "PENDING_EXTENSION_RERUN_H2",
            "failure_disposition": "sensitivity_only_without_replacement_if_rerun_h2_z_lt_4_or_intercept_gt_1.2",
        })
    write_tsv(args.h2_out, h2_fields, h2_rows)

    print(
        f"SOURCE_CONTRACTS_OK traits={len(panel)} schemas={len(schemas)} "
        f"status={policy['analysis_status']} panel_sha256={policy['panel_sha256']}"
    )


if __name__ == "__main__":
    main()
