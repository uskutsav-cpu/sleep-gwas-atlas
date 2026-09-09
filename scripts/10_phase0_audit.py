#!/usr/bin/env python3
"""Derive evidence-backed readiness gates for every locked panel trait.

The manifest declares source selection only. This audit derives later stages
from source-registry evidence, metadata, local files, and (optionally) the h2
QC table. A missing artifact is a visible blocker, never an implicit status.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "panel_version", "trait_id", "label", "domain", "type",
    "phenotype_definition", "source_note", "dataset_version", "source_id",
    "raw_file", "ncase", "ncontrol", "n_total", "pop_prev", "build",
    "source_status", "pmid", "doi", "ancestry", "pop_prev_citation",
}
SOURCE_COLUMNS = {
    "source_id", "trait_ids", "source_page_url", "download_url", "access",
    "archive_sha256", "raw_files", "pmid", "ancestry_reported",
    "build_status", "acquisition_status",
}
SCHEMA_COLUMNS = {
    "source_id", "trait_id", "schema_status", "file_format", "variant_id", "chromosome",
    "position", "effect_allele", "other_allele", "effect",
    "effect_convention", "standard_error", "p_value", "eaf", "info",
    "sample_size", "notes",
}
VARIANT_MAPPING_COLUMNS = {
    "source_id", "trait_id", "strategy", "map_path", "map_bytes",
    "map_sha256", "input_identity", "output_build", "notes",
}
VARIANT_MAPPING_STRATEGIES = {"BY_COORD_ALLELES", "BY_RSID_ALLELES"}
VARIANT_MAP_SCHEMA = "atlas.hm3-grch37-variant-map.v1"
LIFTOVER_COLUMNS = {
    "source_id", "trait_id", "strategy", "chain_path", "chain_bytes",
    "chain_md5", "chain_sha256", "source_build", "output_build", "notes",
}
HM3_PREFILTER_COLUMNS = {
    "source_id", "trait_id", "strategy", "snp_column", "allowlist_path",
    "rationale",
}
HM3_PREFILTER_STRATEGY = "HAPMAP3_RSID_ALLOWLIST"
NONJOURNAL_CITATION_COLUMNS = {
    "source_id", "trait_id", "citation_type", "title", "publisher",
    "release_date", "source_url", "notes",
}
NONJOURNAL_CITATION_TYPE = "PUBLIC_DATA_RELEASE"
LIFTOVER_STRATEGY = "UCSC_CHAIN_POINT"
MISSING_TEXT = {"", "na", "nan", "none", "null", "unresolved", "unregistered"}


def populated(value):
    if value is None:
        return False
    normalized = str(value).strip().lower()
    return normalized not in MISSING_TEXT and not normalized.startswith("unresolved")


def numeric(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def manifest_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_public_sources(path, trait_ids):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: source registry not found: {path}")
    sources = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = SOURCE_COLUMNS.difference(sources.columns)
    if missing:
        raise SystemExit(f"ERROR: source registry missing columns: {sorted(missing)}")
    by_trait = {}
    for _, source in sources.iterrows():
        mapped_ids = [item.strip() for item in source["trait_ids"].split(",") if item.strip()]
        for trait_id in mapped_ids:
            if trait_id not in trait_ids:
                continue  # source registry may retain archived candidate sources
            if trait_id in by_trait:
                raise SystemExit(f"ERROR: multiple source rows claim trait_id {trait_id}")
            by_trait[trait_id] = source
    return by_trait


def load_nonjournal_citations(path, selected_sources):
    """Load explicit primary-source citations for releases without a paper.

    This is intentionally narrow: a citation must name the exact selected
    source/trait pair and describe an official public data release.  It never
    supplies a surrogate PMID or relaxes any other source-readiness check.
    """
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: non-journal citation registry not found: {path}")
    citations = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = NONJOURNAL_CITATION_COLUMNS.difference(citations.columns)
    if missing:
        raise SystemExit(
            f"ERROR: non-journal citation registry missing columns: {sorted(missing)}"
        )
    duplicate_traits = citations.duplicated(subset=["trait_id"], keep=False)
    if duplicate_traits.any():
        duplicates = citations.loc[
            duplicate_traits, ["source_id", "trait_id"]
        ].to_dict("records")
        raise SystemExit(f"ERROR: duplicate non-journal citation rows: {duplicates}")
    invalid = []
    for _, citation in citations.iterrows():
        selected_source = selected_sources.get(citation["trait_id"])
        if selected_source != citation["source_id"]:
            invalid.append({
                "source_id": citation["source_id"],
                "trait_id": citation["trait_id"],
                "selected_source_id": selected_source or "UNSELECTED_TRAIT",
            })
        if citation["citation_type"] != NONJOURNAL_CITATION_TYPE:
            invalid.append({
                "trait_id": citation["trait_id"],
                "invalid_citation_type": citation["citation_type"],
            })
        for field in ["title", "publisher", "release_date", "source_url"]:
            if not populated(citation[field]):
                invalid.append({
                    "trait_id": citation["trait_id"],
                    "missing_field": field,
                })
        if populated(citation["source_url"]) and not citation["source_url"].startswith("https://"):
            invalid.append({
                "trait_id": citation["trait_id"],
                "invalid_source_url": citation["source_url"],
            })
    if invalid:
        raise SystemExit(f"ERROR: invalid non-journal citations: {invalid}")
    return {row["trait_id"]: row for _, row in citations.iterrows()}


def load_source_schemas(path, selected_sources):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: GWAS schema registry not found: {path}")
    schemas = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = SCHEMA_COLUMNS.difference(schemas.columns)
    if missing:
        raise SystemExit(f"ERROR: GWAS schema registry missing columns: {sorted(missing)}")
    duplicate_pairs = schemas.duplicated(subset=["source_id", "trait_id"], keep=False)
    if duplicate_pairs.any():
        duplicates = schemas.loc[duplicate_pairs, ["source_id", "trait_id"]].to_dict("records")
        raise SystemExit(f"ERROR: duplicate source/trait schema rows: {duplicates}")
    invalid_pairs = []
    for _, schema in schemas.iterrows():
        selected_source = selected_sources.get(schema["trait_id"])
        if selected_source != schema["source_id"]:
            invalid_pairs.append({
                "source_id": schema["source_id"],
                "trait_id": schema["trait_id"],
                "selected_source_id": selected_source or "UNSELECTED_TRAIT",
            })
    if invalid_pairs:
        raise SystemExit(f"ERROR: source schemas do not match selected trait sources: {invalid_pairs}")
    allowed = {
        "SCHEMA_PENDING",
        "SCHEMA_VERIFIED",
        "SCHEMA_VERIFIED_REQUIRES_VARIANT_MAPPING",
    }
    invalid = schemas.loc[~schemas["schema_status"].isin(allowed), ["source_id", "schema_status"]]
    if len(invalid):
        raise SystemExit(f"ERROR: invalid source schema status: {invalid.to_dict('records')}")
    return {(row["source_id"], row["trait_id"]): row for _, row in schemas.iterrows()}


def load_variant_mappings(path, selected_sources):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: variant-mapping registry not found: {path}")
    mappings = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = VARIANT_MAPPING_COLUMNS.difference(mappings.columns)
    if missing:
        raise SystemExit(f"ERROR: variant-mapping registry missing columns: {sorted(missing)}")
    duplicate_traits = mappings.duplicated(subset=["trait_id"], keep=False)
    if duplicate_traits.any():
        duplicates = mappings.loc[duplicate_traits, ["source_id", "trait_id"]].to_dict("records")
        raise SystemExit(f"ERROR: duplicate variant-mapping trait rows: {duplicates}")
    invalid = []
    for _, mapping in mappings.iterrows():
        selected_source = selected_sources.get(mapping["trait_id"])
        if selected_source != mapping["source_id"]:
            invalid.append({
                "source_id": mapping["source_id"],
                "trait_id": mapping["trait_id"],
                "selected_source_id": selected_source or "UNSELECTED_TRAIT",
            })
        if mapping["strategy"] not in VARIANT_MAPPING_STRATEGIES:
            invalid.append({
                "trait_id": mapping["trait_id"],
                "invalid_strategy": mapping["strategy"],
            })
        if mapping["output_build"] != "GRCh37/hg19":
            invalid.append({
                "trait_id": mapping["trait_id"],
                "invalid_output_build": mapping["output_build"],
            })
        if not mapping["map_bytes"].isdigit() or int(mapping["map_bytes"]) <= 0:
            invalid.append({
                "trait_id": mapping["trait_id"],
                "invalid_map_bytes": mapping["map_bytes"],
            })
        if not re.fullmatch(r"[0-9a-f]{64}", mapping["map_sha256"]):
            invalid.append({
                "trait_id": mapping["trait_id"],
                "invalid_map_sha256": mapping["map_sha256"],
            })
    if invalid:
        raise SystemExit(f"ERROR: invalid variant-mapping plans: {invalid}")
    return {row["trait_id"]: row for _, row in mappings.iterrows()}


def variant_map_issues(mapping):
    path = mapping["map_path"]
    provenance_path = f"{path}.provenance.json"
    issues = []
    if not os.path.isfile(path):
        issues.append("variant_mapping_reference_missing")
        return issues
    if not os.path.isfile(provenance_path):
        issues.append("variant_mapping_provenance_missing")
        return issues
    try:
        with open(provenance_path, encoding="utf-8") as handle:
            provenance = json.load(handle)
    except (OSError, json.JSONDecodeError):
        issues.append("variant_mapping_provenance_invalid")
        return issues
    if provenance.get("schema_version") != VARIANT_MAP_SCHEMA:
        issues.append("variant_mapping_schema_invalid")
    if provenance.get("genome_build") != "GRCh37/hg19":
        issues.append("variant_mapping_build_invalid")
    actual_sha256 = manifest_sha256(path)
    actual_bytes = os.path.getsize(path)
    if actual_sha256 != mapping["map_sha256"]:
        issues.append("variant_mapping_registered_checksum_mismatch")
    if actual_bytes != int(mapping["map_bytes"]):
        issues.append("variant_mapping_registered_byte_count_mismatch")
    if provenance.get("map_sha256") != actual_sha256:
        issues.append("variant_mapping_checksum_mismatch")
    if provenance.get("map_bytes") != actual_bytes:
        issues.append("variant_mapping_byte_count_mismatch")
    if not numeric(provenance.get("mapped_rows")):
        issues.append("variant_mapping_row_count_missing")
    return issues


def load_liftover_plans(path, selected_sources):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: liftover-plan registry not found: {path}")
    plans = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = LIFTOVER_COLUMNS.difference(plans.columns)
    if missing:
        raise SystemExit(f"ERROR: liftover-plan registry missing columns: {sorted(missing)}")
    duplicate_traits = plans.duplicated(subset=["trait_id"], keep=False)
    if duplicate_traits.any():
        duplicates = plans.loc[duplicate_traits, ["source_id", "trait_id"]].to_dict("records")
        raise SystemExit(f"ERROR: duplicate liftover trait rows: {duplicates}")
    invalid = []
    for _, plan in plans.iterrows():
        selected_source = selected_sources.get(plan["trait_id"])
        if selected_source != plan["source_id"]:
            invalid.append({
                "source_id": plan["source_id"],
                "trait_id": plan["trait_id"],
                "selected_source_id": selected_source or "UNSELECTED_TRAIT",
            })
        if plan["strategy"] != LIFTOVER_STRATEGY:
            invalid.append({"trait_id": plan["trait_id"], "invalid_strategy": plan["strategy"]})
        if plan["source_build"] != "GRCh38/hg38" or plan["output_build"] != "GRCh37/hg19":
            invalid.append({
                "trait_id": plan["trait_id"],
                "invalid_build_pair": f"{plan['source_build']}->{plan['output_build']}",
            })
        if not plan["chain_bytes"].isdigit() or int(plan["chain_bytes"]) <= 0:
            invalid.append({"trait_id": plan["trait_id"], "invalid_chain_bytes": plan["chain_bytes"]})
        if not re.fullmatch(r"[0-9a-f]{32}", plan["chain_md5"]):
            invalid.append({"trait_id": plan["trait_id"], "invalid_chain_md5": plan["chain_md5"]})
        if not re.fullmatch(r"[0-9a-f]{64}", plan["chain_sha256"]):
            invalid.append({"trait_id": plan["trait_id"], "invalid_chain_sha256": plan["chain_sha256"]})
    if invalid:
        raise SystemExit(f"ERROR: invalid liftover plans: {invalid}")
    return {row["trait_id"]: row for _, row in plans.iterrows()}


def liftover_chain_issues(plan):
    path = plan["chain_path"]
    if not os.path.isfile(path):
        return ["liftover_chain_missing"]
    issues = []
    if os.path.getsize(path) != int(plan["chain_bytes"]):
        issues.append("liftover_chain_byte_count_mismatch")
    if file_md5(path) != plan["chain_md5"]:
        issues.append("liftover_chain_md5_mismatch")
    if manifest_sha256(path) != plan["chain_sha256"]:
        issues.append("liftover_chain_sha256_mismatch")
    return issues


def load_hm3_prefilter_plans(path, selected_sources):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: HapMap3 prefilter-plan registry not found: {path}")
    plans = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = HM3_PREFILTER_COLUMNS.difference(plans.columns)
    if missing:
        raise SystemExit(
            f"ERROR: HapMap3 prefilter-plan registry missing columns: {sorted(missing)}"
        )
    duplicate_traits = plans.duplicated(subset=["trait_id"], keep=False)
    if duplicate_traits.any():
        duplicates = plans.loc[
            duplicate_traits, ["source_id", "trait_id"]
        ].to_dict("records")
        raise SystemExit(f"ERROR: duplicate HapMap3 prefilter trait rows: {duplicates}")
    invalid = []
    for _, plan in plans.iterrows():
        selected_source = selected_sources.get(plan["trait_id"])
        if selected_source != plan["source_id"]:
            invalid.append({
                "source_id": plan["source_id"],
                "trait_id": plan["trait_id"],
                "selected_source_id": selected_source or "UNSELECTED_TRAIT",
            })
        if plan["strategy"] != HM3_PREFILTER_STRATEGY:
            invalid.append({
                "trait_id": plan["trait_id"],
                "invalid_strategy": plan["strategy"],
            })
        if not populated(plan["snp_column"]) or not populated(plan["allowlist_path"]):
            invalid.append({
                "trait_id": plan["trait_id"],
                "missing_prefilter_contract": True,
            })
    if invalid:
        raise SystemExit(f"ERROR: invalid HapMap3 prefilter plans: {invalid}")
    return {row["trait_id"]: row for _, row in plans.iterrows()}


def retained_prefilter_issues(plan, source, trait, harmonized_dir):
    """Validate a retained prefilter after its registered large raw is removed."""
    trait_id = trait["trait_id"]
    output = os.path.join(harmonized_dir, ".prefilter", f"{trait_id}.hm3.tsv.gz")
    provenance_path = os.path.join(
        harmonized_dir, f"{trait_id}.prefilter.provenance.json"
    )
    issues = []
    if source is None or source.get("archive_member", "") != "DIRECT_GZIP":
        issues.append("prefilter_source_not_registered_direct_gzip")
        return issues
    if source.get("raw_files", "") != trait["raw_file"]:
        issues.append("prefilter_source_raw_file_mismatch")
    if not os.path.isfile(output):
        issues.append("retained_prefilter_missing")
        return issues
    if not os.path.isfile(provenance_path):
        issues.append("retained_prefilter_provenance_missing")
        return issues
    try:
        with open(provenance_path, encoding="utf-8") as handle:
            provenance = json.load(handle)
    except (OSError, json.JSONDecodeError):
        issues.append("retained_prefilter_provenance_invalid")
        return issues
    required = {
        "strategy", "input_bytes", "input_sha256", "snp_column",
        "allowlist_path", "allowlist_sha256", "source_rows", "retained_rows",
        "output_bytes", "output_sha256",
    }
    if required.difference(provenance):
        issues.append("retained_prefilter_provenance_incomplete")
        return issues
    if provenance["strategy"] != plan["strategy"]:
        issues.append("retained_prefilter_strategy_mismatch")
    if provenance["snp_column"] != plan["snp_column"]:
        issues.append("retained_prefilter_snp_column_mismatch")
    try:
        if int(provenance["input_bytes"]) != int(source["archive_bytes"]):
            issues.append("retained_prefilter_source_byte_count_mismatch")
    except (TypeError, ValueError):
        issues.append("retained_prefilter_source_byte_count_invalid")
    if str(provenance["input_sha256"]).lower() != str(source["archive_sha256"]).lower():
        issues.append("retained_prefilter_source_checksum_mismatch")
    allowlist = plan["allowlist_path"]
    if os.path.realpath(str(provenance["allowlist_path"])) != os.path.realpath(allowlist):
        issues.append("retained_prefilter_allowlist_path_mismatch")
    if not os.path.isfile(allowlist):
        issues.append("retained_prefilter_allowlist_missing")
    elif manifest_sha256(allowlist) != str(provenance["allowlist_sha256"]).lower():
        issues.append("retained_prefilter_allowlist_checksum_mismatch")
    actual_bytes = os.path.getsize(output)
    actual_sha256 = manifest_sha256(output)
    try:
        if actual_bytes != int(provenance["output_bytes"]):
            issues.append("retained_prefilter_output_byte_count_mismatch")
        source_rows = int(provenance["source_rows"])
        retained_rows = int(provenance["retained_rows"])
        if source_rows <= 0 or retained_rows <= 0 or retained_rows > source_rows:
            issues.append("retained_prefilter_row_counts_invalid")
    except (TypeError, ValueError):
        issues.append("retained_prefilter_numeric_provenance_invalid")
    if actual_sha256 != str(provenance["output_sha256"]).lower():
        issues.append("retained_prefilter_output_checksum_mismatch")
    return issues


def load_h2(path):
    if not path:
        return {}
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: h2 table not found: {path}")
    h2 = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if not {"trait", "verdict"}.issubset(h2.columns):
        raise SystemExit("ERROR: h2 table must contain trait and verdict columns")
    if h2["trait"].duplicated().any():
        raise SystemExit("ERROR: h2 table contains duplicate trait rows")
    return {row["trait"]: row for _, row in h2.iterrows()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--lock", default="config/analysis_panel.lock.json")
    parser.add_argument("--sources", default="config/public_gwas_sources.tsv")
    parser.add_argument("--schemas", default="config/gwas_schemas.tsv")
    parser.add_argument("--variant-mappings", default="config/variant_mapping_plans.tsv")
    parser.add_argument("--liftover-plans", default="config/liftover_plans.tsv")
    parser.add_argument(
        "--hm3-prefilter-plans", default="config/hm3_prefilter_plans.tsv"
    )
    parser.add_argument(
        "--nonjournal-citations", default="config/nonjournal_source_citations.tsv"
    )
    parser.add_argument("--h2")
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--harmonized-dir", default="data/harmonized")
    parser.add_argument("--munged-dir", default="data/munged")
    parser.add_argument("--out", default="results/tables/trait_readiness.tsv")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    subprocess.run(
        [sys.executable, "scripts/00_validate_panel.py", "--manifest", args.config,
         "--lock", args.lock],
        check=True,
    )
    config = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
    missing = REQUIRED_COLUMNS.difference(config.columns)
    if missing:
        raise SystemExit(f"ERROR: analysis panel missing columns: {sorted(missing)}")
    sources_by_trait = load_public_sources(args.sources, set(config["trait_id"]))
    selected_sources = dict(zip(config["trait_id"], config["source_id"]))
    nonjournal_by_trait = load_nonjournal_citations(
        args.nonjournal_citations, selected_sources
    )
    required_nonjournal_traits = {
        row["trait_id"]
        for _, row in config.iterrows()
        if row["source_status"] == "SOURCE_VERIFIED" and not populated(row["pmid"])
    }
    if set(nonjournal_by_trait) != required_nonjournal_traits:
        raise SystemExit(
            "ERROR: non-journal citations must exactly cover verified sources without a PMID: "
            f"expected {sorted(required_nonjournal_traits)}, got {sorted(nonjournal_by_trait)}"
        )
    schemas_by_source_trait = load_source_schemas(args.schemas, selected_sources)
    mappings_by_trait = load_variant_mappings(args.variant_mappings, selected_sources)
    liftover_by_trait = load_liftover_plans(args.liftover_plans, selected_sources)
    hm3_prefilter_by_trait = load_hm3_prefilter_plans(
        args.hm3_prefilter_plans, selected_sources
    )
    required_mapping_traits = {
        trait_id
        for (_, trait_id), schema in schemas_by_source_trait.items()
        if schema["schema_status"] == "SCHEMA_VERIFIED_REQUIRES_VARIANT_MAPPING"
    }
    if set(mappings_by_trait) != required_mapping_traits:
        raise SystemExit(
            "ERROR: variant-mapping plans must exactly cover mapping-required schemas: "
            f"expected {sorted(required_mapping_traits)}, got {sorted(mappings_by_trait)}"
        )
    required_liftover_traits = {
        row["trait_id"]
        for _, row in config.iterrows()
        if row["source_status"] == "SOURCE_VERIFIED" and row["build"] in {"hg38", "GRCh38"}
    }
    if set(liftover_by_trait) != required_liftover_traits:
        raise SystemExit(
            "ERROR: liftover plans must exactly cover verified hg38 sources: "
            f"expected {sorted(required_liftover_traits)}, got {sorted(liftover_by_trait)}"
        )
    h2_by_trait = load_h2(args.h2)
    panel_hash = manifest_sha256(args.config)

    rows = []
    for _, trait in config.iterrows():
        trait_id = trait["trait_id"]
        source = sources_by_trait.get(trait_id)
        schema = schemas_by_source_trait.get((trait["source_id"], trait_id))
        mapping = mappings_by_trait.get(trait_id)
        mapping_issues = variant_map_issues(mapping) if mapping is not None else []
        liftover = liftover_by_trait.get(trait_id)
        liftover_issues = liftover_chain_issues(liftover) if liftover is not None else []

        source_issues = []
        if trait["source_status"] != "SOURCE_VERIFIED":
            source_issues.append("source_status_pending")
        if source is None:
            source_issues.append("source_not_registered")
        else:
            if trait["source_id"] != source["source_id"]:
                source_issues.append("manifest_source_id_mismatch")
            if source["access"].upper() != "PUBLIC":
                source_issues.append("source_not_public")
            if source["acquisition_status"].upper() != "SOURCE_VERIFIED_PUBLIC":
                source_issues.append("source_acquisition_not_verified")
            if not populated(source["source_page_url"]) or not populated(source["download_url"]):
                source_issues.append("source_url_missing")
            if not populated(source["archive_sha256"]):
                source_issues.append("source_checksum_missing")
            raw_files = {item.strip() for item in source["raw_files"].split(",")}
            if trait["raw_file"] not in raw_files:
                source_issues.append("source_registry_raw_file_mismatch")
        source_verified = not source_issues

        harmonization_issues = []
        if not source_verified:
            harmonization_issues.append("source_not_verified")
        if not populated(trait["phenotype_definition"]):
            harmonization_issues.append("phenotype_definition_unresolved")
        if not populated(trait["dataset_version"]):
            harmonization_issues.append("dataset_version_unresolved")
        if not populated(trait["pmid"]) and trait_id not in nonjournal_by_trait:
            harmonization_issues.append("primary_publication_unresolved")
        if trait["ancestry"].upper() != "EUR":
            harmonization_issues.append("eur_subset_unresolved")
        rsid_mapping_sets_build = (
            mapping is not None and mapping["strategy"] == "BY_RSID_ALLELES"
        )
        liftover_sets_build = liftover is not None
        if (
            liftover_sets_build
            and (
                source is None
                or source["build_status"].upper()
                != "HEADER_VALIDATED_HG38_REQUIRES_LIFTOVER"
            )
        ):
            harmonization_issues.append("liftover_source_build_not_validated_hg38")
        if (
            trait["build"] not in {"hg19", "GRCh37"}
            and not rsid_mapping_sets_build
            and not liftover_sets_build
        ):
            harmonization_issues.append("not_hg19_requires_explicit_build_decision")
        mapping_validates_grch37_coordinate = (
            mapping is not None
            and mapping["strategy"] == "BY_COORD_ALLELES"
            and trait["build"] in {"hg19", "GRCh37"}
        )
        if (
            source is not None
            and source["build_status"].upper() != "HEADER_VALIDATED_HG19"
            and not rsid_mapping_sets_build
            and not mapping_validates_grch37_coordinate
            and not liftover_sets_build
        ):
            harmonization_issues.append("source_build_not_header_validated_hg19")
        harmonization_issues.extend(liftover_issues)
        if schema is None:
            harmonization_issues.append("source_schema_not_registered")
            schema_status = "UNREGISTERED"
        else:
            schema_status = schema["schema_status"]
            if schema_status == "SCHEMA_VERIFIED_REQUIRES_VARIANT_MAPPING":
                harmonization_issues.extend(mapping_issues)
            elif schema_status != "SCHEMA_VERIFIED":
                harmonization_issues.append("source_schema_not_verified")
        if trait["type"] == "binary":
            if not numeric(trait["ncase"]) or not numeric(trait["ncontrol"]):
                harmonization_issues.append("binary_ncase_or_ncontrol_unresolved")
        elif not numeric(trait["n_total"]):
            harmonization_issues.append("continuous_n_total_unresolved")
        raw_path = os.path.join(args.raw_dir, trait["raw_file"])
        if not os.path.isfile(raw_path):
            prefilter = hm3_prefilter_by_trait.get(trait_id)
            if prefilter is None:
                harmonization_issues.append("registered_raw_file_not_materialized")
            else:
                harmonization_issues.extend(
                    retained_prefilter_issues(
                        prefilter, source, trait, args.harmonized_dir
                    )
                )
        harmonization_ready = not harmonization_issues

        ldsc_issues = []
        if not harmonization_ready:
            ldsc_issues.append("harmonization_not_ready")
        qc_path = os.path.join(args.harmonized_dir, f"{trait_id}.qc.txt")
        munged_path = os.path.join(args.munged_dir, f"{trait_id}.sumstats.gz")
        if not os.path.isfile(qc_path):
            ldsc_issues.append("harmonization_qc_ledger_missing")
        if not os.path.isfile(munged_path):
            ldsc_issues.append("hapmap3_sumstats_missing")
        ldsc_ready = not ldsc_issues

        liability_issues = []
        if not ldsc_ready:
            liability_issues.append("ldsc_not_ready")
        if trait["type"] == "binary":
            prevalence = numeric(trait["pop_prev"])
            if prevalence is None or not 0 < prevalence < 1:
                liability_issues.append("population_prevalence_unresolved")
            if not populated(trait["pop_prev_citation"]):
                liability_issues.append("population_prevalence_citation_missing")
        liability_h2_ready = not liability_issues

        phase1_issues = []
        if not liability_h2_ready:
            phase1_issues.append("liability_h2_not_ready")
        h2_row = h2_by_trait.get(trait_id)
        if h2_row is None:
            phase1_issues.append("h2_qc_not_run" if args.h2 else "h2_table_not_supplied")
            h2_verdict = "UNASSESSED"
            h2_reason = "UNASSESSED"
        else:
            h2_verdict = h2_row["verdict"]
            h2_reason = h2_row.get("qc_reason", "") or "UNSPECIFIED"
            if h2_verdict != "PASS":
                phase1_issues.append("h2_qc_not_pass")
        phase1_pass = not phase1_issues

        if phase1_pass:
            stage = "PHASE1_PASS"
        elif liability_h2_ready:
            stage = "LIABILITY_H2_READY"
        elif ldsc_ready:
            stage = "LDSC_READY"
        elif harmonization_ready:
            stage = "HARMONIZATION_READY"
        elif source_verified:
            stage = "SOURCE_VERIFIED"
        else:
            stage = "SOURCE_PENDING"

        rows.append({
            "panel_version": trait["panel_version"],
            "manifest_sha256": panel_hash,
            "trait_id": trait_id,
            "label": trait["label"],
            "domain": trait["domain"],
            "type": trait["type"],
            "declared_source_status": trait["source_status"],
            "schema_status": schema_status,
            "source_verified": source_verified,
            "harmonization_ready": harmonization_ready,
            "ldsc_ready": ldsc_ready,
            "liability_h2_ready": liability_h2_ready,
            "phase1_pass": phase1_pass,
            "readiness_stage": stage,
            "h2_verdict": h2_verdict,
            "h2_qc_reason": h2_reason,
            "source_issues": ";".join(source_issues) if source_issues else "none",
            "harmonization_issues": ";".join(harmonization_issues) if harmonization_issues else "none",
            "ldsc_issues": ";".join(ldsc_issues) if ldsc_issues else "none",
            "liability_h2_issues": ";".join(liability_issues) if liability_issues else "none",
            "phase1_issues": ";".join(phase1_issues) if phase1_issues else "none",
        })

    audit = pd.DataFrame(rows)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    audit.to_csv(args.out, sep="\t", index=False)
    declaration_mismatches = audit[
        (audit["declared_source_status"] == "SOURCE_VERIFIED") & ~audit["source_verified"]
    ]
    print(f"Wrote readiness ledger: {args.out}")
    for stage in ["source_verified", "harmonization_ready", "ldsc_ready", "liability_h2_ready", "phase1_pass"]:
        print(f"{stage}: {int(audit[stage].sum())} / {len(audit)}")
    print(f"SOURCE_VERIFIED declarations without evidence: {len(declaration_mismatches)}")
    if args.strict and len(declaration_mismatches):
        bad = ", ".join(declaration_mismatches["trait_id"].tolist())
        raise SystemExit(f"ERROR: source-verification declaration mismatch: {bad}")


if __name__ == "__main__":
    main()
