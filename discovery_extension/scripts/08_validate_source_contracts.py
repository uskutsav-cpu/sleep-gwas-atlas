#!/usr/bin/env python3
"""Validate extension source contracts without claiming unrun analyses."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path("discovery_extension")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    panel_check = subprocess.run(
        ["python3", str(ROOT / "scripts/05_validate_extension_panel.py")],
        check=False, text=True, capture_output=True,
    )
    if panel_check.returncode:
        fail(f"panel validation failed:\n{panel_check.stdout}{panel_check.stderr}")

    panel_path = ROOT / "config/candidate_traits.tsv"
    panel = read_tsv(panel_path)
    ids = [row["extension_trait_id"] for row in panel]
    sources = read_tsv(ROOT / "config/extension_gwas_sources.tsv")
    if [row["extension_trait_id"] for row in sources] != ids:
        fail("source contract is not a one-to-one ordered copy of the locked panel")
    by_id = {row["extension_trait_id"]: row for row in panel}
    for source in sources:
        locked = by_id[source["extension_trait_id"]]
        expected = {
            "source_url": locked["source_url"],
            "source_checksum": locked["checksum"],
            "source_file_size_bytes": locked["source_file_size_bytes"],
            "tabix_url": locked["source_tabix_url"],
            "tabix_checksum": locked["source_tabix_checksum"],
            "tabix_size_bytes": locked["source_tabix_size_bytes"],
            "build": "GRCh37",
            "ancestry": "EUR",
        }
        for field, value in expected.items():
            if source[field] != value:
                fail(f"source mismatch for {source['extension_trait_id']} field {field}")

    schemas = {row["schema_id"]: row for row in read_tsv(ROOT / "config/extension_gwas_schemas.tsv")}
    expected_schemas = {"PANUKBB_EUR_CONTINUOUS_V1", "PANUKBB_EUR_BINARY_V1"}
    if set(schemas) != expected_schemas:
        fail("exact continuous/binary Pan-UKB schema pair is required")
    for schema in schemas.values():
        for field, value in {
            "status": "HEADER_RANGE_VERIFIED",
            "chromosome": "chr", "position": "pos", "effect_allele": "alt",
            "other_allele": "ref", "beta": "beta_EUR", "standard_error": "se_EUR",
            "p_value_source": "neglog10_pval_EUR", "low_confidence": "low_confidence_EUR",
            "imputation_info": "JOIN_VARIANT_MANIFEST.info",
        }.items():
            if schema[field] != value:
                fail(f"invalid schema {schema['schema_id']} field {field}")

    mappings = read_tsv(ROOT / "config/extension_variant_mapping_plan.tsv")
    if len(mappings) != 1:
        fail("exactly one shared variant mapping contract is required")
    mapping = mappings[0]
    if mapping["join_key"] != "exact(chrom,pos,ref,alt)" or mapping["build"] != "GRCh37":
        fail("variant identity mapping must use the exact GRCh37 coordinate/allele key")
    hm3 = Path(mapping["ldsc_reference"])
    if not hm3.is_file() or sha256(hm3) != mapping["ldsc_reference_sha256"]:
        fail("pinned HapMap3 mapping is missing or changed")
    if int(mapping["source_size_bytes"]) != 2701503051 or int(mapping["tabix_size_bytes"]) != 2261519:
        fail("variant-manifest object sizes differ from verified HTTP metadata")

    policy = json.loads((ROOT / "config/extension_harmonization_policy.json").read_text())
    if policy["panel_sha256"] != sha256(panel_path):
        fail("harmonization policy is bound to a different panel hash")
    if policy["primary_h2_gate"] != {
        "LDSC_intercept_maximum_inclusive": 1.2,
        "failure_disposition": "exclude trait from primary extension rg/FDR family; retain in sensitivity ledger without replacement",
        "h2_z_minimum_inclusive": 4.0,
        "rerun_required_after_harmonization": True,
    }:
        fail("primary h2 gate differs from the locked extension contract")
    filters = "\n".join(policy["filters_in_order"])
    for required in ("INFO > 0.9", "MAF > 0.01", "MHC chr6:25000000-34000000", "strand-ambiguous"):
        if required not in filters:
            fail(f"harmonization policy is missing: {required}")
    if not policy["primary_rg_family"]["never_merge_with_core_family"]:
        fail("extension FDR family is not isolated from the immutable core")

    pilots = json.loads((ROOT / "provenance/panukbb/schema_pilots.json").read_text())
    pilot_by_id = {row["schema_id"]: row for row in pilots["probes"]}
    required_headers = {
        "PANUKBB_EUR_CONTINUOUS_V1": {"chr", "pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"},
        "PANUKBB_EUR_BINARY_V1": {"chr", "pos", "ref", "alt", "af_cases_EUR", "af_controls_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"},
        "PANUKBB_FULL_VARIANT_QC_V1": {"chrom", "pos", "ref", "alt", "rsid", "varid", "high_quality", "info"},
    }
    if set(pilot_by_id) != set(required_headers):
        fail("schema pilot registry is incomplete")
    for schema_id, required in required_headers.items():
        actual = set(pilot_by_id[schema_id]["header"].split("\t"))
        if not required.issubset(actual):
            fail(f"schema pilot {schema_id} lacks required fields {sorted(required - actual)}")

    h2 = read_tsv(ROOT / "results/extension_source_h2_qc.tsv")
    if [row["extension_trait_id"] for row in h2] != ids:
        fail("source h2 ledger is not one-to-one with the lock")
    if any(row["source_precheck"] != "PASS" for row in h2):
        fail("a locked trait did not pass the documented prospective source h2 screen")
    if any(row["primary_rg_eligibility"] != "PENDING_EXTENSION_RERUN_H2" for row in h2):
        fail("h2 ledger incorrectly claims a completed extension rerun")

    preflight = json.loads((ROOT / "provenance/acquisition_preflight.json").read_text())
    expected_total = (
        sum(int(row["source_file_size_bytes"]) for row in panel)
        + sum(int(row["source_tabix_size_bytes"]) for row in panel)
        + 2701503051 + 2261519
    )
    if preflight["compressed_source_bytes"] != expected_total:
        fail("storage preflight total does not match the locked source contract")
    if preflight["download_started"]:
        fail("preflight claims a download started even though no acquisition ledger exists")
    acquisition = read_tsv(ROOT / "results/acquisition_plan.tsv")
    if len(acquisition) != len(panel) + 2:
        fail("acquisition plan must contain 100 phenotype files plus two shared references")

    print(
        f"SOURCE_CONTRACTS_VALID traits={len(panel)} schemas={len(schemas)} "
        f"source_gib={preflight['compressed_source_gib']} status={preflight['status']}"
    )


if __name__ == "__main__":
    main()
