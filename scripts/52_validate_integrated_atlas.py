#!/usr/bin/env python3
"""Fail closed unless the ten canonical atlas tables are complete and linked."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


SHA256 = re.compile(r"^[0-9a-f]{64}$")
MISSING = {"", "NA"}
MOLECULAR_COVERAGE_FIELDS = [
    "coverage_id", "locus_id", "pair_id", "trait_id", "modality",
    "source_family_id", "required_unit_count", "completed_unit_count",
    "analyzed_unit_count", "no_evidence_unit_count", "access_blocked_unit_count",
    "not_applicable_unit_count", "coverage_status", "evidence_path",
    "evidence_sha256", "provenance_id",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            fail(f"invalid or duplicate header fields in {path}")
        rows = list(reader)
    return fields, rows


def populated(value: str) -> bool:
    return value not in MISSING


def integer(value: str, field: str, identity: str, minimum: int = 0) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        fail(f"{field} is not an integer for {identity}")
        raise AssertionError from exc
    if parsed < minimum:
        fail(f"{field} is below {minimum} for {identity}")
    return parsed


def number(value: str, field: str, identity: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"{field} is not numeric for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed):
        fail(f"{field} is not finite for {identity}")
    return parsed


def probability(value: str, field: str, identity: str) -> None:
    if value in MISSING:
        return
    parsed = number(value, field, identity)
    if not 0 <= parsed <= 1:
        fail(f"{field} is outside [0,1] for {identity}")


def unique_keys(name: str, rows: list[dict[str, str]], key_fields: list[str]) -> set[tuple[str, ...]]:
    observed: set[tuple[str, ...]] = set()
    for index, row in enumerate(rows, start=2):
        key = tuple(row[field] for field in key_fields)
        if any(not populated(value) for value in key):
            fail(f"missing primary-key value in {name} line {index}")
        if key in observed:
            fail(f"duplicate primary key {key!r} in {name}")
        observed.add(key)
    return observed


def required(value: str, field: str, identity: str) -> str:
    if not populated(value):
        fail(f"missing {field} for {identity}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--atlas-dir", default="results/atlas")
    parser.add_argument("--schema", default="config/atlas_table_schema.json")
    parser.add_argument("--policy", default="config/downstream_analysis_policy.json")
    parser.add_argument(
        "--molecular-coverage", default="results/tables/molecular_locus_coverage.tsv",
    )
    parser.add_argument("--interpretation-policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--interpretation-manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--interpretation-manifest-lock", default="results/tables/interpretation_task_manifest.lock.json")
    parser.add_argument("--interpretation-coverage", default="results/tables/interpretation_coverage.tsv")
    parser.add_argument("--interpretation-provenance", default="results/atlas/interpretation.provenance.json")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    atlas_dir = root / args.atlas_dir
    schema = json.loads((root / args.schema).read_text(encoding="utf-8"))
    policy = json.loads((root / args.policy).read_text(encoding="utf-8"))
    interpretation_policy_path = root / args.interpretation_policy
    interpretation_policy = json.loads(interpretation_policy_path.read_text(encoding="utf-8"))
    definitions = schema.get("tables", {})
    expected_names = policy["integrated_atlas"]["tables"]
    if list(definitions) != expected_names or len(expected_names) != 10:
        fail("atlas schema and downstream policy do not define the same ordered ten-table family")
    if definitions["edges.tsv"]["fields"] != policy["integrated_atlas"]["edge_fields"]:
        fail("edge fields differ between the atlas schema and downstream policy")

    tables: dict[str, list[dict[str, str]]] = {}
    forbidden = [marker.upper() for marker in schema["forbidden_markers"]]
    probability_fields = set(schema["numeric_probability_fields"])
    for name, definition in definitions.items():
        fields, rows = read_tsv(atlas_dir / name)
        if fields != definition["fields"]:
            fail(f"{name} header differs from the frozen schema")
        expected_rows = definition.get("expected_rows")
        if expected_rows is not None and len(rows) != expected_rows:
            fail(f"{name} has {len(rows)} rows; expected exactly {expected_rows}")
        if len(rows) < definition.get("minimum_rows", 0):
            fail(f"{name} has no real evidence rows")
        unique_keys(name, rows, definition["primary_key"])
        for line, row in enumerate(rows, start=2):
            identity = "/".join(row[field] for field in definition["primary_key"])
            payload = "\t".join(row.values()).upper()
            if any(marker in payload for marker in forbidden):
                fail(f"forbidden synthetic/placeholder marker in {name} line {line}")
            required(row["provenance_id"], "provenance_id", identity) if "provenance_id" in row else None
            for field in probability_fields.intersection(row):
                probability(row[field], field, identity)
        tables[name] = rows

    panel_fields, panel = read_tsv(root / "config/analysis_panel.tsv")
    del panel_fields
    trait_ids = [row["trait_id"] for row in panel]
    sleep_ids = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
    non_sleep_ids = [row["trait_id"] for row in panel if row["domain"] != "sleep"]
    if len(trait_ids) != 45 or len(sleep_ids) != 12 or len(non_sleep_ids) != 33:
        fail("analysis panel no longer has the locked 45/12/33 shape")
    if [row["trait_id"] for row in tables["traits.tsv"]] != trait_ids:
        fail("traits.tsv is not in exact locked panel order")
    for row in tables["traits.tsv"]:
        identity = row["trait_id"]
        if row["source_status"] != "SOURCE_VERIFIED":
            fail(f"unverified source in traits.tsv for {identity}")
        for field in ("harmonized_sha256", "munged_sha256"):
            if not SHA256.fullmatch(row[field]):
                fail(f"invalid {field} for {identity}")

    expected_pairs = [(sleep, other) for sleep in sleep_ids for other in non_sleep_ids]
    pair_rows = tables["trait_pairs.tsv"]
    if [(row["sleep_trait"], row["non_sleep_trait"]) for row in pair_rows] != expected_pairs:
        fail("trait_pairs.tsv is not the exact locked 396-pair order")
    pair_ids = set()
    for row in pair_rows:
        expected = f"{row['sleep_trait']}__{row['non_sleep_trait']}"
        if row["pair_id"] != expected:
            fail(f"noncanonical pair_id {row['pair_id']!r}")
        pair_ids.add(expected)

    traits = set(trait_ids)
    loci = tables["loci.tsv"]
    locus_ids = {row["locus_id"] for row in loci}
    for row in loci:
        identity = row["locus_id"]
        if row["sleep_trait"] not in set(sleep_ids) or row["non_sleep_trait"] not in set(non_sleep_ids):
            fail(f"out-of-panel trait link for locus {identity}")
        if f"{row['sleep_trait']}__{row['non_sleep_trait']}" not in pair_ids:
            fail(f"locus {identity} is not linked to a locked pair")
        chromosome = integer(row["chromosome"], "chromosome", identity, 1)
        start = integer(row["start_bp"], "start_bp", identity, 1)
        end = integer(row["end_bp"], "end_bp", identity, 1)
        if chromosome > 22 or start > end:
            fail(f"invalid autosomal interval for locus {identity}")

    locus_by_id = {row["locus_id"]: row for row in loci}
    shared_path = root / "results/atlas/shared_loci.tsv"
    _, shared = read_tsv(shared_path)
    shared_primary = [row for row in shared if row.get("analysis_tier") == "PRIMARY_PHASE1"]
    shared_by_id = {row.get("shared_locus_id", ""): row for row in shared_primary}
    if "" in shared_by_id or len(shared_by_id) != len(shared_primary) or set(shared_by_id) != locus_ids:
        fail("loci.tsv is not the exact primary cross-method shared-locus family")
    for identity, locus in locus_by_id.items():
        upstream = shared_by_id[identity]
        if (
            locus["sleep_trait"] != upstream.get("sleep_trait")
            or locus["non_sleep_trait"] != upstream.get("non_sleep_trait")
            or locus["analysis_tier"] != upstream.get("analysis_tier")
            or locus["chromosome"] != upstream.get("CHR")
            or locus["start_bp"] != upstream.get("START")
            or locus["end_bp"] != upstream.get("STOP")
        ):
            fail(f"canonical locus fields drifted from shared_loci.tsv for {identity}")

    variants = tables["variants.tsv"]
    variant_locus = {(row["locus_id"], row["variant_id"]) for row in variants}
    variant_ids = {row["variant_id"] for row in variants}
    loci_with_variants = set()
    for row in variants:
        identity = f"{row['locus_id']}/{row['variant_id']}"
        if row["locus_id"] not in locus_ids:
            fail(f"orphan variant {identity}")
        chromosome = integer(row["chromosome"], "chromosome", identity, 1)
        position = integer(row["position_bp"], "position_bp", identity, 1)
        locus = locus_by_id[row["locus_id"]]
        if (
            chromosome > 22
            or chromosome != int(locus["chromosome"])
            or not int(locus["start_bp"]) <= position <= int(locus["end_bp"])
            or row["effect_allele"] not in {"A", "C", "G", "T"}
            or row["other_allele"] not in {"A", "C", "G", "T"}
            or row["effect_allele"] == row["other_allele"]
        ):
            fail(f"invalid allele/coordinate record for {identity}")
        loci_with_variants.add(row["locus_id"])
    primary_loci = {row["locus_id"] for row in loci if row["analysis_tier"] == "PRIMARY_PHASE1"}
    if not primary_loci or not primary_loci.issubset(loci_with_variants):
        fail("not every primary cross-method shared locus has fine-mapped variants")

    genes = tables["genes.tsv"]
    gene_locus = {(row["locus_id"], row["gene_id"]) for row in genes}
    gene_ids = {row["gene_id"] for row in genes}
    loci_with_genes = set()
    for row in genes:
        identity = f"{row['locus_id']}/{row['gene_id']}"
        if row["locus_id"] not in locus_ids:
            fail(f"orphan gene {identity}")
        stream_count = integer(row["evidence_stream_count"], "evidence_stream_count", identity, 1)
        required(row["claim_limit"], "claim_limit", identity)
        for field in ("eqtl_status", "sqtl_status", "pqtl_pwas_status", "twas_status"):
            if row[field] not in policy["molecular_integration"]["allowed_outcomes"]:
                fail(f"invalid {field} for {identity}")
        statuses = [row[field] for field in ("eqtl_status", "sqtl_status", "pqtl_pwas_status", "twas_status")]
        if stream_count != statuses.count("SUPPORTED"):
            fail(f"evidence_stream_count differs from supported modalities for {identity}")
        if row["colocalization_status"] not in policy["molecular_integration"]["allowed_outcomes"]:
            fail(f"invalid colocalization_status for {identity}")
        loci_with_genes.add(row["locus_id"])

    coverage_path = root / args.molecular_coverage
    coverage_fields, molecular_coverage = read_tsv(coverage_path)
    if coverage_fields != MOLECULAR_COVERAGE_FIELDS:
        fail("molecular_locus_coverage.tsv header differs from the frozen schema")
    unique_keys("molecular_locus_coverage.tsv", molecular_coverage, ["coverage_id"])
    required_modalities = set(policy["molecular_integration"]["required_modalities"])
    observed_coverage: set[tuple[str, str]] = set()
    for row in molecular_coverage:
        identity = row["coverage_id"]
        if row["locus_id"] not in locus_ids:
            fail(f"orphan molecular coverage row {identity}")
        if row["modality"] not in required_modalities:
            fail(f"unknown molecular coverage modality for {identity}")
        required_units = integer(row["required_unit_count"], "required_unit_count", identity, 1)
        completed_units = integer(row["completed_unit_count"], "completed_unit_count", identity, 0)
        component_counts = [integer(row[field], field, identity, 0) for field in (
            "analyzed_unit_count", "no_evidence_unit_count", "access_blocked_unit_count",
            "not_applicable_unit_count",
        )]
        if (
            row["coverage_status"] != "COMPLETE" or completed_units != required_units
            or any(value > completed_units for value in component_counts)
        ):
            fail(f"incomplete molecular coverage accounting for {identity}")
        evidence_path = root / row["evidence_path"]
        if not evidence_path.is_file() or not SHA256.fullmatch(row["evidence_sha256"]) or sha256(evidence_path) != row["evidence_sha256"]:
            fail(f"molecular coverage evidence differs for {identity}")
        observed_coverage.add((row["locus_id"], row["modality"]))
    expected_coverage = {(locus_id, modality) for locus_id in primary_loci for modality in required_modalities}
    if observed_coverage != expected_coverage:
        fail("molecular coverage does not span every primary locus and required modality")

    interpretation_manifest_path = root / args.interpretation_manifest
    interpretation_lock_path = root / args.interpretation_manifest_lock
    interpretation_coverage_path = root / args.interpretation_coverage
    interpretation_fields, interpretation_tasks = read_tsv(interpretation_manifest_path)
    interpretation_lock = json.loads(interpretation_lock_path.read_text(encoding="utf-8"))
    if (
        interpretation_fields != interpretation_policy["task_manifest_fields"]
        or interpretation_lock.get("policy_sha256") != sha256(interpretation_policy_path)
        or interpretation_lock.get("task_manifest_sha256") != sha256(interpretation_manifest_path)
        or interpretation_lock.get("task_ids_in_locked_order")
        != [row["task_id"] for row in interpretation_tasks]
        or interpretation_lock.get("interpretation_results_accessed_before_task_lock") is not False
    ):
        fail("interpretation task family differs from its pre-result lock")
    interpretation_coverage_fields, interpretation_coverage = read_tsv(interpretation_coverage_path)
    if interpretation_coverage_fields != interpretation_policy["coverage_fields"]:
        fail("interpretation_coverage.tsv header differs from the frozen policy")
    unique_keys("interpretation_coverage.tsv", interpretation_coverage, ["coverage_id"])
    expected_interpretation_coverage = {"ICOV__" + row["task_id"] for row in interpretation_tasks}
    if {row["coverage_id"] for row in interpretation_coverage} != expected_interpretation_coverage:
        fail("interpretation coverage is not the exact locked task family")
    for row in interpretation_coverage:
        identity = row["coverage_id"]
        required_units = integer(row["required_unit_count"], "required_unit_count", identity, 1)
        completed_units = integer(row["completed_unit_count"], "completed_unit_count", identity, 0)
        components = [integer(row[field], field, identity, 0) for field in (
            "analyzed_unit_count", "no_evidence_unit_count", "access_blocked_unit_count",
            "not_applicable_unit_count",
        )]
        if (
            row["coverage_status"] != "COMPLETE" or completed_units != required_units
            or sum(components) != completed_units
        ):
            fail(f"incomplete interpretation coverage accounting for {identity}")
        if components[2] != 0:
            fail(f"access-blocked interpretation task remains unresolved: {identity}")
        evidence_path = root / row["evidence_path"]
        if (
            not evidence_path.is_file() or not SHA256.fullmatch(row["evidence_sha256"])
            or sha256(evidence_path) != row["evidence_sha256"]
        ):
            fail(f"interpretation coverage evidence differs for {identity}")
    interpretation_provenance_path = root / args.interpretation_provenance
    interpretation_provenance = json.loads(interpretation_provenance_path.read_text(encoding="utf-8"))
    if (
        interpretation_provenance.get("policy_sha256") != sha256(interpretation_policy_path)
        or interpretation_provenance.get("task_manifest_sha256") != sha256(interpretation_manifest_path)
        or interpretation_provenance.get("outputs", {}).get(args.interpretation_coverage)
        != sha256(interpretation_coverage_path)
    ):
        fail("interpretation aggregate provenance differs from its locked inputs/coverage")
    for name in ("regulatory_elements.tsv", "cell_types.tsv", "pathways.tsv", "causal_tests.tsv"):
        relative = f"results/atlas/{name}"
        if interpretation_provenance.get("outputs", {}).get(relative) != sha256(atlas_dir / name):
            fail(f"{name} differs from interpretation aggregate provenance")

    regulatory = tables["regulatory_elements.tsv"]
    regulatory_ids = {row["regulatory_element_id"] for row in regulatory}
    for row in regulatory:
        identity = row["regulatory_element_id"]
        if (row["locus_id"], row["variant_id"]) not in variant_locus:
            fail(f"regulatory element {identity} references an absent locus/variant")
        if (row["locus_id"], row["target_gene_id"]) not in gene_locus:
            fail(f"regulatory element {identity} references an absent locus/gene")
    required_layers = set(policy["regulatory_mapping"]["required_layers"])
    required_contexts = set(policy["regulatory_mapping"]["required_context_domains"])
    if any(row["element_type"] not in required_layers for row in regulatory):
        fail("regulatory_elements.tsv contains an unlocked regulatory layer")
    if any(row["context_domain"] not in required_contexts for row in regulatory):
        fail("regulatory_elements.tsv contains an unlocked context domain")

    node_ids = traits | locus_ids | variant_ids | gene_ids | regulatory_ids
    for name in ("cell_types.tsv", "pathways.tsv"):
        for row in tables[name]:
            identity = row[definitions[name]["primary_key"][0]]
            if row["trait_or_locus_id"] not in traits | locus_ids:
                fail(f"{name} row {identity} references an absent trait/locus")
            if name == "cell_types.tsv":
                node_ids.add(row["cell_type_id"])
            else:
                node_ids.add(row["pathway_id"])
    if any(row["method"] not in set(policy["cell_types"]["required_strategies"]) for row in tables["cell_types.tsv"]):
        fail("cell_types.tsv contains an unlocked method")
    if any(row["domain"] not in set(policy["cell_types"]["required_domains"]) for row in tables["cell_types.tsv"]):
        fail("cell_types.tsv contains an unlocked biological domain")
    if any(row["resource"] not in set(policy["pathways"]["required_resources"]) for row in tables["pathways.tsv"]):
        fail("pathways.tsv contains an unlocked resource")

    for row in tables["causal_tests.tsv"]:
        identity = row["causal_test_id"]
        if row["exposure"] not in traits or row["outcome"] not in traits:
            fail(f"causal test {identity} references an absent trait")
        instruments = integer(row["n_instruments"], "n_instruments", identity, 1)
        minimum_f = number(row["f_statistic_min"], "f_statistic_min", identity)
        if instruments < policy["causal_inference"]["minimum_instruments_for_primary_multi_instrument_inference"]:
            fail(f"causal test {identity} has too few primary instruments")
        if minimum_f < policy["causal_inference"]["minimum_F_statistic"]:
            fail(f"causal test {identity} has a weak instrument")
    allowed_causal_methods = set(policy["causal_inference"]["methods"][:4]) | {"CAUSE", "LHC_MR"}
    if any(row["method"] not in allowed_causal_methods for row in tables["causal_tests.tsv"]):
        fail("causal_tests.tsv contains an unlocked estimator")
    if any(row["direction"] not in set(policy["causal_inference"]["directions"]) for row in tables["causal_tests.tsv"]):
        fail("causal_tests.tsv contains an unlocked direction")

    for row in tables["edges.tsv"]:
        identity = row["edge_id"]
        if row["source"] not in node_ids or row["target"] not in node_ids:
            fail(f"edge {identity} references an absent node")
        required(row["relationship"], "relationship", identity)
        required(row["method"], "method", identity)
        required(row["evidence_level"], "evidence_level", identity)
        required(row["dataset"], "dataset", identity)
        required(row["version"], "version", identity)

    if not args.quiet:
        print(
            "ATLAS_SCHEMA_OK "
            + " ".join(f"{name}={len(tables[name])}" for name in expected_names)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
