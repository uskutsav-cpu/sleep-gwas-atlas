#!/usr/bin/env python3
"""Validate and integrate locked molecular-QTL colocalization and corrected TWAS."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from liftover_chain import load_chain
import molecular_contract


COLOC_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset1_id",
    "dataset2_id", "molecular_feature_id", "gene_id", "gene_symbol", "modality",
    "tissue_cell_context", "coloc_method", "p1", "p2", "p12", "prior_role",
    "signal1", "signal2", "hit1", "hit2", "nsnps", "PP_H0", "PP_H1", "PP_H2",
    "PP_H3", "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant",
    "top_shared_variant_PP_H4", "fine_mapping_qc", "prior_robust", "analysis_status",
    "single_signal_fallback_justification", "source_family_id", "source_release",
    "provenance_id", "claim_limit",
]
ENGINE_COLOC_FIELDS = [
    "comparison_id", "pair_id", "locus_id", "comparison_type", "dataset1_id",
    "dataset2_id", "molecular_feature_id", "tissue_cell_context", "coloc_method",
    "p1", "p2", "p12", "prior_role", "signal1", "signal2", "hit1", "hit2",
    "nsnps", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "PP_H4_over_PP_H3",
    "top_shared_variant", "top_shared_variant_PP_H4", "fine_mapping_qc",
    "analysis_status", "single_signal_fallback_justification", "claim_limit",
]
EVIDENCE_FIELDS = [
    "evidence_id", "locus_id", "pair_id", "trait_id", "gene_id", "gene_symbol",
    "modality", "context", "source_family_id", "dataset_id", "feature_id",
    "evidence_type", "zscore", "p_value", "fdr", "pp_h4", "pp_h4_over_pp_h3",
    "prior_robust", "analysis_status", "status_reason", "evidence_path",
    "evidence_sha256", "provenance_id", "claim_limit",
]
COVERAGE_FIELDS = [
    "coverage_id", "locus_id", "pair_id", "trait_id", "modality",
    "source_family_id", "required_unit_count", "completed_unit_count",
    "analyzed_unit_count", "no_evidence_unit_count", "access_blocked_unit_count",
    "not_applicable_unit_count", "coverage_status", "evidence_path",
    "evidence_sha256", "provenance_id",
]
GENE_FIELDS = [
    "gene_id", "gene_symbol", "locus_id", "evidence_stream_count", "eqtl_status",
    "sqtl_status", "pqtl_pwas_status", "twas_status", "colocalization_status",
    "context", "evidence_level", "claim_limit", "provenance_id",
]
TWAS_FIELDS = molecular_contract.TWAS_RESULT_FIELDS
TWAS_COVERAGE_FIELDS = molecular_contract.TWAS_COVERAGE_FIELDS
GENE_ATTRIBUTE = re.compile(r'(?:^|;\s*)([A-Za-z_]+) "([^"]*)"')
MISSING = {"", "NA"}


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
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def versionless_gene(value: str) -> str:
    value = value.strip()
    return value.split(".", 1)[0] if value not in MISSING else "NA"


def normalize_modality(value: str) -> str:
    normalized = value.strip().lower().replace("-", "_")
    if normalized in {"eqtl", "ge"}:
        return "eQTL"
    if normalized in {"sqtl", "leafcutter"}:
        return "sQTL"
    if normalized in {"pqtl", "pwas", "aptamer", "pqtl_or_pwas"}:
        return "pQTL_or_PWAS"
    if normalized == "twas":
        return "TWAS"
    fail(f"unsupported molecular modality: {value}")
    raise AssertionError


def coverage_modalities(value: str) -> list[str]:
    if value.strip().lower().replace("-", "_") == "eqtl_or_sqtl":
        return ["eQTL", "sQTL"]
    return [normalize_modality(value)]


def probability(value: str, field: str, identity: str, missing_ok: bool = True) -> float | None:
    if value in MISSING:
        if missing_ok:
            return None
        fail(f"missing {field} for {identity}")
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"nonnumeric {field} for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        fail(f"invalid {field} for {identity}")
    return parsed


def ratio(value: str, identity: str) -> float | None:
    if value in MISSING:
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"nonnumeric PP_H4_over_PP_H3 for {identity}")
        raise AssertionError from exc
    if math.isnan(parsed) or parsed < 0:
        fail(f"invalid PP_H4_over_PP_H3 for {identity}")
    return parsed


def robust_signal_pairs(
    rows: list[dict[str, str]], grid: list[float], minimum_h4: float, minimum_ratio: float,
) -> set[tuple[str, str]]:
    """Return signal pairs that pass at every p12; never pool different signals."""
    expected = {float(value) for value in grid}
    if rows and {float(row["p12"]) for row in rows} != expected:
        fail("colocalization rows do not cover the exact locked prior grid")
    by_pair: dict[tuple[str, str], dict[float, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        key = (row["signal1"], row["signal2"])
        p12 = float(row["p12"])
        if p12 in by_pair[key]:
            fail(f"duplicate colocalization signal/prior row: {key}/{p12}")
        by_pair[key][p12] = row
    robust: set[tuple[str, str]] = set()
    for key, values in by_pair.items():
        if set(values) != expected:
            continue
        passing = True
        for row in values.values():
            h4 = probability(row["PP_H4"], "PP_H4", f"{key}")
            h4_h3 = ratio(row["PP_H4_over_PP_H3"], f"{key}")
            passing &= (
                row["coloc_method"] == "COLOC_SUSIE"
                and row["fine_mapping_qc"] == "PASS"
                and h4 is not None and h4 >= minimum_h4
                and h4_h3 is not None and h4_h3 >= minimum_ratio
            )
        if passing:
            robust.add(key)
    return robust


def select_primary_row(
    rows: list[dict[str, str]], primary_p12: float, robust: set[tuple[str, str]],
) -> dict[str, str] | None:
    candidates = [row for row in rows if float(row["p12"]) == primary_p12]
    if rows and not candidates:
        fail("primary p12 colocalization row is absent")
    def key(row: dict[str, str]) -> tuple[int, float, str, str]:
        pair = (row["signal1"], row["signal2"])
        value = probability(row["PP_H4"], "PP_H4", row["comparison_id"])
        return (int(pair in robust), value if value is not None else -1.0, row["signal1"], row["signal2"])
    return max(candidates, key=key) if candidates else None


def parse_gene_annotation(path: Path) -> dict[str, dict[str, object]]:
    genes: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chromosome = fields[0][3:] if fields[0].lower().startswith("chr") else fields[0]
            if not chromosome.isdigit() or not 1 <= int(chromosome) <= 22:
                continue
            attributes = dict(GENE_ATTRIBUTE.findall(fields[8]))
            gene = versionless_gene(attributes.get("gene_id", "NA"))
            if gene == "NA" or gene in genes:
                fail(f"missing or duplicate versionless gene at annotation line {line_number}: {gene}")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                fail(f"invalid gene coordinates at annotation line {line_number}")
                raise AssertionError from exc
            if start <= 0 or end < start:
                fail(f"invalid gene interval at annotation line {line_number}")
            genes[gene] = {
                "gene_id": gene, "gene_symbol": attributes.get("gene_name", "NA") or "NA",
                "chromosome": int(chromosome), "start": start, "end": end,
            }
    if not genes:
        fail("gene annotation contains no autosomal genes")
    return genes


def overlap_gene_loci(
    genes: dict[str, dict[str, object]], loci: list[dict[str, str]], chain,
) -> tuple[dict[str, list[str]], dict[str, tuple[int, int, int]]]:
    """Map loci to GRCh38 and return exact interval-overlap genes."""
    by_chromosome: dict[int, list[dict[str, object]]] = defaultdict(list)
    for gene in genes.values():
        by_chromosome[int(gene["chromosome"])].append(gene)
    overlap: dict[str, list[str]] = {}
    mapped_loci: dict[str, tuple[int, int, int]] = {}
    for locus in loci:
        chromosome, start, end = int(locus["chromosome"]), int(locus["start_bp"]), int(locus["end_bp"])
        start_status, mapped_start = chain.map_point(chromosome, start)
        end_status, mapped_end = chain.map_point(chromosome, end)
        if start_status != "mapped" or end_status != "mapped" or mapped_start is None or mapped_end is None:
            fail(f"locus endpoints do not both map to GRCh38: {locus['locus_id']}")
        if mapped_start[0] != mapped_end[0] or mapped_start[2] != mapped_end[2] or mapped_start[1] > mapped_end[1]:
            fail(f"locus endpoint mapping is not a single forward GRCh38 interval: {locus['locus_id']}")
        target_chr, target_start, target_end = mapped_start[0], mapped_start[1], mapped_end[1]
        mapped_loci[locus["locus_id"]] = (target_chr, target_start, target_end)
        overlap[locus["locus_id"]] = sorted(
            str(gene["gene_id"]) for gene in by_chromosome[target_chr]
            if int(gene["start"]) <= target_end and int(gene["end"]) >= target_start
        )
    return overlap, mapped_loci


def checked_output_hash(provenance: dict[str, object], relative: str, path: Path) -> None:
    outputs = provenance.get("outputs", {})
    if not isinstance(outputs, dict) or outputs.get(relative) != sha256(path):
        fail(f"artifact differs from provenance: {relative}")


def status_from_coverage(rows: list[dict[str, object]]) -> str:
    if any(int(row["analyzed_unit_count"]) > 0 or int(row["no_evidence_unit_count"]) > 0 for row in rows):
        return "NO_EVIDENCE_FOUND"
    if any(int(row["access_blocked_unit_count"]) > 0 for row in rows):
        return "ACCESS_BLOCKED"
    return "NOT_APPLICABLE"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--feature-manifest", default="results/tables/molecular_feature_manifest.tsv")
    parser.add_argument("--feature-lock", default="results/tables/molecular_feature_manifest.lock.json")
    parser.add_argument("--search-coverage", default="results/tables/molecular_search_coverage.tsv")
    parser.add_argument("--search-plan", default="results/tables/molecular_search_plan.tsv")
    parser.add_argument("--search-plan-lock", default="results/tables/molecular_search_plan.lock.json")
    parser.add_argument("--molecular-run-dir", default="results/molecular/coloc_runs")
    parser.add_argument("--trait-coloc", default="results/tables/trait_trait_colocalization.tsv")
    parser.add_argument("--fine-mapping-provenance", default="results/atlas/fine_mapping.provenance.json")
    parser.add_argument("--twas", default="results/tables/twas.tsv")
    parser.add_argument("--twas-coverage", default="results/tables/twas_coverage.tsv")
    parser.add_argument("--twas-provenance", default="results/tables/twas.provenance.json")
    parser.add_argument("--twas-models", default="results/tables/twas_model_registry.tsv")
    parser.add_argument("--twas-models-lock", default="results/tables/twas_model_registry.lock.json")
    parser.add_argument("--twas-phi-exclusions", default="results/tables/twas_model_phi_exclusions.tsv")
    parser.add_argument("--loci", default="results/atlas/loci.tsv")
    parser.add_argument("--traits", default="results/atlas/traits.tsv")
    parser.add_argument("--molecular-coloc-out", default="results/tables/molecular_qtl_colocalization.tsv")
    parser.add_argument("--all-coloc-out", default="results/tables/colocalization.tsv")
    parser.add_argument("--coverage-out", default="results/tables/molecular_locus_coverage.tsv")
    parser.add_argument("--evidence-out", default="results/tables/molecular_evidence.tsv")
    parser.add_argument("--genes-out", default="results/atlas/genes.tsv")
    parser.add_argument("--provenance-out", default="results/atlas/molecular.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, policy = molecular_contract.load_policy(root, args.policy)
    claim_limit = policy["claim_limit"]
    required_modalities = policy["required_modalities"]
    if required_modalities != ["eQTL", "sQTL", "pQTL_or_PWAS", "TWAS"]:
        fail("molecular integration modalities differ from the locked family")

    paths = {name: root / getattr(args, name) for name in (
        "feature_manifest", "feature_lock", "search_coverage", "search_plan",
        "search_plan_lock", "trait_coloc", "fine_mapping_provenance", "twas",
        "twas_coverage", "twas_provenance", "twas_models", "twas_models_lock",
        "twas_phi_exclusions", "loci", "traits",
    )}
    loci, _, _, fine_mapping_provenance = molecular_contract.validate_fine_mapping(root)
    _, traits = read_tsv(paths["traits"])
    if any(row["analysis_tier"] != "PRIMARY_PHASE1" for row in loci):
        fail("molecular integration locus family is not the primary fine-mapped family")
    locus_by_id = {row["locus_id"]: row for row in loci}
    trait_by_id = {row["trait_id"]: row for row in traits}
    if len(traits) != 45 or len(trait_by_id) != 45 or len(locus_by_id) != len(loci):
        fail("locus or trait identities are incomplete or duplicated")
    for locus in loci:
        if locus["sleep_trait"] not in trait_by_id or locus["non_sleep_trait"] not in trait_by_id:
            fail(f"locus references a trait outside the locked panel: {locus['locus_id']}")

    annotation_asset = next(
        (asset for asset in policy["metadata_assets"] if asset["id"] == policy["twas"]["gene_annotation_asset_id"]),
        None,
    )
    if annotation_asset is None:
        fail("TWAS gene annotation asset is not in the locked metadata family")
    annotation_path = root / annotation_asset["path"]
    if (
        not annotation_path.is_file() or annotation_path.stat().st_size != annotation_asset["bytes"]
        or sha256(annotation_path) != annotation_asset["sha256"]
    ):
        fail("GTEx v8 GENCODE v26 gene annotation differs from its exact pin")
    chain_spec = policy["reference_build"]
    chain, chain_provenance = load_chain(
        root / chain_spec["chain_path"], chain_spec["chain_sha256"], chain_spec["chain_bytes"],
    )
    genes = parse_gene_annotation(annotation_path)
    genes_by_locus, mapped_loci = overlap_gene_loci(genes, loci, chain)
    loci_by_gene_trait: dict[tuple[str, str], list[str]] = defaultdict(list)
    for locus in loci:
        for gene in genes_by_locus[locus["locus_id"]]:
            loci_by_gene_trait[(gene, locus["sleep_trait"])].append(locus["locus_id"])
            loci_by_gene_trait[(gene, locus["non_sleep_trait"])].append(locus["locus_id"])

    # Validate the exact molecular search/feature family before reading a result.
    plan, _ = molecular_contract.validate_search_plan(
        root, paths["search_plan"], paths["search_plan_lock"], policy_path,
        root / "results/tables/molecular_preflight.json",
    )
    feature_manifest, feature_lock = molecular_contract.validate_feature_family(
        root, paths["feature_manifest"], paths["feature_lock"], policy_path,
        root / "results/tables/molecular_preflight.json",
    )
    _, search_coverage = molecular_contract.read_tsv(
        paths["search_coverage"], molecular_contract.SEARCH_COVERAGE_FIELDS,
    )
    plan_by_task = {row["search_task_id"]: row for row in plan}

    primary_p12 = float(policy["colocalization"]["p12_primary"])
    p12_grid = [float(value) for value in policy["colocalization"]["p12_sensitivity"]]
    minimum_h4 = float(policy["colocalization"]["shared_signal_min_pp_h4"])
    minimum_h4_h3 = float(policy["colocalization"]["shared_signal_min_pp_h4_over_pp_h3"])
    molecular_coloc: list[dict[str, object]] = []
    evidence: list[dict[str, object]] = []
    support_by_gene: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    context_by_gene: dict[tuple[str, str], set[str]] = defaultdict(set)
    symbols_by_gene: dict[str, set[str]] = defaultdict(set)
    qtl_run_hashes: dict[str, str] = {}
    for feature in feature_manifest:
        comparison = feature["comparison_id"]
        run = root / args.molecular_run_dir / comparison
        provenance_path = run / "provenance.json"
        if not provenance_path.is_file():
            fail(f"molecular colocalization run is absent: {comparison}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        task_path = root / "results/molecular/tasks" / f"{comparison}.tsv"
        task_lock_path = root / "results/molecular/tasks" / f"{comparison}.lock.json"
        task_lock = json.loads(task_lock_path.read_text(encoding="utf-8"))
        if (
            provenance.get("schema_version") != "atlas-v1.0-molecular-coloc-run.2"
            or provenance.get("comparison_id") != comparison
            or provenance.get("policy_sha256") != sha256(policy_path)
            or provenance.get("task_sha256") != sha256(task_path)
            or provenance.get("task_lock_sha256") != sha256(task_lock_path)
            or provenance.get("manifest_sha256") != sha256(paths["feature_manifest"])
            or provenance.get("manifest_lock_sha256") != sha256(paths["feature_lock"])
            or provenance.get("script_sha256") != molecular_contract.script_hashes(root, "qtl")
            or task_lock.get("schema_version") != "atlas-v1.0-molecular-coloc-task.2"
            or task_lock.get("task_sha256") != sha256(task_path)
            or task_lock.get("script_sha256") != molecular_contract.script_hashes(root, "qtl")
        ):
            fail(f"molecular colocalization provenance differs: {comparison}")
        for filename, record in provenance.get("outputs", {}).items():
            output_path = run / filename
            if not output_path.is_file() or sha256(output_path) != record.get("sha256"):
                fail(f"molecular run output differs from provenance: {comparison}/{filename}")
        if set(provenance.get("outputs", {})) != set(molecular_contract.fine_mapping_contract.ENGINE_OUTPUT_FIELDS):
            fail(f"molecular colocalization output family is incomplete: {comparison}")
        diagnostics_fields, diagnostics = read_tsv(run / "diagnostics.tsv")
        if len(diagnostics) != 2 or any(row["model_converged"] != "TRUE" for row in diagnostics):
            fail(f"molecular SuSiE models are incomplete or nonconverged: {comparison}")
        coloc_fields, coloc_rows = read_tsv(run / "colocalization.tsv")
        if coloc_fields != ENGINE_COLOC_FIELDS or any(row["comparison_id"] != comparison for row in coloc_rows):
            fail(f"molecular colocalization schema or identity drifted: {comparison}")
        if any(
            provenance["outputs"][filename].get("rows")
            != len(read_tsv(run / filename)[1])
            for filename in molecular_contract.fine_mapping_contract.ENGINE_OUTPUT_FIELDS
        ):
            fail(f"molecular colocalization row counts differ from provenance: {comparison}")
        robust = robust_signal_pairs(coloc_rows, p12_grid, minimum_h4, minimum_h4_h3)
        selected = select_primary_row(coloc_rows, primary_p12, robust)
        gene_id = versionless_gene(feature["gene_id"])
        gene_symbol = feature["gene_symbol"] or "NA"
        if gene_id != "NA" and gene_symbol != "NA":
            symbols_by_gene[gene_id].add(gene_symbol)
        provenance_id = "MQC:" + sha256(provenance_path)
        for row in coloc_rows:
            pair = (row["signal1"], row["signal2"])
            normalized = {
                **row, "gene_id": gene_id, "gene_symbol": gene_symbol,
                "modality": normalize_modality(feature["modality"]),
                "prior_robust": "TRUE" if pair in robust else "FALSE",
                "source_family_id": feature["source_family_id"],
                "source_release": feature["exact_release"], "provenance_id": provenance_id,
            }
            molecular_coloc.append({field: normalized[field] for field in COLOC_FIELDS})
        status = "SUPPORTED" if robust else "NO_EVIDENCE_FOUND"
        selected_h4 = selected["PP_H4"] if selected else "NA"
        selected_ratio = selected["PP_H4_over_PP_H3"] if selected else "NA"
        evidence_path = run / "colocalization.tsv"
        evidence.append({
            "evidence_id": f"MQTL__{comparison}", "locus_id": feature["locus_id"],
            "pair_id": feature["pair_id"], "trait_id": feature["trait_id"],
            "gene_id": gene_id, "gene_symbol": gene_symbol,
            "modality": normalize_modality(feature["modality"]), "context": feature["context"],
            "source_family_id": feature["source_family_id"], "dataset_id": feature["dataset_id"],
            "feature_id": feature["feature_id"], "evidence_type": "QTL_COLOCALIZATION",
            "zscore": "NA", "p_value": "NA", "fdr": "NA", "pp_h4": selected_h4,
            "pp_h4_over_pp_h3": selected_ratio, "prior_robust": "TRUE" if robust else "FALSE",
            "analysis_status": status,
            "status_reason": "SAME_SIGNAL_PAIR_PASSES_ALL_LOCKED_PRIORS" if robust else "NO_PRIOR_ROBUST_SHARED_SIGNAL",
            "evidence_path": str(evidence_path.relative_to(root)), "evidence_sha256": sha256(evidence_path),
            "provenance_id": provenance_id, "claim_limit": claim_limit,
        })
        if status == "SUPPORTED" and gene_id != "NA":
            modality = normalize_modality(feature["modality"])
            support_by_gene[(feature["locus_id"], gene_id, modality)].add(comparison)
            context_by_gene[(feature["locus_id"], gene_id)].add(feature["context"])
        qtl_run_hashes[comparison] = sha256(provenance_path)

    # QTL source coverage is exact per locked locus/source family.
    qtl_coverage_groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in search_coverage:
        task = plan_by_task.get(row["search_task_id"])
        if task is None or row["locus_id"] != task["locus_id"] or row["source_family_id"] != task["source_family_id"]:
            fail(f"molecular search coverage differs from plan: {row['search_task_id']}")
        for modality in coverage_modalities(row["modality"]):
            qtl_coverage_groups[(row["locus_id"], row["source_family_id"], modality)].append(row)
    coverage: list[dict[str, object]] = []
    for (locus_id, source, modality), rows in sorted(qtl_coverage_groups.items()):
        locus = locus_by_id[locus_id]
        outcomes = [row["search_outcome"] for row in rows]
        analyzed = sum(value == "DATA_READY" for value in outcomes)
        no_evidence = sum(value in {"NO_FEATURES_IN_LOCUS", "NO_ANALYZABLE_FEATURES"} for value in outcomes)
        access_blocked = outcomes.count("ACCESS_BLOCKED")
        not_applicable = outcomes.count("NOT_APPLICABLE")
        if analyzed + no_evidence + access_blocked + not_applicable != len(rows):
            fail(f"invalid source terminal outcome at {locus_id}/{source}")
        evidence_path = paths["search_coverage"]
        coverage.append({
            "coverage_id": f"COV__{locus_id}__{source}__{modality}", "locus_id": locus_id,
            "pair_id": f"{locus['sleep_trait']}__{locus['non_sleep_trait']}", "trait_id": "ALL",
            "modality": modality, "source_family_id": source,
            "required_unit_count": len(rows), "completed_unit_count": len(rows),
            "analyzed_unit_count": analyzed, "no_evidence_unit_count": no_evidence,
            "access_blocked_unit_count": access_blocked, "not_applicable_unit_count": not_applicable,
            "coverage_status": "COMPLETE", "evidence_path": str(evidence_path.relative_to(root)),
            "evidence_sha256": sha256(evidence_path), "provenance_id": "MQLOCK:" + sha256(paths["feature_lock"]),
        })

    # Validate and normalize trait-trait colocalization into the same schema.
    checked_output_hash(fine_mapping_provenance, args.trait_coloc, paths["trait_coloc"])
    trait_fields, trait_coloc = read_tsv(paths["trait_coloc"])
    if trait_fields != ENGINE_COLOC_FIELDS:
        fail("trait-trait colocalization schema drifted")
    trait_coloc_by_comparison: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in trait_coloc:
        trait_coloc_by_comparison[row["comparison_id"]].append(row)
    trait_coloc_normalized: list[dict[str, object]] = []
    trait_provenance_id = "FM:" + sha256(paths["fine_mapping_provenance"])
    for comparison, rows in trait_coloc_by_comparison.items():
        robust = robust_signal_pairs(rows, p12_grid, minimum_h4, minimum_h4_h3)
        for row in rows:
            normalized = {
                **row, "gene_id": "NA", "gene_symbol": "NA", "modality": "TRAIT_TRAIT",
                "prior_robust": "TRUE" if (row["signal1"], row["signal2"]) in robust else "FALSE",
                "source_family_id": "ATLAS_FINE_MAPPING", "source_release": "atlas-v1.0",
                "provenance_id": trait_provenance_id,
            }
            trait_coloc_normalized.append({field: normalized[field] for field in COLOC_FIELDS})

    # Validate model and TWAS provenance, then stream the potentially large result table.
    models, model_lock = molecular_contract.validate_model_registry(
        root, paths["twas_models"], paths["twas_models_lock"], policy_path,
        paths["twas_phi_exclusions"],
    )
    model_by_id = {row["model_id"]: row for row in models}
    contexts_by_family: dict[str, set[str]] = defaultdict(set)
    for model in models:
        contexts_by_family[model["model_family"]].add(model["context"])
    expected_contexts = int(policy["twas"]["phi_model_source"]["expected_context_count"])
    if any(len(contexts) != expected_contexts for contexts in contexts_by_family.values()):
        fail("TWAS model families do not contain the exact locked context count")
    phi_fields, phi_exclusions = read_tsv(paths["twas_phi_exclusions"])
    del phi_fields
    if model_lock.get("phi_exclusions_sha256") != sha256(paths["twas_phi_exclusions"]):
        fail("TWAS phi exclusions differ from the model lock")

    twas_provenance = json.loads(paths["twas_provenance"].read_text(encoding="utf-8"))
    if (
        twas_provenance.get("schema_version") != "atlas-v1.0-twas-canonical.2"
        or twas_provenance.get("policy_sha256") != sha256(policy_path)
        or twas_provenance.get("script_sha256") != molecular_contract.script_hashes(root, "twas")
    ):
        fail("TWAS canonical provenance differs from the sealed production contract")
    checked_output_hash(twas_provenance, args.twas, paths["twas"])
    checked_output_hash(twas_provenance, args.twas_coverage, paths["twas_coverage"])
    twas_coverage_fields, twas_coverage = read_tsv(paths["twas_coverage"])
    if twas_coverage_fields != TWAS_COVERAGE_FIELDS:
        fail("TWAS coverage schema drifted")
    completed_twas: dict[tuple[str, str], dict[str, str]] = {}
    ineligible_traits: dict[str, dict[str, str]] = {}
    for row in twas_coverage:
        if row["analysis_status"] == "COMPLETED":
            key = (row["trait_id"], row["model_id"])
            if key in completed_twas or row["model_id"] not in model_by_id:
                fail(f"duplicate or unknown completed TWAS unit: {key}")
            if row["provenance_sha256"] != sha256(root / row["provenance_path"]):
                fail(f"TWAS run provenance differs: {row['run_id']}")
            completed_twas[key] = row
        elif row["analysis_status"] == "NOT_APPLICABLE" and row["model_id"] == "NA":
            if row["trait_id"] in ineligible_traits:
                fail(f"duplicate terminal TWAS trait: {row['trait_id']}")
            ineligible_traits[row["trait_id"]] = row
        else:
            fail(f"invalid TWAS coverage outcome: {row['trait_id']}/{row['model_id']}")
    expected_eligible_traits = {
        trait for trait in trait_by_id if trait not in ineligible_traits
    }
    expected_completed = {(trait, model) for trait in expected_eligible_traits for model in model_by_id}
    if set(completed_twas) != expected_completed:
        fail("TWAS coverage does not contain the exact eligible trait-by-model family")

    twas_unit_state: dict[tuple[str, str, str], dict[str, int]] = {}
    for locus in loci:
        for trait in (locus["sleep_trait"], locus["non_sleep_trait"]):
            for model_id in model_by_id:
                key = (locus["locus_id"], trait, model_id)
                twas_unit_state[key] = {"observed": 0, "supported": 0}
    twas_file_hash = sha256(paths["twas"])
    twas_coverage_hash = sha256(paths["twas_coverage"])
    seen_twas_units: set[tuple[str, str]] = set()
    current_twas_unit: tuple[str, str] | None = None
    current_unit_genes: set[str] = set()
    twas_row_count = 0
    with paths["twas"].open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if (reader.fieldnames or []) != TWAS_FIELDS:
            fail("TWAS aggregate schema drifted")
        for row in reader:
            twas_row_count += 1
            unit_key = (row["trait_id"], row["model_id"])
            if unit_key != current_twas_unit:
                if unit_key in seen_twas_units:
                    fail(f"noncontiguous repeated TWAS unit: {unit_key}")
                seen_twas_units.add(unit_key)
                current_twas_unit = unit_key
                current_unit_genes.clear()
            gene_id = versionless_gene(row["gene_id"])
            if gene_id in current_unit_genes:
                fail(f"duplicate TWAS gene within run: {unit_key}/{gene_id}")
            current_unit_genes.add(gene_id)
            if unit_key not in completed_twas:
                fail(f"TWAS row is outside the completed family: {row['twas_id']}")
            if row["status"] not in {"SUPPORTED", "NO_EVIDENCE_FOUND"}:
                fail(f"invalid TWAS result status: {row['twas_id']}")
            p_value = probability(row["p_value"], "p_value", row["twas_id"], missing_ok=False)
            fdr = probability(row["fdr"], "fdr", row["twas_id"], missing_ok=False)
            del p_value, fdr
            model = model_by_id[row["model_id"]]
            if row["context"] != model["context"] or row["model_family"] != model["model_family"]:
                fail(f"TWAS model metadata drifted: {row['twas_id']}")
            if gene_id != "NA" and row["gene_symbol"] not in MISSING:
                symbols_by_gene[gene_id].add(row["gene_symbol"])
            for locus_id in loci_by_gene_trait.get((gene_id, row["trait_id"]), []):
                locus = locus_by_id[locus_id]
                unit = twas_unit_state[(locus_id, row["trait_id"], row["model_id"])]
                unit["observed"] += 1
                unit["supported"] += int(row["status"] == "SUPPORTED")
                evidence.append({
                    "evidence_id": f"TWASLOC__{locus_id}__{row['twas_id']}", "locus_id": locus_id,
                    "pair_id": f"{locus['sleep_trait']}__{locus['non_sleep_trait']}",
                    "trait_id": row["trait_id"], "gene_id": gene_id,
                    "gene_symbol": genes.get(gene_id, {}).get("gene_symbol", row["gene_symbol"]),
                    "modality": "TWAS", "context": row["context"],
                    "source_family_id": row["model_family"], "dataset_id": row["model_id"],
                    "feature_id": gene_id, "evidence_type": "VARIANCE_CONTROLLED_SPREDIXCAN",
                    "zscore": row["zscore"], "p_value": row["p_value"], "fdr": row["fdr"],
                    "pp_h4": "NA", "pp_h4_over_pp_h3": "NA", "prior_robust": "NA",
                    "analysis_status": row["status"],
                    "status_reason": "BH_FDR_LE_0.05_WITHIN_TRAIT_MODEL_FAMILY" if row["status"] == "SUPPORTED" else "BH_FDR_GT_0.05_WITHIN_TRAIT_MODEL_FAMILY",
                    "evidence_path": args.twas, "evidence_sha256": twas_file_hash,
                    "provenance_id": row["provenance_id"], "claim_limit": claim_limit,
                })
                if row["status"] == "SUPPORTED":
                    support_by_gene[(locus_id, gene_id, "TWAS")].add(row["twas_id"])
                    context_by_gene[(locus_id, gene_id)].add(row["context"])
    if twas_row_count != int(twas_provenance.get("result_row_count", -1)):
        fail("TWAS aggregate row count differs from provenance")

    # Preserve gene-specific missing-phi exclusions when the excluded gene overlaps a locus.
    twas_hash = sha256(paths["twas_phi_exclusions"])
    for exclusion in phi_exclusions:
        model_id = exclusion["model_id"]
        if model_id not in model_by_id or exclusion["exclusion_reason"] != "MODEL_FEATURE_MISSING_PHI":
            fail("invalid TWAS model-feature phi exclusion")
        gene_id = versionless_gene(exclusion["gene_id"])
        for trait in expected_eligible_traits:
            for locus_id in loci_by_gene_trait.get((gene_id, trait), []):
                locus = locus_by_id[locus_id]
                model = model_by_id[model_id]
                evidence.append({
                    "evidence_id": f"TWASPHI__{locus_id}__{trait}__{model_id}__{gene_id}",
                    "locus_id": locus_id, "pair_id": f"{locus['sleep_trait']}__{locus['non_sleep_trait']}",
                    "trait_id": trait, "gene_id": gene_id,
                    "gene_symbol": genes.get(gene_id, {}).get("gene_symbol", "NA"), "modality": "TWAS",
                    "context": model["context"], "source_family_id": model["model_family"],
                    "dataset_id": model_id, "feature_id": gene_id, "evidence_type": "MODEL_ELIGIBILITY",
                    "zscore": "NA", "p_value": "NA", "fdr": "NA", "pp_h4": "NA",
                    "pp_h4_over_pp_h3": "NA", "prior_robust": "NA", "analysis_status": "NOT_APPLICABLE",
                    "status_reason": "MODEL_FEATURE_MISSING_PHI", "evidence_path": args.twas_phi_exclusions,
                    "evidence_sha256": twas_hash, "provenance_id": "TWMODEL:" + sha256(paths["twas_models_lock"]),
                    "claim_limit": claim_limit,
                })

    # Aggregate exact 49-context TWAS coverage for each trait represented at each locus.
    for locus in loci:
        locus_id = locus["locus_id"]
        pair_id = f"{locus['sleep_trait']}__{locus['non_sleep_trait']}"
        for trait in (locus["sleep_trait"], locus["non_sleep_trait"]):
            for family, contexts in sorted(contexts_by_family.items()):
                family_models = [model["model_id"] for model in models if model["model_family"] == family]
                if len(family_models) != len(contexts):
                    fail(f"duplicate TWAS contexts in model family: {family}")
                if trait in ineligible_traits:
                    analyzed = no_evidence = access_blocked = 0
                    not_applicable = len(family_models)
                else:
                    states = [twas_unit_state[(locus_id, trait, model_id)] for model_id in family_models]
                    analyzed = len(states)
                    no_evidence = sum(state["supported"] == 0 for state in states)
                    access_blocked = not_applicable = 0
                coverage.append({
                    "coverage_id": f"COV__{locus_id}__{trait}__{family}", "locus_id": locus_id,
                    "pair_id": pair_id, "trait_id": trait, "modality": "TWAS",
                    "source_family_id": family, "required_unit_count": len(family_models),
                    "completed_unit_count": len(family_models), "analyzed_unit_count": analyzed,
                    "no_evidence_unit_count": no_evidence, "access_blocked_unit_count": access_blocked,
                    "not_applicable_unit_count": not_applicable, "coverage_status": "COMPLETE",
                    "evidence_path": args.twas_coverage, "evidence_sha256": twas_coverage_hash,
                    "provenance_id": "TWAS:" + sha256(paths["twas_provenance"]),
                })

    # Require every locus to have complete coverage for all four modalities.
    coverage_by_locus_modality: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in coverage:
        coverage_by_locus_modality[(str(row["locus_id"]), str(row["modality"]))].append(row)
    expected_coverage_keys = {(locus_id, modality) for locus_id in locus_by_id for modality in required_modalities}
    if set(coverage_by_locus_modality) != expected_coverage_keys:
        missing = sorted(expected_coverage_keys - set(coverage_by_locus_modality))
        extra = sorted(set(coverage_by_locus_modality) - expected_coverage_keys)
        fail(f"molecular locus coverage differs from exact locus-by-modality family: missing={missing[:5]} extra={extra[:5]}")
    if any(row["coverage_status"] != "COMPLETE" for row in coverage):
        fail("molecular locus coverage contains an incomplete unit")

    # Build only genuinely supported genes; absence remains in coverage, never a fake gene row.
    molecular_input_digest = hashlib.sha256(json.dumps({
        "feature_lock": sha256(paths["feature_lock"]), "qtl_runs": qtl_run_hashes,
        "twas": sha256(paths["twas"]), "twas_provenance": sha256(paths["twas_provenance"]),
        "annotation": sha256(annotation_path), "chain": chain_provenance["chain_sha256"],
        "loci": sha256(paths["loci"]),
    }, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    molecular_provenance_id = "MOL:" + molecular_input_digest
    gene_rows: list[dict[str, object]] = []
    supported_gene_keys = sorted({(locus, gene) for locus, gene, _modality in support_by_gene})
    qtl_modalities = {"eQTL", "sQTL", "pQTL_or_PWAS"}
    for locus_id, gene_id in supported_gene_keys:
        status: dict[str, str] = {}
        for modality in required_modalities:
            if support_by_gene.get((locus_id, gene_id, modality)):
                status[modality] = "SUPPORTED"
                continue
            evaluated = any(
                row["locus_id"] == locus_id and row["gene_id"] == gene_id
                and row["modality"] == modality and row["analysis_status"] == "NO_EVIDENCE_FOUND"
                for row in evidence
            )
            if evaluated:
                status[modality] = "NO_EVIDENCE_FOUND"
            else:
                status[modality] = status_from_coverage(coverage_by_locus_modality[(locus_id, modality)])
        stream_count = sum(status[modality] == "SUPPORTED" for modality in required_modalities)
        if stream_count < 1:
            fail(f"unsupported gene reached canonical gene table: {locus_id}/{gene_id}")
        symbol = str(genes.get(gene_id, {}).get("gene_symbol", "NA"))
        if symbol == "NA" and symbols_by_gene[gene_id]:
            symbol = sorted(symbols_by_gene[gene_id])[0]
        qtl_evaluated = [status[modality] for modality in qtl_modalities]
        coloc_status = "SUPPORTED" if "SUPPORTED" in qtl_evaluated else (
            "NO_EVIDENCE_FOUND" if "NO_EVIDENCE_FOUND" in qtl_evaluated else
            "ACCESS_BLOCKED" if "ACCESS_BLOCKED" in qtl_evaluated else "NOT_APPLICABLE"
        )
        gene_rows.append({
            "gene_id": gene_id, "gene_symbol": symbol, "locus_id": locus_id,
            "evidence_stream_count": stream_count, "eqtl_status": status["eQTL"],
            "sqtl_status": status["sQTL"], "pqtl_pwas_status": status["pQTL_or_PWAS"],
            "twas_status": status["TWAS"], "colocalization_status": coloc_status,
            "context": ";".join(sorted(context_by_gene[(locus_id, gene_id)])) or "NA",
            "evidence_level": "HIGH_CONFIDENCE_CONVERGENT_GENE" if stream_count >= 2 else "SINGLE_STREAM_SUPPORTED_GENE",
            "claim_limit": claim_limit, "provenance_id": molecular_provenance_id,
        })

    molecular_coloc.sort(key=lambda row: (
        str(row["locus_id"]), str(row["comparison_id"]), float(str(row["p12"])),
        str(row["signal1"]), str(row["signal2"]),
    ))
    trait_coloc_normalized.sort(key=lambda row: (
        str(row["locus_id"]), str(row["comparison_id"]), float(str(row["p12"])),
        str(row["signal1"]), str(row["signal2"]),
    ))
    all_coloc = trait_coloc_normalized + molecular_coloc
    evidence.sort(key=lambda row: (str(row["locus_id"]), str(row["modality"]), str(row["trait_id"]), str(row["evidence_id"])))
    coverage.sort(key=lambda row: (str(row["locus_id"]), str(row["modality"]), str(row["trait_id"]), str(row["source_family_id"])))
    payloads = {
        root / args.molecular_coloc_out: table_text(COLOC_FIELDS, molecular_coloc),
        root / args.all_coloc_out: table_text(COLOC_FIELDS, all_coloc),
        root / args.coverage_out: table_text(COVERAGE_FIELDS, coverage),
        root / args.evidence_out: table_text(EVIDENCE_FIELDS, evidence),
        root / args.genes_out: table_text(GENE_FIELDS, gene_rows),
    }
    provenance = {
        "schema_version": "atlas-v1.0-molecular-integration.2", "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "feature_manifest_sha256": sha256(paths["feature_manifest"]),
        "feature_lock_sha256": sha256(paths["feature_lock"]), "search_plan_sha256": sha256(paths["search_plan"]),
        "search_plan_lock_sha256": sha256(paths["search_plan_lock"]), "qtl_run_provenance_sha256": qtl_run_hashes,
        "trait_trait_colocalization_sha256": sha256(paths["trait_coloc"]),
        "fine_mapping_provenance_sha256": sha256(paths["fine_mapping_provenance"]),
        "twas_sha256": sha256(paths["twas"]), "twas_provenance_sha256": sha256(paths["twas_provenance"]),
        "twas_model_registry_sha256": sha256(paths["twas_models"]),
        "twas_model_registry_lock_sha256": sha256(paths["twas_models_lock"]),
        "gene_annotation": {"path": annotation_asset["path"], "bytes": annotation_asset["bytes"], "sha256": annotation_asset["sha256"]},
        "liftover_chain": chain_provenance,
        "mapped_loci_grch38": {key: list(value) for key, value in sorted(mapped_loci.items())},
        "molecular_qtl_colocalization_rows": len(molecular_coloc), "all_colocalization_rows": len(all_coloc),
        "molecular_evidence_rows": len(evidence), "molecular_locus_coverage_rows": len(coverage),
        "supported_gene_rows": len(gene_rows), "molecular_provenance_id": molecular_provenance_id,
        "zero_locus_qtl_not_applicable": len(loci) == 0,
        "script_sha256": {
            "qtl": molecular_contract.script_hashes(root, "qtl"),
            "twas": molecular_contract.script_hashes(root, "twas"),
        },
        "outputs": {str(path.relative_to(root)): hashlib.sha256(text.encode()).hexdigest() for path, text in payloads.items()},
        "claim_limit": claim_limit,
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        for path, text in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                fail(f"molecular aggregate differs from deterministic recomputation: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            fail("molecular integration provenance drifted")
    else:
        if any(path.exists() for path in [*payloads, provenance_path]):
            fail("molecular integration outputs already exist; refusing overwrite")
        for path, text in payloads.items():
            atomic_text(path, text)
        atomic_text(provenance_path, provenance_text)
    if not args.quiet:
        print(
            f"MOLECULAR_INTEGRATION_OK qtl_coloc={len(molecular_coloc)} "
            f"evidence={len(evidence)} coverage={len(coverage)} genes={len(gene_rows)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
