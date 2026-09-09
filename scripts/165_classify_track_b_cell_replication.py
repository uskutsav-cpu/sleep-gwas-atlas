#!/usr/bin/env python3
"""Classify frozen Allen-to-GSE67835 broad-cell replication.

This consumes only a complete, contract-bound 7+7 MAGMA enrichment run.  It
never treats the GSE hybrid column as aligned, never treats the shared GSE
neurons column as two independent replications, and never emits a pericyte or
mouse-derived claim.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import stat
import sys
from pathlib import Path, PurePosixPath
from typing import Mapping, Sequence

sys.dont_write_bytecode = True

ROOT_AT_IMPORT = Path(__file__).resolve().parents[1]
ENRICHMENT_REL = "scripts/164_run_track_b_cell_class_enrichment.py"
ENRICHMENT_SHA256 = "9d5bb622b5c27c6f38b1cf6aff4402bdcfa847e748465266207f07e6dee06619"
SELF_REL = "scripts/165_classify_track_b_cell_replication.py"
REPLICATION_FIELDS = [
    "target_trait_id", "allen_class", "gse_class", "alignment_status",
    "replication_resolution", "discovery_status", "discovery_beta",
    "discovery_bh_fdr", "discovery_direction_pass", "replication_status",
    "replication_beta", "replication_bh_fdr", "replication_direction_pass",
    "discovery_eligible", "direction_concordant", "both_atlas_fdr_pass",
    "replication_classification", "independent_replication_credit", "claim_limit",
]


def _bootstrap_hash(path: Path) -> str:
    before_path = os.lstat(path)
    if stat.S_ISLNK(before_path.st_mode) or not stat.S_ISREG(before_path.st_mode):
        raise RuntimeError("enrichment implementation is not a real regular file")
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        digest = hashlib.sha256()
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_size,
        value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before_path) != identity(before) or identity(before) != identity(after) or identity(after) != identity(final):
        raise RuntimeError("enrichment implementation drifted during classifier bootstrap")
    return digest.hexdigest()


if _bootstrap_hash(ROOT_AT_IMPORT / ENRICHMENT_REL) != ENRICHMENT_SHA256:
    raise RuntimeError("enrichment implementation differs from the classifier bootstrap pin")
_SPEC = importlib.util.spec_from_file_location("track_b_cell_class_enrichment_v4", ROOT_AT_IMPORT / ENRICHMENT_REL)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError("cannot load the frozen enrichment implementation")
E = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E
_SPEC.loader.exec_module(E)


def _parse_float(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise E.ContractError(f"nonnumeric {label}") from exc
    if not math.isfinite(parsed):
        raise E.ContractError(f"nonfinite {label}")
    return parsed


def _validate_normalized_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    dataset_id: str,
    trait_id: str,
    classes: Sequence[str],
    minimum_gene_rows: int,
) -> None:
    if len(rows) != 7 or [row.get("class_id") for row in rows] != list(classes):
        raise E.ContractError(f"{dataset_id} does not retain its exact ordered seven-class family")
    p_for_bh: list[float] = []
    for row, class_id in zip(rows, classes):
        if (
            row.get("dataset_id") != dataset_id or row.get("target_trait_id") != trait_id
            or row.get("source_column") != class_id
            or row.get("direction_pass") not in {"YES", "NO"}
        ):
            raise E.ContractError(f"{dataset_id}/{class_id} normalized identity is invalid")
        status = str(row.get("analysis_status"))
        p_for = _parse_float(str(row.get("p_value_for_bh")), f"{dataset_id}/{class_id} BH input")
        if not 0 <= p_for <= 1:
            raise E.ContractError(f"{dataset_id}/{class_id} has an invalid BH input")
        p_for_bh.append(p_for)
        if status == "VALID":
            beta = _parse_float(str(row.get("beta")), f"{dataset_id}/{class_id} beta")
            se = _parse_float(str(row.get("se")), f"{dataset_id}/{class_id} SE")
            z_value = _parse_float(str(row.get("z")), f"{dataset_id}/{class_id} z")
            p_value = _parse_float(str(row.get("p_value")), f"{dataset_id}/{class_id} p")
            try:
                n_genes = int(str(row.get("n_genes")))
            except ValueError as exc:
                raise E.ContractError(f"{dataset_id}/{class_id} has an invalid gene count") from exc
            if (
                se <= 0 or not 0 <= p_value <= 1 or n_genes < minimum_gene_rows
                or p_for != p_value
                or not math.isclose(z_value, beta / se, rel_tol=1e-14, abs_tol=1e-14)
                or row.get("direction_pass") != ("YES" if beta > 0 else "NO")
            ):
                raise E.ContractError(f"{dataset_id}/{class_id} valid-row statistics are inconsistent")
        else:
            if not status.startswith("FAILED_") or p_for != 1.0 or row.get("direction_pass") != "NO":
                raise E.ContractError(f"{dataset_id}/{class_id} failed row entered inference")
    adjusted = E.bh_adjust(list(classes), p_for_bh)
    for row in rows:
        expected = format(adjusted[str(row["class_id"])], ".17g")
        if row.get("bh_fdr") != expected:
            raise E.ContractError(f"{dataset_id} BH correction is not the exact seven-class family")


def _directory_file_family(root: Path, relative: str, expected: Sequence[str]) -> None:
    relative = E.safe_relative(relative)
    E._check_parents(root, f"{relative}/sentinel")
    info = os.lstat(root / relative)
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise E.ContractError("enrichment input is not a real directory")
    observed = sorted(entry.name for entry in os.scandir(root / relative))
    if observed != sorted(expected):
        raise E.ContractError("enrichment directory has an incomplete or unknown file family")
    for name in expected:
        E.stable_file(root, f"{relative}/{name}")


def validate_enrichment_family(
    root: Path,
    enrichment_relative: str,
    policy: Mapping[str, object],
) -> tuple[dict[str, object], list[dict[str, str]], list[dict[str, str]], dict[str, E.FileRecord]]:
    production = policy.get("production_outputs")
    gate_policy = policy.get("upstream_gate")
    if not isinstance(production, dict) or not isinstance(gate_policy, dict):
        raise E.ContractError("V4 runtime policy lacks its production/gate sections")
    expected_files = production["enrichment_required_files"]
    if not isinstance(expected_files, list):
        raise E.ContractError("V4 enrichment file family is malformed")
    _directory_file_family(root, enrichment_relative, expected_files)
    provenance_rel = f"{enrichment_relative}/P14_P16_ENRICHMENT.provenance.json"
    provenance = E.read_canonical_json(root, provenance_rel)
    required = {
        "schema_version", "terminal_status", "pair_id", "target_trait_id",
        "upstream_gate_path", "upstream_gate_sha256", "eligible_family",
        "contract_lock_path", "contract_lock_sha256", "policy_path", "policy_sha256",
        "genes_raw_path", "genes_raw_sha256", "gene_results_provenance_path",
        "gene_results_provenance_sha256", "allen_aggregation", "runs",
        "multiple_testing", "valid_result_rows", "total_result_rows", "output_sha256",
        "claim_limit",
    }
    if (
        set(provenance) != required
        or provenance.get("schema_version") != "track-b-p14-p16-enrichment-run.1"
        or provenance.get("terminal_status") not in {
            "COMPLETED_ALL_14_ROWS_VALID", "COMPLETED_WITH_FAILURES_PRESERVED",
        }
        or provenance.get("total_result_rows") != 14
        or provenance.get("contract_lock_path") != E.CONTRACT_LOCK_REL
        or provenance.get("contract_lock_sha256") != E.stable_file(root, E.CONTRACT_LOCK_REL).sha256
        or provenance.get("policy_path") != E.POLICY_REL
        or provenance.get("policy_sha256") != E.stable_file(root, E.POLICY_REL).sha256
        or provenance.get("multiple_testing") != {
            "allen_denominator": 7, "gse_denominator": 7,
            "average_is_hypothesis": False, "failed_or_na_p_for_bh": 1.0,
        }
    ):
        raise E.ContractError("enrichment provenance schema/model is invalid")
    pair_id, trait_id = str(provenance.get("pair_id")), str(provenance.get("target_trait_id"))
    expected_dir = str(production["enrichment_dir_template"]).format(pair_id=pair_id, trait_id=trait_id)
    if enrichment_relative != expected_dir:
        raise E.ContractError("enrichment input is outside its exact pair/trait namespace")
    gate = E.validate_upstream_gate(root, str(provenance["upstream_gate_path"]), policy)
    if (
        gate["gate_sha256"] != provenance.get("upstream_gate_sha256")
        or gate["pair_id"] != pair_id or gate["target_trait_id"] != trait_id
        or gate["eligible_family"] != provenance.get("eligible_family")
    ):
        raise E.ContractError("enrichment provenance is not bound to the exact current upstream gate")
    genes_rel, genes_record, genes_prov_rel, genes_prov_record = E.validate_gene_results(
        root, trait_id, str(gate["target_source_id"]), policy,
    )
    if (
        genes_rel != provenance.get("genes_raw_path")
        or genes_record.sha256 != provenance.get("genes_raw_sha256")
        or genes_prov_rel != provenance.get("gene_results_provenance_path")
        or genes_prov_record.sha256 != provenance.get("gene_results_provenance_sha256")
    ):
        raise E.ContractError("enrichment provenance is not bound to current gene results")
    output_hashes = provenance.get("output_sha256")
    expected_hashed = set(expected_files) - {"P14_P16_ENRICHMENT.provenance.json"}
    if not isinstance(output_hashes, dict) or set(output_hashes) != expected_hashed:
        raise E.ContractError("enrichment provenance output family is incomplete")
    records: dict[str, E.FileRecord] = {}
    for name in expected_files:
        record = E.stable_file(root, f"{enrichment_relative}/{name}")
        records[name] = record
        if name != "P14_P16_ENRICHMENT.provenance.json" and record.sha256 != output_hashes[name]:
            raise E.ContractError(f"enrichment output differs from provenance: {name}")
    allen_rel = f"{enrichment_relative}/P14_ALLEN_BROAD_CLASS_RESULTS.tsv"
    gse_rel = f"{enrichment_relative}/P16_GSE67835_CLASS_RESULTS.tsv"
    allen_rows = E.parse_tsv_payload(E.read_stable_bytes(root, allen_rel), E.RESULT_FIELDS, "Allen enrichment")
    gse_rows = E.parse_tsv_payload(E.read_stable_bytes(root, gse_rel), E.RESULT_FIELDS, "GSE67835 enrichment")
    minimum_rows = int(gate_policy["minimum_gene_rows"])
    _validate_normalized_rows(
        allen_rows, dataset_id="Allen_Human_MTG_level2", trait_id=trait_id,
        classes=E.CLASS_ORDER, minimum_gene_rows=minimum_rows,
    )
    _validate_normalized_rows(
        gse_rows, dataset_id="GSE67835_Human_Cortex_woFetal", trait_id=trait_id,
        classes=E.GSE_ORDER, minimum_gene_rows=minimum_rows,
    )
    runs = provenance.get("runs")
    if (
        not isinstance(runs, list) or len(runs) != 2
        or any(not isinstance(row, dict) for row in runs)
        or [row.get("dataset_id") for row in runs] != [
            "Allen_Human_MTG_level2", "GSE67835_Human_Cortex_woFetal",
        ]
    ):
        raise E.ContractError("enrichment provenance lacks the exact two-atlas run family")
    run_by_id = {str(row["dataset_id"]): row for row in runs if isinstance(row, dict)}
    families = (
        (
            "Allen_Human_MTG_level2", E.CLASS_ORDER,
            "ARITHMETIC_MEAN_ALL_MAPPED_SUBTYPES_MATH_FSUM_SOURCE_ORDER",
            E.EXPECTED_CLASS_COUNTS, allen_rows,
        ),
        (
            "GSE67835_Human_Cortex_woFetal", E.GSE_ORDER,
            "SOURCE_CLASS_COLUMN_NO_AGGREGATION", {name: 1 for name in E.GSE_ORDER}, gse_rows,
        ),
    )
    specs = E._dataset_specs(policy)
    for dataset_id, classes, aggregation, counts, normalized in families:
        run = run_by_id[dataset_id]
        failure_code = str(run.get("failure_code"))
        if failure_code != "NONE" and not failure_code.startswith("FAILED_"):
            raise E.ContractError(f"invalid run failure state: {dataset_id}")
        if (
            run.get("result_rows") != 7
            or run.get("valid_rows") != sum(row["analysis_status"] == "VALID" for row in normalized)
            or run.get("source_text_sha256") != specs[dataset_id]["text_sha256"]
        ):
            raise E.ContractError(f"enrichment run metadata is inconsistent: {dataset_id}")
        if dataset_id == "Allen_Human_MTG_level2":
            aggregate = provenance.get("allen_aggregation")
            if not isinstance(aggregate, dict) or run.get("matrix_sha256") != aggregate.get("sha256"):
                raise E.ContractError("Allen aggregate identity is inconsistent")
        elif run.get("matrix_sha256") != specs[dataset_id]["text_sha256"]:
            raise E.ContractError("GSE matrix identity is inconsistent")
        reconstructed = E.normalize_gsa(
            root / enrichment_relative / f"{dataset_id}.gsa.out",
            classes, dataset_id, trait_id, aggregation, counts,
            minimum_gene_rows=minimum_rows,
            forced_failure=None if failure_code == "NONE" else failure_code,
        )
        if E.tsv_bytes(E.RESULT_FIELDS, reconstructed) != E.tsv_bytes(E.RESULT_FIELDS, normalized):
            raise E.ContractError(f"normalized enrichment rows do not reproduce raw MAGMA output: {dataset_id}")
    valid_rows = sum(row["analysis_status"] == "VALID" for row in [*allen_rows, *gse_rows])
    expected_terminal = "COMPLETED_ALL_14_ROWS_VALID" if valid_rows == 14 else "COMPLETED_WITH_FAILURES_PRESERVED"
    if provenance.get("valid_result_rows") != valid_rows or provenance.get("terminal_status") != expected_terminal:
        raise E.ContractError("enrichment terminal status does not match its preserved rows")
    return provenance, allen_rows, gse_rows, records


def _significant_positive(row: Mapping[str, str], alpha: float) -> bool:
    return (
        row.get("analysis_status") == "VALID"
        and row.get("direction_pass") == "YES"
        and _parse_float(str(row.get("beta")), "effect") > 0
        and _parse_float(str(row.get("bh_fdr")), "FDR") <= alpha
    )


def classify_replication(
    trait_id: str,
    allen_rows: Sequence[Mapping[str, str]],
    gse_rows: Sequence[Mapping[str, str]],
    alignment: Sequence[Mapping[str, str]],
    *,
    alpha: float = 0.05,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not math.isfinite(alpha) or alpha != 0.05:
        raise E.ContractError("replication alpha differs from the frozen threshold")
    if len(alignment) != 7 or [row.get("gse_class") for row in alignment] != list(E.GSE_ORDER):
        raise E.ContractError("replication alignment is not the complete ordered GSE family")
    allen = {str(row["class_id"]): row for row in allen_rows}
    gse = {str(row["class_id"]): row for row in gse_rows}
    if set(allen) != set(E.CLASS_ORDER) or set(gse) != set(E.GSE_ORDER):
        raise E.ContractError("replication input lacks an exact seven-class atlas family")
    output: list[dict[str, object]] = []
    claim = "Cross-atlas corrected directional support is not causal-cell, subtype, pericyte, or locus-mediation evidence."
    for link in alignment:
        gse_class = str(link["gse_class"])
        replication = gse[gse_class]
        targets = [] if link["allen_broad_classes"] == "NA" else str(link["allen_broad_classes"]).split(";")
        if gse_class == "hybrid":
            classification = (
                "UNALIGNED_HYBRID_DENOMINATOR_RETAINED"
                if replication["analysis_status"] == "VALID"
                else "UNAVAILABLE_UNALIGNED_HYBRID_TEST_FAILED_DENOMINATOR_RETAINED"
            )
            output.append({
                "target_trait_id": trait_id, "allen_class": "NA", "gse_class": gse_class,
                "alignment_status": link["alignment_status"],
                "replication_resolution": link["replication_resolution"],
                "discovery_status": "NOT_ALIGNED", "discovery_beta": "NA",
                "discovery_bh_fdr": "NA", "discovery_direction_pass": "NA",
                "replication_status": replication["analysis_status"],
                "replication_beta": replication["beta"], "replication_bh_fdr": replication["bh_fdr"],
                "replication_direction_pass": replication["direction_pass"],
                "discovery_eligible": "NA", "direction_concordant": "NA",
                "both_atlas_fdr_pass": "NA", "replication_classification": classification,
                "independent_replication_credit": 0, "claim_limit": claim,
            })
            continue
        if not targets:
            raise E.ContractError(f"aligned GSE class lacks an Allen target: {gse_class}")
        for allen_class in targets:
            discovery = allen.get(allen_class)
            if discovery is None:
                raise E.ContractError(f"alignment names an absent Allen class: {allen_class}")
            discovery_eligible = _significant_positive(discovery, alpha)
            replication_valid = replication["analysis_status"] == "VALID"
            if discovery["analysis_status"] != "VALID":
                classification = "UNAVAILABLE_DISCOVERY_TEST_FAILED"
            elif not discovery_eligible:
                classification = "NOT_ELIGIBLE_DISCOVERY_FDR_OR_DIRECTION"
            elif not replication_valid:
                classification = "UNAVAILABLE_REPLICATION_TEST_FAILED"
            else:
                discovery_beta = _parse_float(str(discovery["beta"]), "Allen beta")
                replication_beta = _parse_float(str(replication["beta"]), "GSE beta")
                direction_ok = discovery_beta > 0 and replication_beta > 0
                replication_fdr = _parse_float(str(replication["bh_fdr"]), "GSE FDR")
                if not direction_ok:
                    classification = "FAILED_REPLICATION_DIRECTION"
                elif replication_fdr > alpha:
                    classification = "FAILED_REPLICATION_FDR"
                elif link["alignment_status"] == "ONE_TO_ONE":
                    classification = "REPLICATED_CELL_CLASS"
                elif link["alignment_status"] == "MANY_TO_ONE" and gse_class == "neurons":
                    classification = "PARTIAL_COARSE_NEURON"
                else:
                    raise E.ContractError("alignment resolution is not supported by the frozen classifier")
            direction_concordant = "NA"
            both_fdr = "NA"
            if discovery["analysis_status"] == "VALID" and replication_valid:
                direction_concordant = "YES" if (
                    _parse_float(str(discovery["beta"]), "Allen beta") > 0
                    and _parse_float(str(replication["beta"]), "GSE beta") > 0
                ) else "NO"
                both_fdr = "YES" if (
                    _parse_float(str(discovery["bh_fdr"]), "Allen FDR") <= alpha
                    and _parse_float(str(replication["bh_fdr"]), "GSE FDR") <= alpha
                ) else "NO"
            credit = 1 if classification == "REPLICATED_CELL_CLASS" else 0
            output.append({
                "target_trait_id": trait_id, "allen_class": allen_class, "gse_class": gse_class,
                "alignment_status": link["alignment_status"],
                "replication_resolution": link["replication_resolution"],
                "discovery_status": discovery["analysis_status"], "discovery_beta": discovery["beta"],
                "discovery_bh_fdr": discovery["bh_fdr"],
                "discovery_direction_pass": discovery["direction_pass"],
                "replication_status": replication["analysis_status"],
                "replication_beta": replication["beta"], "replication_bh_fdr": replication["bh_fdr"],
                "replication_direction_pass": replication["direction_pass"],
                "discovery_eligible": "YES" if discovery_eligible else "NO",
                "direction_concordant": direction_concordant, "both_atlas_fdr_pass": both_fdr,
                "replication_classification": classification,
                "independent_replication_credit": credit, "claim_limit": claim,
            })
    if len(output) != 8 or sum(row["gse_class"] == "hybrid" for row in output) != 1:
        raise E.ContractError("replication output is not the exact seven-target/eight-alignment family")
    partial_rows = [row for row in output if row["replication_classification"] == "PARTIAL_COARSE_NEURON"]
    summary = {
        "schema_version": "track-b-p16-cell-class-replication-summary.1",
        "target_trait_id": trait_id,
        "classification_row_count": 8,
        "allen_bh_denominator": 7,
        "gse_bh_denominator": 7,
        "hybrid_tested_and_denominator_retained": True,
        "hybrid_aligned": False,
        "independent_one_to_one_replication_count": sum(int(row["independent_replication_credit"]) for row in output),
        "partial_coarse_neuron_allen_row_count": len(partial_rows),
        "partial_coarse_neuron_shared_gse_target_count": 1 if partial_rows else 0,
        "many_to_one_double_credit_allowed": False,
        "classification_counts": {
            name: sum(row["replication_classification"] == name for row in output)
            for name in sorted({str(row["replication_classification"]) for row in output})
        },
        "pericyte_claim_allowed": False,
        "mouse_substitution_allowed": False,
        "claim_limit": claim,
    }
    return output, summary


def run_classifier(root: Path, enrichment_relative: str, out_relative: str) -> dict[str, object]:
    policy = E.validate_runtime_contract(root)
    enrichment_relative = E.safe_relative(enrichment_relative)
    provenance, allen_rows, gse_rows, input_records = validate_enrichment_family(
        root, enrichment_relative, policy,
    )
    pair_id, trait_id = str(provenance["pair_id"]), str(provenance["target_trait_id"])
    production = policy["production_outputs"]
    assert isinstance(production, dict)
    expected_out = str(production["replication_dir_template"]).format(pair_id=pair_id, trait_id=trait_id)
    if E.safe_relative(out_relative) != expected_out:
        raise E.ContractError(f"replication output must use the frozen pair/trait namespace: {expected_out}")
    try:
        os.lstat(root / expected_out)
    except FileNotFoundError:
        pass
    else:
        raise E.ContractError("replication output directory already exists; refusing replacement")
    alignment = E.parse_tsv_payload(
        E.read_stable_bytes(root, E.ALIGNMENT_REL), E.ALIGNMENT_FIELDS, "Allen/GSE alignment",
    )
    if alignment != E.parse_tsv_payload(
        E.tsv_bytes(E.ALIGNMENT_FIELDS, E.alignment_rows()), E.ALIGNMENT_FIELDS, "expected alignment",
    ):
        raise E.ContractError("Allen/GSE alignment semantics differ from the frozen classifier")
    rows, summary = classify_replication(trait_id, allen_rows, gse_rows, alignment)
    table_payload = E.tsv_bytes(REPLICATION_FIELDS, rows)
    summary_payload = E.canonical_json(summary)
    output_payloads: dict[str, bytes] = {
        "P16_CLASS_REPLICATION.tsv": table_payload,
        "P16_CLASS_REPLICATION_SUMMARY.json": summary_payload,
    }
    classifier_provenance = {
        "schema_version": "track-b-p16-cell-class-replication-run.1",
        "terminal_status": "COMPLETE_WITH_FAILED_AND_UNALIGNED_ROWS_PRESERVED",
        "pair_id": pair_id, "target_trait_id": trait_id,
        "contract_lock_path": E.CONTRACT_LOCK_REL,
        "contract_lock_sha256": E.stable_file(root, E.CONTRACT_LOCK_REL).sha256,
        "enrichment_directory": enrichment_relative,
        "enrichment_provenance_sha256": input_records["P14_P16_ENRICHMENT.provenance.json"].sha256,
        "enrichment_input_sha256": {
            name: record.sha256 for name, record in sorted(input_records.items())
        },
        "alignment_path": E.ALIGNMENT_REL,
        "alignment_sha256": E.stable_file(root, E.ALIGNMENT_REL).sha256,
        "tested_statistic": "MAGMA_V1.10_GENE_PROPERTY_BETA_OVER_SE",
        "thresholds": {
            "allen_bh_denominator": 7, "gse_bh_denominator": 7,
            "fdr_alpha": 0.05, "direction": "BETA_GT_ZERO_IN_BOTH_ATLASES",
        },
        "hybrid_rule": "TESTED_AND_GSE_BH_DENOMINATOR_RETAINED_BUT_UNALIGNED",
        "neurons_rule": "EXC_AND_INH_SHARE_ONE_GSE_NEURONS_TARGET_PARTIAL_ONLY_NO_DOUBLE_CREDIT",
        "output_sha256": {
            name: hashlib.sha256(payload).hexdigest() for name, payload in sorted(output_payloads.items())
        },
        "claim_limit": summary["claim_limit"],
    }
    output_payloads["P16_CLASS_REPLICATION.provenance.json"] = E.canonical_json(classifier_provenance)
    if set(output_payloads) != set(production["replication_required_files"]):
        raise E.ContractError("replication output family differs from the frozen production family")
    E.publish_directory_no_replace(
        root, out_relative, output_payloads,
        completion_name="P16_CLASS_REPLICATION.provenance.json",
    )
    return classifier_provenance


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT_AT_IMPORT)
    parser.add_argument("--enrichment-dir", required=True, help="exact no-replace enrichment output directory")
    parser.add_argument("--out-dir", required=True, help="exact no-replace replication output directory")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        raise E.ContractError("repository root is not a directory")
    result = run_classifier(root, args.enrichment_dir, args.out_dir)
    print(f"P16_CLASS_REPLICATION_OK pair={result['pair_id']} trait={result['target_trait_id']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except E.ContractError as exc:
        raise SystemExit(f"ERROR: {exc}")
