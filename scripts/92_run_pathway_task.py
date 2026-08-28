#!/usr/bin/env python3
"""Run one locked public Reactome or GO pathway task without post-result choices."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

from pathway_sources import (
    hypergeometric_right_tail,
    load_go_sets,
    load_magma_gene_sets,
    load_msigdb_symbol_sets,
    load_reactome_sets,
    load_unique_gencode_symbols,
)


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


def table_text(fields: list[str], rows: list[dict[str, str]]) -> str:
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


def check_file(path: Path, expected: dict[str, object], label: str) -> None:
    if (
        not path.is_file() or path.stat().st_size != expected["bytes"]
        or sha256(path) != expected["sha256"]
    ):
        fail(f"{label} differs from its exact source pin")


def source_sets(
    root: Path, bundle: dict[str, object], source_id: str, symbol_map: dict[str, str],
    minimum_size: int, maximum_size: int,
) -> tuple[dict[str, tuple[str, set[str]]], dict[str, object], list[Path]]:
    resource = bundle["resources"][source_id]
    if source_id == "REACTOME":
        archive = root / resource["path"]
        check_file(archive, resource, "Reactome v97 archive")
        sets, observed = load_reactome_sets(archive, symbol_map, minimum_size, maximum_size)
        files = [archive]
        expected = {
            key: resource[key] for key in (
                "archive_member", "archive_member_bytes", "archive_member_sha256",
                "source_pathways", "source_symbols", "mapped_genes", "eligible_sets",
                "eligible_gene_memberships",
            )
        }
    elif source_id == "GO":
        ontology = root / resource["ontology"]["path"]
        annotation = root / resource["annotation"]["path"]
        check_file(ontology, resource["ontology"], "GO ontology")
        check_file(annotation, resource["annotation"], "GO human annotation")
        sets, observed = load_go_sets(
            ontology, annotation, symbol_map, minimum_size, maximum_size,
        )
        files = [ontology, annotation]
        expected = {
            "ontology_data_version": resource["ontology"]["data_version"],
            "active_ontology_terms": resource["ontology"]["active_terms"],
            "gaf_rows": resource["annotation"]["rows"],
            "gaf_date_generated": resource["annotation"]["date_generated"],
            "gaf_go_version": resource["annotation"]["go_version"],
            **{
                key: resource[key] for key in (
                    "excluded_not_rows", "excluded_unmapped_rows", "direct_annotated_terms",
                    "propagated_terms", "mapped_genes", "eligible_sets",
                    "eligible_gene_memberships",
                )
            },
        }
    elif source_id == "MSIGDB":
        source = root / resource["path"]
        check_file(source, resource, "MSigDB v2026.1 human symbols GMT")
        sets, observed = load_msigdb_symbol_sets(
            source, symbol_map, minimum_size, maximum_size,
        )
        files = [source]
        expected = {
            key: resource[key] for key in (
                "source_sets", "source_symbols", "mapped_genes", "eligible_sets",
                "eligible_gene_memberships",
            )
        }
    elif source_id == "MAGMA_GENE_SETS":
        source = root / resource["path"]
        check_file(source, resource, "FUMA MSigDB v2023.1Hs MAGMA GMT")
        sets, observed = load_magma_gene_sets(
            source, set(symbol_map.values()), minimum_size, maximum_size,
        )
        files = [source]
        expected = {
            key: resource[key] for key in (
                "source_sets", "source_gene_ids", "mapped_genes", "eligible_sets",
                "eligible_gene_memberships",
            )
        }
    else:
        fail(f"unsupported automatic public pathway source: {source_id}")
    observed["eligible_gene_union"] = len(set().union(*(genes for _, genes in sets.values())))
    expected["eligible_gene_union"] = resource["eligible_gene_union"]
    if observed != expected:
        fail(f"{source_id} pathway family differs from its pre-result manifest counts")
    return sets, observed, files


def fmt(value: float) -> str:
    return format(value, ".12g")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--provenance-out", required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["public_pathway_sources"]
    task_manifest_path = root / args.manifest
    _, tasks = read_tsv(task_manifest_path)
    selected_tasks = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected_tasks) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected_tasks[0]
    if (
        task["analysis_family"] != "pathway" or task["method"] != "competitive_gene_set"
        or task["source_id"] not in spec["source_ids"]
    ):
        fail("task is not supported by the public pathway adapter")
    bundle_path = root / spec["component_manifest"]
    if not bundle_path.is_file() or sha256(bundle_path) != spec["component_manifest_sha256"]:
        fail("public pathway source manifest differs from its policy pin")
    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    resource = bundle["resources"][task["source_id"]]
    if task["source_release"] != resource["release"]:
        fail("pathway task release differs from its source manifest")
    gencode = root / spec["gencode_path"]
    check_file(gencode, bundle["identifier_mapping"], "GENCODE identifier reference")
    try:
        symbol_map, gencode_observed = load_unique_gencode_symbols(gencode)
    except (OSError, UnicodeError, ValueError) as exc:
        fail(f"GENCODE identifier reference is invalid: {exc}")
    expected_gencode = {
        key: bundle["identifier_mapping"][key] for key in (
            "gene_rows", "distinct_symbols", "unique_symbols", "ambiguous_symbols",
        )
    }
    if gencode_observed != expected_gencode:
        fail("GENCODE exact-symbol mapping family differs from its pre-result manifest")
    minimum_size, maximum_size = map(int, bundle["set_size_range"])
    try:
        sets, source_observed, source_files = source_sets(
            root, bundle, task["source_id"], symbol_map, minimum_size, maximum_size,
        )
    except (OSError, UnicodeError, ValueError) as exc:
        fail(f"{task['source_id']} source is invalid: {exc}")
    genes_path, loci_path = root / "results/atlas/genes.tsv", root / "results/atlas/loci.tsv"
    _, gene_rows = read_tsv(genes_path)
    _, locus_rows = read_tsv(loci_path)
    locus_ids = {row["locus_id"] for row in locus_rows}
    high_confidence: dict[str, set[str]] = defaultdict(set)
    excluded_unmapped: dict[str, int] = defaultdict(int)
    for row in gene_rows:
        if row["evidence_level"] != "HIGH_CONFIDENCE_CONVERGENT_GENE":
            continue
        locus_id, gene_id, symbol = row["locus_id"], row["gene_id"], row["gene_symbol"]
        if locus_id not in locus_ids:
            fail(f"high-confidence gene references an absent locus: {locus_id}")
        mapped = symbol_map.get(symbol)
        if mapped is None:
            excluded_unmapped[locus_id] += 1
            continue
        if mapped != gene_id.split(".", 1)[0]:
            fail(f"high-confidence gene ID/symbol differs from exact GENCODE mapping: {locus_id}/{gene_id}")
        high_confidence[locus_id].add(mapped)
    universe_genes = set().union(*(genes for _, genes in sets.values()))
    universe = len(universe_genes)
    if universe != resource["eligible_gene_union"]:
        fail("pathway background universe differs from its source manifest")
    fields = policy["pathways"]["canonical_fields"]
    rows: list[dict[str, str]] = []
    tested_loci: dict[str, int] = {}
    for locus_id in sorted(high_confidence):
        selected = high_confidence[locus_id] & universe_genes
        excluded_unmapped[locus_id] += len(high_confidence[locus_id] - universe_genes)
        if not selected:
            continue
        tested_loci[locus_id] = len(selected)
        for pathway_id, (pathway_name, gene_set) in sorted(sets.items()):
            contributing = sorted(selected & gene_set)
            overlap, set_size, selected_size = len(contributing), len(gene_set), len(selected)
            p_value = hypergeometric_right_tail(overlap, universe, set_size, selected_size)
            effect = overlap * universe / (selected_size * set_size)
            rows.append({
                "pathway_id": pathway_id, "pathway_name": pathway_name,
                "resource": {
                    "REACTOME": "Reactome", "MAGMA_GENE_SETS": "MAGMA_gene_sets",
                }.get(task["source_id"], task["source_id"]),
                "gene_set_size": str(set_size),
                "contributing_genes": ";".join(contributing) if contributing else "NA",
                "trait_or_locus_id": locus_id, "effect": fmt(effect),
                "p_value": fmt(p_value), "fdr": "NA",
                "evidence_level": "UNADJUSTED_ONE_SIDED_HYPERGEOMETRIC",
                "dataset_version": task["source_release"],
                "provenance_id": f"PATHSRC__{task['source_id']}__{task['source_release']}",
            })
    if any(list(row) != fields for row in rows):
        fail("public pathway result differs from the locked canonical schema")
    rows.sort(key=lambda row: (row["trait_or_locus_id"], row["pathway_id"]))
    payload = table_text(fields, rows)
    out_path, provenance_path = root / args.out, root / args.provenance_out
    if out_path.exists() or provenance_path.exists():
        fail("public pathway task output already exists; refusing overwrite")
    atomic_text(out_path, payload)
    status = "COMPLETED" if rows else "NO_EVIDENCE_FOUND"
    reason = (
        f"Tested all {len(sets)} prespecified {task['source_id']} sets in "
        f"{len(tested_loci)} loci with source-mapped high-confidence genes"
        if rows else
        "No high-confidence convergent atlas gene mapped exactly into the locked pathway background"
    )
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"], "terminal_status": status, "terminal_reason": reason,
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(task_manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"], "source_id": task["source_id"],
        "source_release": task["source_release"], "source_manifest_sha256": sha256(bundle_path),
        "source_files": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in source_files
        ],
        "gencode_path": str(gencode), "gencode_sha256": sha256(gencode),
        "genes_path": str(genes_path), "genes_sha256": sha256(genes_path),
        "loci_path": str(loci_path), "loci_sha256": sha256(loci_path),
        "source_validation": source_observed, "background_gene_count": universe,
        "tested_loci_selected_gene_counts": tested_loci,
        "excluded_high_confidence_gene_counts": dict(sorted(excluded_unmapped.items())),
        "test": bundle["enrichment"]["test"], "effect": bundle["enrichment"]["effect"],
        "result_rows": len(rows), "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "claim_limit": bundle["claim_limit"],
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"PATHWAY_TASK_OK task={task['task_id']} status={status} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
