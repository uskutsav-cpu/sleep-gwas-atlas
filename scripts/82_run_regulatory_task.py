#!/usr/bin/env python3
"""Run a checksum-locked GENCODE or SCREEN regulatory-overlap task."""
from __future__ import annotations

import argparse
from bisect import bisect_right
import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

from liftover_chain import load_chain


MISSING = {"", "NA"}
SCREEN_METADATA = re.compile(
    r'\{"name":"([^"]*)","collection":"Core","ontology":"([^"]*)",'
    r'"lifeStage":"([^"]*)","sampleType":"([^"]*)","displayName":"([^"]*)",'
    r'"assays":\[\{"id":"dnase-[^"]*","assay":"dnase","url":"[^"]*",'
    r'"experimentAccession":"(ENCSR[^"]*)"'
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


def populated(value: str) -> bool:
    return value not in MISSING


def versionless_gene(value: str) -> str:
    return value.split(".", 1)[0]


def attributes(value: str) -> dict[str, str]:
    return {key: text for key, text in re.findall(r'(\S+) "([^"]*)";', value)}


def evidence_id(*values: str) -> str:
    digest = hashlib.sha256("\x1f".join(values).encode()).hexdigest()[:24]
    return "AUTOREG::" + digest


def eligible_variants(rows: list[dict[str, str]], locus_id: str) -> list[dict[str, str]]:
    selected = []
    for row in rows:
        if row["locus_id"] != locus_id or row["qc_status"] != "PASS":
            continue
        if not any(populated(row[field]) for field in (
            "credible_set_sleep", "credible_set_non_sleep", "shared_signal_posterior",
        )):
            continue
        selected.append(row)
    return selected


def lifted_variants(
    variants: list[dict[str, str]], chain,
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    mapped: list[dict[str, object]] = []
    exclusions: list[dict[str, str]] = []
    for row in variants:
        status, target = chain.map_point(row["chromosome"], row["position_bp"])
        if status != "mapped" or target is None:
            exclusions.append({"variant_id": row["variant_id"], "reason": status})
            continue
        mapped.append({
            "variant": row, "chromosome_grch38": int(target[0]),
            "position_grch38": int(target[1]), "target_strand": str(target[2]),
        })
    return mapped, exclusions


def promoter_rows(
    task: dict[str, str], policy: dict[str, object], mapped: list[dict[str, object]],
    genes: list[dict[str, str]], gtf_path: Path,
) -> list[dict[str, str]]:
    supported = {
        versionless_gene(row["gene_id"]): row
        for row in genes if row["locus_id"] == task["locus_id"]
    }
    if not supported or not mapped:
        return []
    by_chromosome: dict[int, list[dict[str, object]]] = {}
    for item in mapped:
        by_chromosome.setdefault(int(item["chromosome_grch38"]), []).append(item)
    window = policy["promoter_mapping"]
    rows: list[dict[str, str]] = []
    seen_genes: set[str] = set()
    with gtf_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                fail(f"invalid GENCODE GTF row at line {line_number}")
            if fields[2] != "gene":
                continue
            annotation = attributes(fields[8])
            gene_id = versionless_gene(annotation.get("gene_id", ""))
            if gene_id not in supported:
                continue
            try:
                chromosome = int(fields[0].removeprefix("chr"))
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            if chromosome not in by_chromosome or fields[6] not in {"+", "-"}:
                continue
            if gene_id in seen_genes:
                fail(f"duplicate GENCODE gene record for supported gene {gene_id}")
            seen_genes.add(gene_id)
            tss = start if fields[6] == "+" else end
            if fields[6] == "+":
                promoter_start = max(1, tss - int(window["upstream_bp"]))
                promoter_end = tss + int(window["downstream_bp"])
            else:
                promoter_start = max(1, tss - int(window["downstream_bp"]))
                promoter_end = tss + int(window["upstream_bp"])
            element = f"GENCODEv26::{gene_id}::chr{chromosome}:{promoter_start}-{promoter_end}:{fields[6]}"
            for item in by_chromosome[chromosome]:
                position = int(item["position_grch38"])
                if not promoter_start <= position <= promoter_end:
                    continue
                variant = item["variant"]
                rows.append({
                    "regulatory_evidence_id": evidence_id(task["task_id"], str(variant["variant_id"]), element, gene_id),
                    "regulatory_element_id": element, "variant_id": str(variant["variant_id"]),
                    "locus_id": task["locus_id"], "element_type": "promoter",
                    "annotation": f"TSS_WINDOW_UPSTREAM_{window['upstream_bp']}_DOWNSTREAM_{window['downstream_bp']}",
                    "effect": "NA", "p_value": "NA",
                    "biosample": "CONTEXT_INDEPENDENT_GENCODE_ANNOTATION", "tissue": "NA",
                    "cell_type": "NA", "context_domain": task["domain"],
                    "target_gene_id": gene_id, "link_method": "GENCODE_TSS_PROMOTER_OVERLAP",
                    "source_dataset": task["source_id"], "source_version": task["source_release"],
                    "evidence_level": "GENOMIC_PROMOTER_ANNOTATION", "provenance_id": "ADAPTER_PENDING",
                })
    return rows


def parse_screen_metadata(path: Path) -> dict[str, dict[str, str]]:
    rows = SCREEN_METADATA.findall(path.read_text(encoding="utf-8"))
    metadata = {
        accession: {
            "name": name.replace("\\'", "'"), "ontology": ontology,
            "life_stage": life_stage, "sample_type": sample_type,
            "display_name": display_name.replace("\\'", "'"),
        }
        for name, ontology, life_stage, sample_type, display_name, accession in rows
    }
    if len(rows) != len(metadata):
        fail("SCREEN Core Collection metadata contains duplicate DNase accessions")
    return metadata


def screen_components(root: Path, policy: dict[str, object]) -> tuple[dict[str, Path], dict[str, object]]:
    spec = policy["screen_registry_v4"]
    manifest_path = root / spec["component_manifest"]
    if sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("SCREEN component manifest differs from the policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components: dict[str, Path] = {}
    for component in manifest["components"]:
        path = root / component["path"]
        if (
            not path.is_file() or path.stat().st_size != component["bytes"]
            or sha256(path) != component["sha256"]
        ):
            fail(f"SCREEN component differs from its exact pin: {component['component_id']}")
        components[component["component_id"]] = path
    return components, manifest


def screen_rows(
    task: dict[str, str], policy: dict[str, object], mapped: list[dict[str, object]],
    coordinate_path: Path, matrix_path: Path, metadata_path: Path,
) -> list[dict[str, str]]:
    if not mapped:
        return []
    metadata = parse_screen_metadata(metadata_path)
    domain_ontologies = set(policy["screen_registry_v4"]["domain_ontologies"][task["domain"]])
    by_chromosome: dict[int, list[tuple[int, dict[str, object]]]] = {}
    for item in mapped:
        by_chromosome.setdefault(int(item["chromosome_grch38"]), []).append((int(item["position_grch38"]), item))
    for values in by_chromosome.values():
        values.sort(key=lambda value: (value[0], str(value[1]["variant"]["variant_id"])))
    overlaps: dict[str, list[tuple[dict[str, object], str]]] = {}
    with coordinate_path.open(encoding="ascii") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                fail(f"invalid SCREEN coordinate row at line {line_number}")
            try:
                chromosome = int(fields[0].removeprefix("chr"))
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            candidates = by_chromosome.get(chromosome, [])
            if not candidates:
                continue
            positions = [value[0] for value in candidates]
            left, right = bisect_right(positions, start), bisect_right(positions, end)
            if left == right:
                continue
            ccre = fields[4]
            overlaps.setdefault(ccre, []).extend((item, fields[5]) for _, item in candidates[left:right])
    if not overlaps:
        return []
    layer = "enhancer" if task["source_id"] == "ENCODE_CCRE_V4_ENHANCER" else "open_chromatin"
    allowed = set(
        policy["screen_registry_v4"][
            "enhancer_classes" if layer == "enhancer" else "open_chromatin_classes"
        ]
    )
    rows: list[dict[str, str]] = []
    with gzip.open(matrix_path, "rt", encoding="ascii", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        selected = [
            (index, accession) for index, accession in enumerate(header[1:], start=1)
            if accession in metadata and metadata[accession]["ontology"] in domain_ontologies
        ]
        if not selected:
            fail(f"SCREEN Core Collection has no frozen biosample for domain {task['domain']}")
        for fields in reader:
            ccre = fields[0]
            if ccre not in overlaps:
                continue
            for index, accession in selected:
                cell_class = fields[index]
                if cell_class not in allowed:
                    continue
                sample = metadata[accession]
                cell_type = sample["display_name"] if sample["sample_type"] != "tissue" else "NA"
                for item, global_class in overlaps[ccre]:
                    variant = item["variant"]
                    rows.append({
                        "regulatory_evidence_id": evidence_id(
                            task["task_id"], str(variant["variant_id"]), ccre, accession, cell_class,
                        ),
                        "regulatory_element_id": ccre, "variant_id": str(variant["variant_id"]),
                        "locus_id": task["locus_id"], "element_type": layer,
                        "annotation": f"CELL_CLASS::{cell_class}::GLOBAL_CLASS::{global_class}",
                        "effect": "NA", "p_value": "NA",
                        "biosample": accession, "tissue": sample["ontology"], "cell_type": cell_type,
                        "context_domain": task["domain"], "target_gene_id": "NA",
                        "link_method": f"SCREEN_CORE_COLLECTION_CLASS::{cell_class}::GLOBAL::{global_class}",
                        "source_dataset": task["source_id"], "source_version": task["source_release"],
                        "evidence_level": "CONTEXT_SPECIFIC_CCRE_OVERLAP", "provenance_id": "ADAPTER_PENDING",
                    })
    return rows


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
    manifest_path = root / args.manifest
    _, tasks = read_tsv(manifest_path)
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    automatic_sources = {policy["promoter_mapping"]["source_id"], *policy["screen_registry_v4"]["source_ids"]}
    if task["analysis_family"] != "regulatory" or task["method"] != "interval_overlap" or task["source_id"] not in automatic_sources:
        fail("task is not supported by the automatic regulatory adapter")
    _, sources = read_tsv(root / policy["source_registry"])
    source = next((row for row in sources if row["source_id"] == task["source_id"]), None)
    if source is None or source["source_status"] != "SOURCE_VERIFIED" or source["exact_release"] != task["source_release"]:
        fail("regulatory source is absent, unverified, or release-mismatched")
    source_path = root / source["local_path"]
    if (
        not source_path.is_file() or source_path.stat().st_size != int(source["expected_bytes"])
        or sha256(source_path) != source["expected_sha256"]
    ):
        fail("regulatory source differs from its exact byte/SHA-256 pin")
    _, variants = read_tsv(root / "results/atlas/variants.tsv")
    _, genes = read_tsv(root / "results/atlas/genes.tsv")
    eligible = eligible_variants(variants, task["locus_id"])
    chain_spec = policy["regulatory_build_harmonization"]
    chain, chain_provenance = load_chain(
        root / chain_spec["chain_path"], chain_spec["chain_sha256"], chain_spec["chain_bytes"],
    )
    mapped, exclusions = lifted_variants(eligible, chain)
    component_hashes: dict[str, str] = {}
    if task["source_id"] == policy["promoter_mapping"]["source_id"]:
        rows = promoter_rows(task, policy, mapped, genes, source_path)
    else:
        components, _ = screen_components(root, policy)
        component_hashes = {identity: sha256(path) for identity, path in sorted(components.items())}
        rows = screen_rows(
            task, policy, mapped, components["GRCH38_CCRE_COORDINATES"],
            components["CORE_COLLECTION_CLASS_MATRIX"],
            components["DOWNLOADS_BIOSAMPLE_METADATA_SNAPSHOT"],
        )
    rows.sort(key=lambda row: (
        row["variant_id"], row["regulatory_element_id"], row["biosample"],
        row["target_gene_id"], row["regulatory_evidence_id"],
    ))
    fields = policy["regulatory_mapping"]["canonical_fields"]
    if any(list(row) != fields for row in rows):
        fail("automatic regulatory row differs from the locked schema")
    payload = table_text(fields, rows)
    out_path = root / args.out
    atomic_text(out_path, payload)
    status = "COMPLETED" if rows else "NO_EVIDENCE_FOUND"
    reason = (
        f"Checksum-locked {task['source_id']} overlap produced {len(rows)} regulatory evidence rows"
        if rows else
        f"No eligible uniquely lifted fine-mapped variant produced a source-qualified {source['layer_or_resource']} overlap in {task['domain']}"
    )
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"], "terminal_status": status, "terminal_reason": reason,
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"], "source_id": task["source_id"],
        "source_release": task["source_release"], "source_sha256": sha256(source_path),
        "screen_component_sha256": component_hashes, "chain": chain_provenance,
        "eligible_variant_count": len(eligible), "mapped_variant_count": len(mapped),
        "liftover_exclusions": exclusions, "result_rows": len(rows),
        "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "claim_limit": policy["regulatory_mapping"]["claim_limit"],
    }
    provenance_path = root / args.provenance_out
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"REGULATORY_ADAPTER_OK task={task['task_id']} status={status} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
