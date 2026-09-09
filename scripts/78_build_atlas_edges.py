#!/usr/bin/env python3
"""Build traceable atlas graph edges and the locked major-conclusion family."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import downstream_contract

EDGE_FIELDS = [
    "edge_id", "source", "target", "relationship", "effect", "direction",
    "p_value", "fdr", "tissue", "cell_type", "method", "evidence_level",
    "dataset", "version", "provenance_id",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def edge_id(prefix: str, *values: str) -> str:
    digest = hashlib.sha256("\x1f".join(values).encode()).hexdigest()[:20]
    return f"EDGE__{prefix}__{digest}"


def edge(
    identity: str, source: str, target: str, relationship: str, *, effect: str = "NA",
    direction: str = "NA", p_value: str = "NA", fdr: str = "NA", tissue: str = "NA",
    cell_type: str = "NA", method: str, evidence_level: str, dataset: str,
    version: str, provenance_id: str,
) -> dict[str, str]:
    return {
        "edge_id": identity, "source": source, "target": target,
        "relationship": relationship, "effect": effect, "direction": direction,
        "p_value": p_value, "fdr": fdr, "tissue": tissue, "cell_type": cell_type,
        "method": method, "evidence_level": evidence_level, "dataset": dataset,
        "version": version, "provenance_id": provenance_id,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--edges-out", default="results/atlas/edges.tsv")
    parser.add_argument("--conclusions-out", default="results/tables/major_conclusions.tsv")
    parser.add_argument("--provenance-out", default="results/atlas/edges.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    paths = {
        name: root / f"results/atlas/{name}.tsv"
        for name in ("traits", "trait_pairs", "loci", "variants", "genes", "regulatory_elements", "cell_types", "pathways", "causal_tests")
    }
    rows = {name: read_tsv(path) for name, path in paths.items()}
    edges: list[dict[str, str]] = []
    conclusions_seed: list[tuple[str, str, str, str, str, Path, str]] = []

    for pair in rows["trait_pairs"]:
        if pair["analysis_tier"] == "PRIMARY_PHASE1" and pair["global_rg_primary_significant"] == "TRUE":
            identity = edge_id("RG", pair["pair_id"])
            edges.append(edge(
                identity, pair["sleep_trait"], pair["non_sleep_trait"], "genetically_correlated_with",
                effect=pair["global_rg"], direction=pair["effect_direction"], p_value=pair["global_rg_p"],
                fdr=pair["global_rg_fdr_primary_372"], method="LDSC_rg",
                evidence_level="ROBUST_GLOBAL_ASSOCIATION", dataset="locked_45_trait_GWAS",
                version="atlas-v1.0", provenance_id="PAIR:" + sha256(paths["trait_pairs"]),
            ))
            conclusions_seed.append((
                "RG__" + pair["pair_id"],
                f"Primary-phase genetic correlation between {pair['sleep_trait']} and {pair['non_sleep_trait']}",
                "ROBUST_GLOBAL_ASSOCIATION", pair["pair_id"],
                f"rg={pair['global_rg']};p={pair['global_rg_p']};fdr={pair['global_rg_fdr_primary_372']}",
                paths["trait_pairs"], identity,
            ))

    for locus in rows["loci"]:
        pair_id = f"{locus['sleep_trait']}__{locus['non_sleep_trait']}"
        for trait in (locus["sleep_trait"], locus["non_sleep_trait"]):
            edges.append(edge(
                edge_id("LOCUS", trait, locus["locus_id"]), trait, locus["locus_id"], "associated_locus",
                direction=locus["effect_direction"], p_value=locus["placo_lead_p"], fdr=locus["conjfdr_min_q"],
                method=locus["method_support"], evidence_level="CROSS_METHOD_SHARED_LOCUS",
                dataset="PLACO_plus_and_conjFDR", version="atlas-v1.0", provenance_id=locus["provenance_id"],
            ))
        conclusions_seed.append((
            "LOCUS__" + locus["locus_id"],
            f"Cross-method shared locus {locus['locus_id']} for {pair_id}",
            "CROSS_METHOD_SHARED_LOCUS", locus["locus_id"],
            f"placo_p={locus['placo_lead_p']};conjfdr={locus['conjfdr_min_q']};direction={locus['effect_direction']}",
            paths["loci"], edge_id("LOCUS", locus["sleep_trait"], locus["locus_id"]),
        ))

    for variant in rows["variants"]:
        edges.append(edge(
            edge_id("VAR", variant["locus_id"], variant["variant_id"]),
            variant["locus_id"], variant["variant_id"], "contains_variant",
            effect=variant["shared_signal_posterior"], method=variant["fine_mapping_method"],
            evidence_level="FINE_MAPPED_VARIANT", dataset=variant["ld_reference"], version="atlas-v1.0",
            provenance_id=variant["provenance_id"],
        ))

    for gene in rows["genes"]:
        identity = edge_id("GENE", gene["locus_id"], gene["gene_id"])
        edges.append(edge(
            identity, gene["locus_id"], gene["gene_id"], "supports_gene",
            tissue=gene["context"], method="molecular_convergence", evidence_level=gene["evidence_level"],
            dataset="molecular_QTL_and_TWAS", version="atlas-v1.0", provenance_id=gene["provenance_id"],
        ))
        if gene["evidence_level"] == "HIGH_CONFIDENCE_CONVERGENT_GENE":
            conclusions_seed.append((
                f"GENE__{gene['locus_id']}__{gene['gene_id']}",
                f"Convergent molecular support for {gene['gene_symbol']} at {gene['locus_id']}",
                gene["evidence_level"], f"{gene['locus_id']}:{gene['gene_id']}",
                f"supported_streams={gene['evidence_stream_count']};contexts={gene['context']}",
                paths["genes"], identity,
            ))

    for regulatory in rows["regulatory_elements"]:
        first = edge_id("REGOVERLAP", regulatory["regulatory_evidence_id"])
        common = {
            "effect": regulatory["effect"], "p_value": regulatory["p_value"],
            "tissue": regulatory["tissue"], "cell_type": regulatory["cell_type"],
            "method": regulatory["link_method"], "evidence_level": regulatory["evidence_level"],
            "dataset": regulatory["source_dataset"], "version": regulatory["source_version"],
            "provenance_id": regulatory["provenance_id"],
        }
        edges.append(edge(first, regulatory["variant_id"], regulatory["regulatory_element_id"], "overlaps_regulatory_element", **common))
        if regulatory["target_gene_id"] not in {"", "NA"}:
            second = edge_id("REGLINK", regulatory["regulatory_evidence_id"], regulatory["target_gene_id"])
            edges.append(edge(second, regulatory["regulatory_element_id"], regulatory["target_gene_id"], "links_to_gene", **common))

    for cell in rows["cell_types"]:
        identity = edge_id("CELL", cell["trait_or_locus_id"], cell["cell_type_id"])
        edges.append(edge(
            identity, cell["trait_or_locus_id"], cell["cell_type_id"], "enriched_in_cell_type",
            effect=cell["effect"], p_value=cell["p_value"], fdr=cell["fdr"], tissue=cell["tissue"],
            cell_type=cell["cell_type"], method=cell["method"], evidence_level=cell["evidence_level"],
            dataset=cell["dataset"], version=cell["version"], provenance_id=cell["provenance_id"],
        ))
        if cell["evidence_level"] == "MULTI_STRATEGY_CELL_TYPE":
            conclusions_seed.append((
                "CELL__" + cell["cell_type_id"],
                f"Multi-strategy cell-type support for {cell['cell_type']} in {cell['trait_or_locus_id']}",
                cell["evidence_level"], cell["cell_type_id"],
                f"effect={cell['effect']};p={cell['p_value']};fdr={cell['fdr']}",
                paths["cell_types"], identity,
            ))

    for pathway in rows["pathways"]:
        identity = edge_id("PATH", pathway["trait_or_locus_id"], pathway["pathway_id"], pathway["resource"])
        edges.append(edge(
            identity, pathway["trait_or_locus_id"], pathway["pathway_id"], "enriched_for_pathway",
            effect=pathway["effect"], p_value=pathway["p_value"], fdr=pathway["fdr"],
            method="competitive_gene_set", evidence_level=pathway["evidence_level"],
            dataset=pathway["resource"], version=pathway["dataset_version"], provenance_id=pathway["provenance_id"],
        ))
        conclusions_seed.append((
            f"PATH__{pathway['resource']}__{pathway['pathway_id']}__{pathway['trait_or_locus_id']}",
            f"FDR-supported {pathway['resource']} pathway {pathway['pathway_name']} for {pathway['trait_or_locus_id']}",
            pathway["evidence_level"], f"{pathway['pathway_id']}:{pathway['trait_or_locus_id']}",
            f"effect={pathway['effect']};p={pathway['p_value']};fdr={pathway['fdr']}",
            paths["pathways"], identity,
        ))

    for causal in rows["causal_tests"]:
        if causal["evidence_level"] != "ROBUST_CAUSAL_EVIDENCE":
            continue
        identity = edge_id("MR", causal["causal_test_id"])
        edges.append(edge(
            identity, causal["exposure"], causal["outcome"], "causal_effect_on",
            effect=causal["effect"], direction=causal["direction"], p_value=causal["p_value"],
            method=causal["method"], evidence_level=causal["evidence_level"],
            dataset="locked_bidirectional_MR", version="atlas-v1.0", provenance_id=causal["provenance_id"],
        ))
        conclusions_seed.append((
            "MR__" + causal["causal_test_id"],
            f"Robust causal-inference evidence from {causal['exposure']} to {causal['outcome']}",
            causal["evidence_level"], causal["causal_test_id"],
            f"method={causal['method']};effect={causal['effect']};p={causal['p_value']}",
            paths["causal_tests"], identity,
        ))

    relationship_rank = {value: index for index, value in enumerate(policy["edge_relationship_order"])}
    if any(row["relationship"] not in relationship_rank for row in edges):
        fail("an atlas edge uses a relationship outside the locked family")
    edges.sort(key=lambda row: (relationship_rank[row["relationship"]], row["source"], row["target"], row["edge_id"]))
    edge_ids = [row["edge_id"] for row in edges]
    if len(edge_ids) != len(set(edge_ids)):
        fail("atlas edge IDs are duplicated")
    edges_text = table_text(EDGE_FIELDS, edges)
    edges_hash = hashlib.sha256(edges_text.encode()).hexdigest()
    conclusions: list[dict[str, str]] = []
    seen_conclusions: set[str] = set()
    for conclusion_id, conclusion, conclusion_type, entity_id, primary, evidence_path, identity in sorted(conclusions_seed):
        if conclusion_id in seen_conclusions:
            continue
        seen_conclusions.add(conclusion_id)
        conclusions.append({
            "conclusion_id": conclusion_id, "conclusion": conclusion,
            "conclusion_type": conclusion_type, "entity_id": entity_id,
            "primary_result": primary, "evidence_path": str(evidence_path.relative_to(root)),
            "evidence_sha256": sha256(evidence_path), "provenance_id": "EDGE:" + edges_hash + ":" + identity,
        })
    conclusion_text = table_text(policy["major_conclusion_fields"], conclusions)
    payloads = {root / args.edges_out: edges_text, root / args.conclusions_out: conclusion_text}
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path),
        "inputs": {str(path.relative_to(root)): sha256(path) for path in paths.values()},
        "edge_rows": len(edges), "major_conclusion_rows": len(conclusions),
        "outputs": {str(path.relative_to(root)): hashlib.sha256(text.encode()).hexdigest() for path, text in payloads.items()},
        "script_sha256": downstream_contract.script_hashes(root, "interpretation"),
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        for path, text in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                fail(f"atlas graph differs from deterministic recomputation: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            fail("atlas graph provenance drifted")
    else:
        if any(path.exists() for path in [*payloads, provenance_path]):
            fail("atlas graph outputs already exist; refusing overwrite")
        for path, text in payloads.items():
            atomic_text(path, text)
        atomic_text(provenance_path, provenance_text)
    if not args.quiet:
        print(f"ATLAS_GRAPH_OK edges={len(edges)} conclusions={len(conclusions)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
