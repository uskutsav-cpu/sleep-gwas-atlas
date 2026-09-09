#!/usr/bin/env python3
"""Build evidence-linked mechanistic graphs and a noncausal ranked synthesis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


CHAIN_FIELDS = [
    "chain_edge_id", "comparison_id", "pair_id", "locus_id", "signal1", "signal2",
    "gene_id", "gene_symbol", "edge_order", "edge_from", "edge_to", "entity_from_ids",
    "entity_to_ids", "edge_status", "evidence_ids", "evidence_provenance", "evidence_classes",
    "source_ids", "contradiction_status", "caveat", "claim_limit",
]
SYNTHESIS_FIELDS = [
    "mechanism_id", "comparison_id", "pair_id", "locus_id", "comparison_type", "signal1", "signal2",
    "gene_id", "gene_symbol", "top_shared_variant", "top_shared_variant_PP_H4",
    "upstream_replication_status", "upstream_shared_signal_status", "molecular_qtl_colocalization_support",
    "credible_variant_regulatory_support", "regulatory_element_target_gene_support",
    "tissue_cell_expression_support", "cell_state_accessibility_support", "pathway_support",
    "model_system_perturbation_support", "supported_evidence_classes", "supported_evidence_ids",
    "source_ids", "contradiction_status", "supported_required_edge_count", "missing_required_edges",
    "strongest_mechanistic_evidence", "confidence_class", "claim_limit",
]
REQUIRED_EDGES = [
    (1, "SLEEP_TRAIT", "SHARED_SIGNAL_MODEL"),
    (2, "SHARED_SIGNAL_MODEL", "CREDIBLE_VARIANT"),
    (3, "CREDIBLE_VARIANT", "REGULATORY_ELEMENT"),
    (4, "REGULATORY_ELEMENT", "TARGET_GENE"),
    (5, "TARGET_GENE", "CELL_TYPE_OR_STATE"),
    (6, "CELL_TYPE_OR_STATE", "PATHWAY"),
    (7, "PATHWAY", "EXTERNAL_PHENOTYPE"),
]
CONFIDENCE_ORDER = {
    "MULTI_LAYER_HUMAN_GENETIC_SUPPORT": 0,
    "COLOCALIZED_MOLECULAR_SUPPORT": 1,
    "ANNOTATION_ONLY_SUPPORT": 2,
    "INSUFFICIENT_OR_CONFLICTING": 3,
}
EVIDENCE_ORDER = [
    "COLOCALIZED_MOLECULAR_QTL", "FINE_MAPPED_REGULATORY_EVIDENCE",
    "TISSUE_CELL_EXPRESSION_OR_ENRICHMENT", "PATHWAY_OR_NETWORK_ANNOTATION",
    "MODEL_SYSTEM_OR_PERTURBATION_EVIDENCE",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def supported(rows: list[dict[str, str]], *, task: str | None = None, gene_id: str | None = None) -> list[dict[str, str]]:
    output = [row for row in rows if row["search_outcome"] == "SUPPORTED"]
    if task is not None:
        output = [row for row in output if row["search_task"] == task]
    if gene_id is not None:
        output = [row for row in output if row["gene_id"] == gene_id]
    return output


def values(rows: list[dict[str, str]], field: str) -> str:
    selected = sorted({row[field] for row in rows if row[field] not in {"", "NA", ".", "PENDING"}})
    return ";".join(selected) if selected else "NONE"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_search_queue.tsv"))
    parser.add_argument("--family-lock", type=Path, default=Path("discovery_extension/config/mechanistic_signal_family.lock.json"))
    parser.add_argument("--evidence", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_evidence.tsv"))
    parser.add_argument("--evidence-provenance", type=Path, default=Path("discovery_extension/provenance/mechanistic_evidence_validation.json"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/mechanistic_annotation_contract.json"))
    parser.add_argument("--fine-mapping", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_colocalization.tsv"))
    parser.add_argument("--fine-mapping-provenance", type=Path, default=Path("discovery_extension/provenance/fine_mapping_colocalization_results.json"))
    parser.add_argument("--chains-out", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_chains.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_synthesis.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/mechanistic_synthesis.json"))
    args = parser.parse_args()

    queue = read_tsv(args.queue)
    lock = json.loads(args.family_lock.read_text(encoding="utf-8"))
    evidence = read_tsv(args.evidence)
    evidence_provenance = json.loads(args.evidence_provenance.read_text(encoding="utf-8"))
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    fine_mapping_provenance = json.loads(args.fine_mapping_provenance.read_text(encoding="utf-8"))
    if sha256(args.queue) != lock.get("queue_template_sha256"):
        raise SystemExit("ERROR: mechanistic queue differs from its pre-search lock")
    if evidence_provenance.get("output_sha256") != sha256(args.evidence):
        raise SystemExit("ERROR: mechanistic evidence differs from its validation provenance")
    if fine_mapping_provenance.get("output_sha256") != sha256(args.fine_mapping) or lock.get("fine_mapping_sha256") != sha256(args.fine_mapping):
        raise SystemExit("ERROR: upstream fine-mapping/colocalization evidence drifted")
    if lock.get("contract_sha256") != sha256(args.contract):
        raise SystemExit("ERROR: mechanistic synthesis contract differs from the search lock")

    queue_by_signal: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in queue:
        key = f"{row['comparison_id']}__{row['signal1']}__{row['signal2']}"
        queue_by_signal[key].append(row)
    evidence_by_signal: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence:
        key = f"{row['comparison_id']}__{row['signal1']}__{row['signal2']}"
        evidence_by_signal[key].append(row)
    if list(queue_by_signal) != lock["signal_keys_in_locked_order"] or set(evidence_by_signal) != set(queue_by_signal):
        raise SystemExit("ERROR: mechanistic synthesis signal family differs from lock/evidence")

    chain_rows: list[dict[str, object]] = []
    synthesis_rows: list[dict[str, object]] = []
    for signal_key in lock["signal_keys_in_locked_order"]:
        tasks = queue_by_signal[signal_key]
        task = tasks[0]
        rows = evidence_by_signal[signal_key]
        supported_rows = supported(rows)
        gene_pairs = sorted({(row["gene_id"], row["gene_symbol"]) for row in supported_rows if row["gene_id"] not in {"", "NA", ".", "PENDING"}})
        if not gene_pairs:
            gene_pairs = [("NO_SUPPORTED_GENE", "NO_SUPPORTED_GENE")]
        for gene_id, gene_symbol in gene_pairs:
            gene_rows = [row for row in supported_rows if row["gene_id"] == gene_id]
            molecular = supported(rows, task="MOLECULAR_QTL", gene_id=gene_id)
            variant_regulatory = supported(rows, task="VARIANT_REGULATORY_ELEMENT")
            regulatory_gene = supported(rows, task="REGULATORY_ELEMENT_TARGET_GENE", gene_id=gene_id)
            tissue = supported(rows, task="TISSUE_CELL_EXPRESSION", gene_id=gene_id)
            accessibility = [row for row in supported(rows, task="CELL_STATE_ACCESSIBILITY") if row["gene_id"] in {gene_id, "NA", "", "."}]
            pathway = supported(rows, task="PATHWAY", gene_id=gene_id)
            model = supported(rows, task="MODEL_SYSTEM_PERTURBATION", gene_id=gene_id)
            contradictions = sorted({row["contradiction_status"] for row in gene_rows if row["contradiction_status"] != "NONE"})
            if not contradictions:
                contradictions = sorted({row["contradiction_status"] for row in rows if row["contradiction_status"] != "NONE"})
            has_contradiction = bool(contradictions)
            if molecular and variant_regulatory and regulatory_gene and tissue and not has_contradiction:
                confidence = "MULTI_LAYER_HUMAN_GENETIC_SUPPORT"
            elif molecular and not has_contradiction:
                confidence = "COLOCALIZED_MOLECULAR_SUPPORT"
            elif gene_rows and not has_contradiction:
                confidence = "ANNOTATION_ONLY_SUPPORT"
            else:
                confidence = "INSUFFICIENT_OR_CONFLICTING"
            evidence_classes = sorted({row["evidence_class"] for row in gene_rows}, key=lambda value: EVIDENCE_ORDER.index(value) if value in EVIDENCE_ORDER else len(EVIDENCE_ORDER))
            strongest = evidence_classes[0] if evidence_classes else "NO_SUPPORTED_GENE_LEVEL_EVIDENCE"

            evidence_for_edge: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
            for row in supported_rows:
                if row["gene_id"] in {gene_id, "NA", "", "."} or row["search_task"] == "VARIANT_REGULATORY_ELEMENT":
                    evidence_for_edge[(row["supports_edge_from"], row["supports_edge_to"])].append(row)
            supported_edge_count = 2
            missing_edges: list[str] = []
            for edge_order, edge_from, edge_to in REQUIRED_EDGES:
                edge_evidence = evidence_for_edge.get((edge_from, edge_to), [])
                if edge_order == 1:
                    edge_status = "SUPPORTED_BY_UPSTREAM_LOCK"
                    evidence_ids = f"UPSTREAM_COLOC:{signal_key}"
                    evidence_provenance_value = f"fine_mapping_sha256={sha256(args.fine_mapping)}"
                    classes, source_ids = "REPLICATED_HUMAN_GENETICS", "UPSTREAM_FINE_MAPPING_COLOCALIZATION"
                    entity_from, entity_to = task["dataset1_id"], signal_key
                    caveat = "Replicated pair selection and shared-signal model are statistical genetic evidence, not a biological mechanism."
                elif edge_order == 2:
                    edge_status = "SUPPORTED_BY_UPSTREAM_LOCK"
                    evidence_ids = f"UPSTREAM_SHARED_VARIANT:{signal_key}"
                    evidence_provenance_value = f"fine_mapping_sha256={sha256(args.fine_mapping)}"
                    classes, source_ids = "FINE_MAPPING_AND_COLOC_POSTERIOR", "UPSTREAM_FINE_MAPPING_COLOCALIZATION"
                    entity_from, entity_to = signal_key, task["top_shared_variant"]
                    caveat = "The top shared-model variant is probabilistic and is not proven causal."
                elif edge_evidence:
                    edge_status = "SUPPORTED"
                    supported_edge_count += 1
                    evidence_ids = values(edge_evidence, "evidence_id")
                    evidence_provenance_value = f"mechanistic_evidence_sha256={sha256(args.evidence)}"
                    classes, source_ids = values(edge_evidence, "evidence_class"), values(edge_evidence, "source_id")
                    entity_lookup = {
                        "CREDIBLE_VARIANT": values(edge_evidence, "variant_id"),
                        "REGULATORY_ELEMENT": values(edge_evidence, "regulatory_element_id"),
                        "TARGET_GENE": gene_id,
                        "CELL_TYPE_OR_STATE": values(edge_evidence, "cell_type_or_state"),
                        "PATHWAY": values(edge_evidence, "pathway_id"),
                        "EXTERNAL_PHENOTYPE": task["dataset2_id"],
                    }
                    entity_from, entity_to = entity_lookup.get(edge_from, "NONE"), entity_lookup.get(edge_to, "NONE")
                    caveat = values(edge_evidence, "caveat")
                else:
                    edge_status = "MISSING"
                    missing_edges.append(f"{edge_from}->{edge_to}")
                    evidence_ids = "NONE"
                    evidence_provenance_value = "NONE"
                    classes, source_ids = "NONE", "NONE"
                    entity_from, entity_to = "NONE", "NONE"
                    caveat = "No validated evidence supports this required edge; the chain remains incomplete."
                chain_rows.append({
                    "chain_edge_id": f"{signal_key}__{gene_id}__E{edge_order}",
                    "comparison_id": task["comparison_id"], "pair_id": task["pair_id"],
                    "locus_id": task["locus_id"], "signal1": task["signal1"], "signal2": task["signal2"],
                    "gene_id": gene_id, "gene_symbol": gene_symbol, "edge_order": edge_order,
                    "edge_from": edge_from, "edge_to": edge_to, "entity_from_ids": entity_from,
                    "entity_to_ids": entity_to, "edge_status": edge_status, "evidence_ids": evidence_ids,
                    "evidence_provenance": evidence_provenance_value, "evidence_classes": classes,
                    "source_ids": source_ids, "contradiction_status": ";".join(contradictions) if contradictions else "NONE",
                    "caveat": caveat, "claim_limit": contract["claim_limit"],
                })
            for edge, edge_evidence in sorted(evidence_for_edge.items()):
                if edge in {(edge_from, edge_to) for _, edge_from, edge_to in REQUIRED_EDGES}:
                    continue
                chain_rows.append({
                    "chain_edge_id": f"{signal_key}__{gene_id}__AUX__{edge[0]}__{edge[1]}",
                    "comparison_id": task["comparison_id"], "pair_id": task["pair_id"],
                    "locus_id": task["locus_id"], "signal1": task["signal1"], "signal2": task["signal2"],
                    "gene_id": gene_id, "gene_symbol": gene_symbol, "edge_order": "AUXILIARY",
                    "edge_from": edge[0], "edge_to": edge[1], "entity_from_ids": "SEE_EVIDENCE",
                    "entity_to_ids": "SEE_EVIDENCE", "edge_status": "SUPPORTED_AUXILIARY",
                    "evidence_ids": values(edge_evidence, "evidence_id"),
                    "evidence_provenance": f"mechanistic_evidence_sha256={sha256(args.evidence)}",
                    "evidence_classes": values(edge_evidence, "evidence_class"),
                    "source_ids": values(edge_evidence, "source_id"),
                    "contradiction_status": ";".join(contradictions) if contradictions else "NONE",
                    "caveat": values(edge_evidence, "caveat"), "claim_limit": contract["claim_limit"],
                })
            synthesis_rows.append({
                "mechanism_id": f"{signal_key}__{gene_id}", "comparison_id": task["comparison_id"],
                "pair_id": task["pair_id"], "locus_id": task["locus_id"],
                "comparison_type": task["comparison_type"], "signal1": task["signal1"], "signal2": task["signal2"],
                "gene_id": gene_id, "gene_symbol": gene_symbol, "top_shared_variant": task["top_shared_variant"],
                "top_shared_variant_PP_H4": task["top_shared_variant_PP_H4"],
                "upstream_replication_status": "ENTRY_GATE_LOCKED_REPLICATED_PRIORITY_LOCUS",
                "upstream_shared_signal_status": "PRIMARY_AND_PRIOR_ROBUST_COLOC_SUSIE",
                "molecular_qtl_colocalization_support": "SUPPORTED" if molecular else "NOT_SUPPORTED",
                "credible_variant_regulatory_support": "SUPPORTED" if variant_regulatory else "NOT_SUPPORTED",
                "regulatory_element_target_gene_support": "SUPPORTED" if regulatory_gene else "NOT_SUPPORTED",
                "tissue_cell_expression_support": "SUPPORTED" if tissue else "NOT_SUPPORTED",
                "cell_state_accessibility_support": "SUPPORTED" if accessibility else "NOT_SUPPORTED",
                "pathway_support": "SUPPORTED" if pathway else "NOT_SUPPORTED",
                "model_system_perturbation_support": "SUPPORTED" if model else "NOT_SUPPORTED",
                "supported_evidence_classes": ";".join(evidence_classes) if evidence_classes else "NONE",
                "supported_evidence_ids": values(gene_rows, "evidence_id"),
                "source_ids": values(gene_rows, "source_id"),
                "contradiction_status": ";".join(contradictions) if contradictions else "NONE",
                "supported_required_edge_count": supported_edge_count,
                "missing_required_edges": ";".join(missing_edges) if missing_edges else "NONE",
                "strongest_mechanistic_evidence": strongest, "confidence_class": confidence,
                "claim_limit": contract["claim_limit"],
            })

    synthesis_rows.sort(key=lambda row: (
        CONFIDENCE_ORDER[str(row["confidence_class"])], -int(row["supported_required_edge_count"]),
        str(row["pair_id"]), str(row["locus_id"]), str(row["gene_id"]),
    ))
    args.chains_out.parent.mkdir(parents=True, exist_ok=True)
    with args.chains_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=CHAIN_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(chain_rows)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=SYNTHESIS_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(synthesis_rows)
    counts = Counter(str(row["confidence_class"]) for row in synthesis_rows)
    provenance = {
        "schema_version": "1.0.0",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "family_lock_sha256": sha256(args.family_lock),
        "validated_evidence_sha256": sha256(args.evidence),
        "fine_mapping_sha256": sha256(args.fine_mapping),
        "signal_count": len(queue_by_signal),
        "candidate_gene_count": len(synthesis_rows),
        "chain_edge_count": len(chain_rows),
        "confidence_counts": dict(sorted(counts.items())),
        "ranking_policy": "ordered confidence class, then number of explicitly supported required edges; no opaque composite score",
        "chains_output": str(args.chains_out), "chains_output_sha256": sha256(args.chains_out),
        "output": str(args.out), "output_sha256": sha256(args.out),
        "warning": contract["claim_limit"],
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"MECHANISTIC_SYNTHESIS_OK signals={len(queue_by_signal)} genes={len(synthesis_rows)} edges={len(chain_rows)} confidence={dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
