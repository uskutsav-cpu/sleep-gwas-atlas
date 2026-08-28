#!/usr/bin/env python3
"""Exercise the locked mechanistic search plan, provenance validator, and graph synthesis."""

from __future__ import annotations

import csv
import hashlib
import json
import runpy
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "discovery_extension"
EVIDENCE_FIELDS = runpy.run_path(str(EXTENSION / "scripts/39_validate_mechanistic_evidence.py"))["EVIDENCE_FIELDS"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str]) -> str:
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    if result.returncode:
        raise AssertionError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}{result.stderr}")
    return result.stdout + result.stderr


class MechanisticAnnotationContractTest(unittest.TestCase):
    def test_evidence_chain_is_complete_only_where_supported(self) -> None:
        with tempfile.TemporaryDirectory(prefix="sleep_gwas_mechanism_synthetic_") as temporary:
            temp = Path(temporary)
            fine_mapping = temp / "fine_mapping_colocalization.tsv"
            fine_mapping_row = {
                "comparison_id": "sleep_test__external_test__locus_1__SLEEP_EQTL__GENE1",
                "pair_id": "sleep_test__external_test", "locus_id": "locus_1",
                "comparison_type": "SLEEP_EQTL", "dataset1_id": "sleep_test", "dataset2_id": "GENE1_eqtl",
                "molecular_feature_id": "GENE1", "signal1": "L1", "signal2": "L1",
                "top_shared_variant": "rs61", "top_shared_variant_PP_H4": "0.91",
                "PP_H4": "0.98", "PP_H4_over_PP_H3": "49", "prior_role": "PRIMARY",
                "coloc_method": "COLOC_SUSIE", "primary_shared_signal_rule_pass": "True",
                "prior_robust": "True", "colocalization_interpretation": "SHARED_SIGNAL_MODEL_SUPPORTED",
                "analysis_status": "COLOC_SUSIE_COMPLETE",
            }
            write_tsv(fine_mapping, list(fine_mapping_row), [fine_mapping_row])
            fine_mapping_provenance = temp / "fine_mapping.json"
            fine_mapping_provenance.write_text(json.dumps({"output_sha256": sha256(fine_mapping)}) + "\n", encoding="utf-8")
            fine_mapping_lock = temp / "fine_mapping.lock.json"
            fine_mapping_lock.write_text("{}\n", encoding="utf-8")

            readiness = temp / "readiness.tsv"
            readiness_provenance = temp / "readiness.json"
            preflight = run([
                "python3", str(EXTENSION / "scripts/37_mechanism_preflight.py"),
                "--fine-mapping", str(fine_mapping), "--fine-mapping-lock", str(fine_mapping_lock),
                "--out", str(readiness), "--provenance-out", str(readiness_provenance),
            ])
            self.assertIn("READY_FOR_RESULT_FREE_SOURCE_PLAN", preflight)

            queue = temp / "queue.tsv"
            family_lock = temp / "family.lock.json"
            queue_provenance = temp / "queue.json"
            queue_output = run([
                "python3", str(EXTENSION / "scripts/38_prepare_mechanism_queue.py"),
                "--fine-mapping", str(fine_mapping), "--fine-mapping-provenance", str(fine_mapping_provenance),
                "--out", str(queue), "--family-lock-out", str(family_lock),
                "--provenance-out", str(queue_provenance),
            ])
            self.assertIn("MECHANISM_QUEUE_OK signals=1 tasks=7", queue_output)

            queue_rows = read_tsv(queue)
            source_for_task = {
                "MOLECULAR_QTL": "EQTL_CATALOGUE",
                "VARIANT_REGULATORY_ELEMENT": "SCREEN_CCRE",
                "REGULATORY_ELEMENT_TARGET_GENE": "OPEN_TARGETS_PLATFORM",
                "TISSUE_CELL_EXPRESSION": "SINGLE_CELL_EXPRESSION_ATLAS",
                "CELL_STATE_ACCESSIBILITY": "ENCODE_DCC",
                "PATHWAY": "REACTOME",
                "MODEL_SYSTEM_PERTURBATION": "BIOSTUDIES_GEO",
            }
            class_for_task = {
                "MOLECULAR_QTL": "COLOCALIZED_MOLECULAR_QTL",
                "VARIANT_REGULATORY_ELEMENT": "FINE_MAPPED_REGULATORY_EVIDENCE",
                "REGULATORY_ELEMENT_TARGET_GENE": "FINE_MAPPED_REGULATORY_EVIDENCE",
                "TISSUE_CELL_EXPRESSION": "TISSUE_CELL_EXPRESSION_OR_ENRICHMENT",
                "CELL_STATE_ACCESSIBILITY": "TISSUE_CELL_EXPRESSION_OR_ENRICHMENT",
                "PATHWAY": "PATHWAY_OR_NETWORK_ANNOTATION",
                "MODEL_SYSTEM_PERTURBATION": "NO_SUPPORTED_EVIDENCE",
            }
            edge_for_task = {
                "MOLECULAR_QTL": ("SHARED_SIGNAL_MODEL", "TARGET_GENE"),
                "VARIANT_REGULATORY_ELEMENT": ("CREDIBLE_VARIANT", "REGULATORY_ELEMENT"),
                "REGULATORY_ELEMENT_TARGET_GENE": ("REGULATORY_ELEMENT", "TARGET_GENE"),
                "TISSUE_CELL_EXPRESSION": ("TARGET_GENE", "CELL_TYPE_OR_STATE"),
                "CELL_STATE_ACCESSIBILITY": ("REGULATORY_ELEMENT", "CELL_TYPE_OR_STATE"),
                "PATHWAY": ("TARGET_GENE", "PATHWAY"),
                "MODEL_SYSTEM_PERTURBATION": ("NA", "NA"),
            }
            evidence_rows: list[dict[str, str]] = []
            for index, task in enumerate(queue_rows, start=1):
                snapshot = temp / f"snapshot_{index}.tsv"
                snapshot.write_text(f"SYNTHETIC TEST SNAPSHOT ONLY\n{task['search_task']}\n", encoding="utf-8")
                outcome = "NO_EVIDENCE_FOUND" if task["search_task"] == "MODEL_SYSTEM_PERTURBATION" else "SUPPORTED"
                edge_from, edge_to = edge_for_task[task["search_task"]]
                row = {field: "NA" for field in EVIDENCE_FIELDS}
                row.update({
                    "evidence_id": f"SYNTHETIC_E{index}", "search_task_id": task["search_task_id"],
                    "comparison_id": task["comparison_id"], "pair_id": task["pair_id"],
                    "locus_id": task["locus_id"], "signal1": task["signal1"], "signal2": task["signal2"],
                    "search_task": task["search_task"], "search_outcome": outcome,
                    "evidence_class": class_for_task[task["search_task"]],
                    "source_id": source_for_task[task["search_task"]],
                    "exact_dataset_release": "SYNTHETIC_RELEASE", "accession": f"SYNTHETIC_ACC_{index}",
                    "source_url": f"synthetic://mechanism/{index}", "source_access_date": "2026-08-28",
                    "source_license_or_terms": "SYNTHETIC_ONLY", "snapshot_path": str(snapshot),
                    "snapshot_sha256": sha256(snapshot), "primary_citation_title": "SYNTHETIC FIXTURE — NOT A REAL CITATION",
                    "species": "Homo sapiens", "biological_context": "SYNTHETIC_TISSUE",
                    "gene_id": "ENSG_SYNTHETIC_1" if outcome == "SUPPORTED" and task["search_task"] not in {"VARIANT_REGULATORY_ELEMENT", "CELL_STATE_ACCESSIBILITY"} else "NA",
                    "gene_symbol": "GENE1" if outcome == "SUPPORTED" and task["search_task"] not in {"VARIANT_REGULATORY_ELEMENT", "CELL_STATE_ACCESSIBILITY"} else "NA",
                    "variant_id": "rs61" if task["search_task"] == "VARIANT_REGULATORY_ELEMENT" else "NA",
                    "regulatory_element_id": "EH_SYNTHETIC_1" if task["search_task"] in {"VARIANT_REGULATORY_ELEMENT", "REGULATORY_ELEMENT_TARGET_GENE", "CELL_STATE_ACCESSIBILITY"} else "NA",
                    "cell_type_or_state": "SYNTHETIC_NEURON" if task["search_task"] in {"TISSUE_CELL_EXPRESSION", "CELL_STATE_ACCESSIBILITY"} else "NA",
                    "pathway_id": "R-HSA-SYNTHETIC" if task["search_task"] == "PATHWAY" else "NA",
                    "pathway_name": "Synthetic pathway" if task["search_task"] == "PATHWAY" else "NA",
                    "colocalization_interpretation": "SHARED_SIGNAL_MODEL_SUPPORTED" if task["search_task"] == "MOLECULAR_QTL" else "NA",
                    "supports_edge_from": edge_from, "supports_edge_to": edge_to,
                    "contradiction_status": "NONE", "evidence_summary": f"SYNTHETIC {outcome} outcome",
                    "caveat": "SYNTHETIC TEST EVIDENCE; NOT A BIOLOGICAL RESULT",
                    "claim_limit": "SYNTHETIC TEST ONLY; NO CAUSAL CLAIM", "curator": "synthetic_test",
                    "curation_date": "2026-08-28",
                })
                evidence_rows.append(row)
            evidence_input = temp / "evidence_input.tsv"
            write_tsv(evidence_input, EVIDENCE_FIELDS, evidence_rows)

            validated_evidence = temp / "mechanistic_evidence.tsv"
            evidence_provenance = temp / "evidence_validation.json"
            validation = run([
                "python3", str(EXTENSION / "scripts/39_validate_mechanistic_evidence.py"),
                "--queue", str(queue), "--family-lock", str(family_lock), "--evidence", str(evidence_input),
                "--out", str(validated_evidence), "--provenance-out", str(evidence_provenance),
            ])
            self.assertIn("MECHANISTIC_EVIDENCE_VALID tasks=7 rows=7", validation)

            chains = temp / "mechanistic_chains.tsv"
            synthesis = temp / "mechanistic_synthesis.tsv"
            synthesis_provenance = temp / "synthesis.json"
            synthesis_output = run([
                "python3", str(EXTENSION / "scripts/40_synthesize_mechanisms.py"),
                "--queue", str(queue), "--family-lock", str(family_lock),
                "--evidence", str(validated_evidence), "--evidence-provenance", str(evidence_provenance),
                "--fine-mapping", str(fine_mapping), "--fine-mapping-provenance", str(fine_mapping_provenance),
                "--chains-out", str(chains), "--out", str(synthesis),
                "--provenance-out", str(synthesis_provenance),
            ])
            self.assertIn("MECHANISTIC_SYNTHESIS_OK signals=1 genes=1", synthesis_output)
            synthesis_rows = read_tsv(synthesis)
            chain_rows = read_tsv(chains)
            self.assertEqual(synthesis_rows[0]["confidence_class"], "MULTI_LAYER_HUMAN_GENETIC_SUPPORT")
            self.assertEqual(synthesis_rows[0]["supported_required_edge_count"], "5")
            self.assertEqual(
                set(synthesis_rows[0]["missing_required_edges"].split(";")),
                {"CELL_TYPE_OR_STATE->PATHWAY", "PATHWAY->EXTERNAL_PHENOTYPE"},
            )
            self.assertEqual(sum(row["edge_status"] == "MISSING" for row in chain_rows), 2)
            self.assertTrue(all("caus" in row["claim_limit"].lower() for row in synthesis_rows))


if __name__ == "__main__":
    unittest.main()
