#!/usr/bin/env python3
"""Freeze the exact downstream interpretation task family before results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import downstream_contract

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


def scope_hash(*values: object) -> str:
    return hashlib.sha256(
        json.dumps(values, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def task_row(
    task_id: str, family: str, unit_id: str, source: dict[str, str], input_hash: str,
    *, domain: str = "NA", locus_id: str = "NA", trait_id: str = "NA",
    pair_id: str = "NA", exposure: str = "NA", outcome: str = "NA",
    direction: str = "NA",
) -> dict[str, str]:
    base = f"results/interpretation/runs/{family}/{task_id}"
    return {
        "task_id": task_id, "analysis_family": family, "unit_id": unit_id,
        "method": source["method"], "domain": domain, "locus_id": locus_id,
        "trait_id": trait_id, "pair_id": pair_id, "exposure": exposure,
        "outcome": outcome, "direction": direction, "source_id": source["source_id"],
        "source_release": source["exact_release"], "input_scope_sha256": input_hash,
        "normalized_result_path": f"{base}/normalized.tsv",
        "provenance_path": f"{base}/provenance.json",
    }


def build_tasks(
    policy: dict[str, object], panel: list[dict[str, str]], pairs: list[dict[str, str]],
    loci: list[dict[str, str]], variants: list[dict[str, str]], genes: list[dict[str, str]],
    sources: list[dict[str, str]],
) -> list[dict[str, str]]:
    by_family = {
        family: [row for row in sources if row["analysis_family"] == family]
        for family in ("regulatory", "cell_type", "pathway", "causal")
    }
    variants_by_locus: dict[str, list[dict[str, str]]] = {}
    genes_by_locus: dict[str, list[dict[str, str]]] = {}
    for locus in loci:
        locus_id = locus["locus_id"]
        variants_by_locus[locus_id] = [row for row in variants if row["locus_id"] == locus_id]
        genes_by_locus[locus_id] = [row for row in genes if row["locus_id"] == locus_id]
    tasks: list[dict[str, str]] = []
    for locus in loci:
        locus_id = locus["locus_id"]
        locus_hash = scope_hash(locus, variants_by_locus[locus_id], genes_by_locus[locus_id])
        for source in by_family["regulatory"]:
            for domain in source["coverage_domains"].split(";"):
                task_id = f"REG__{locus_id}__{source['source_id']}__{domain}"
                tasks.append(task_row(
                    task_id, "regulatory", locus_id, source,
                    scope_hash(locus_hash, source["source_id"], source["exact_release"], domain),
                    domain=domain, locus_id=locus_id,
                    pair_id=f"{locus['sleep_trait']}__{locus['non_sleep_trait']}",
                ))
    for trait in panel:
        trait_id = trait["trait_id"]
        for source in by_family["cell_type"]:
            for domain in source["coverage_domains"].split(";"):
                task_id = f"CELL__{trait_id}__{source['source_id']}__{domain}"
                tasks.append(task_row(
                    task_id, "cell_type", trait_id, source,
                    scope_hash(trait, source["source_id"], source["exact_release"], domain),
                    domain=domain, trait_id=trait_id,
                ))
    high_confidence = [row for row in genes if row["evidence_level"] == "HIGH_CONFIDENCE_CONVERGENT_GENE"]
    gene_hash = scope_hash(high_confidence)
    for source in by_family["pathway"]:
        task_id = f"PATH__{source['source_id']}"
        tasks.append(task_row(
            task_id, "pathway", "HIGH_CONFIDENCE_GENES", source,
            scope_hash(gene_hash, source["source_id"], source["exact_release"]),
        ))
    for pair in pairs:
        pair_id = pair["pair_id"]
        for direction in policy["causal_inference"]["directions"]:
            exposure, outcome = (
                (pair["sleep_trait"], pair["non_sleep_trait"])
                if direction == "sleep_to_non_sleep"
                else (pair["non_sleep_trait"], pair["sleep_trait"])
            )
            for source in by_family["causal"]:
                task_id = f"MR__{pair_id}__{direction}__{source['source_id']}"
                tasks.append(task_row(
                    task_id, "causal", f"{pair_id}__{direction}", source,
                    scope_hash(pair, exposure, outcome, direction, source["source_id"], source["exact_release"]),
                    pair_id=pair_id, exposure=exposure, outcome=outcome, direction=direction,
                ))
    return tasks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--preflight", default="results/tables/interpretation_preflight.json")
    parser.add_argument("--out", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--lock-out", default="results/tables/interpretation_task_manifest.lock.json")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    preflight_path = root / args.preflight
    downstream_path = root / "config/downstream_analysis_policy.json"
    downstream_contract.validate_interpretation_preflight(
        root, policy_path, downstream_path, preflight_path,
    )
    registry_path = root / policy["source_registry"]
    references_path = root / policy["method_references"]
    panel_path = root / "config/analysis_panel.tsv"
    pairs_path = root / "results/atlas/trait_pairs.tsv"
    loci_path = root / "results/atlas/loci.tsv"
    variants_path = root / "results/atlas/variants.tsv"
    genes_path = root / "results/atlas/genes.tsv"
    panel, pairs, loci = read_tsv(panel_path), read_tsv(pairs_path), read_tsv(loci_path)
    variants, genes, sources = read_tsv(variants_path), read_tsv(genes_path), read_tsv(registry_path)
    if len(panel) != policy["expected_traits"] or len(pairs) != policy["expected_sleep_non_sleep_pairs"]:
        fail("interpretation inputs differ from the locked trait/pair family")
    if any(row["source_status"] == "CURATION_REQUIRED" for row in sources):
        fail("source registry still contains pre-result curation blockers")
    run_root = root / "results/interpretation/runs"
    if run_root.exists() and any(path.is_file() for path in run_root.rglob("*")):
        fail("interpretation results exist before the task family is locked")
    tasks = build_tasks(policy, panel, pairs, loci, variants, genes, sources)
    fields = policy["task_manifest_fields"]
    if not tasks or any(list(row) != fields for row in tasks):
        fail("interpretation task schema differs from the locked policy")
    task_ids = [row["task_id"] for row in tasks]
    if len(task_ids) != len(set(task_ids)):
        fail("interpretation task IDs are duplicated")
    text = table_text(fields, tasks)
    counts = {
        family: sum(row["analysis_family"] == family for row in tasks)
        for family in ("regulatory", "cell_type", "pathway", "causal")
    }
    expected_causal = (
        policy["expected_sleep_non_sleep_pairs"]
        * len(policy["causal_inference"]["directions"])
        * len(policy["causal_inference"]["methods"])
    )
    if counts["causal"] != expected_causal:
        fail("causal task family is not exact pair-by-direction-by-method coverage")
    lock = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "preflight_sha256": sha256(preflight_path),
        "source_registry_sha256": sha256(registry_path),
        "method_references_sha256": sha256(references_path),
        "input_sha256": {
            str(path.relative_to(root)): sha256(path)
            for path in (panel_path, pairs_path, loci_path, variants_path, genes_path)
        },
        "task_manifest": args.out, "task_manifest_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "task_ids_in_locked_order": task_ids, "task_count": len(tasks),
        "task_counts_by_family": counts,
        "interpretation_results_accessed_before_task_lock": False,
        "script_sha256": downstream_contract.script_hashes(root, "interpretation"),
    }
    lock_text = json.dumps(lock, indent=2, sort_keys=True) + "\n"
    out_path, lock_path = root / args.out, root / args.lock_out
    if args.validate_only:
        if not out_path.is_file() or out_path.read_text(encoding="utf-8") != text:
            fail("interpretation task manifest differs from deterministic recomputation")
        if not lock_path.is_file() or lock_path.read_text(encoding="utf-8") != lock_text:
            fail("interpretation task lock differs from deterministic recomputation")
    else:
        if out_path.exists() or lock_path.exists():
            fail("interpretation task lock already exists; refusing overwrite")
        atomic_text(out_path, text)
        atomic_text(lock_path, lock_text)
    print(
        f"INTERPRETATION_TASKS_LOCKED total={len(tasks)} "
        + " ".join(f"{family}={count}" for family, count in counts.items())
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
