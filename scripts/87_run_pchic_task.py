#!/usr/bin/env python3
"""Materialize one locked PCHi-C task, including explicit unavailable domains."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--cache")
    parser.add_argument("--cache-provenance")
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
    if (
        task["analysis_family"] != "regulatory" or task["method"] != "chromatin_contact"
        or task["source_id"] != policy["pchic_2016"]["source_id"]
    ):
        fail("task is not supported by the PCHi-C task adapter")
    source_manifest_path = root / policy["pchic_2016"]["component_manifest"]
    if (
        not source_manifest_path.is_file()
        or sha256(source_manifest_path) != policy["pchic_2016"]["component_manifest_sha256"]
    ):
        fail("PCHi-C source manifest differs from its policy pin")
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    if task["source_release"] != source_manifest["release"]:
        fail("PCHi-C task release differs from the source manifest")
    component_hashes: dict[str, str] = {}
    for component in source_manifest["components"]:
        path = root / component["path"]
        if (
            not path.is_file() or path.stat().st_size != component["bytes"]
            or sha256(path) != component["sha256"]
        ):
            fail(f"PCHi-C source component differs from its pin: {component['component_id']}")
        component_hashes[component["component_id"]] = sha256(path)

    fields = policy["regulatory_mapping"]["canonical_fields"]
    rows: list[dict[str, str]] = []
    cache_details: dict[str, object] | None = None
    if task["domain"] in policy["pchic_2016"]["unavailable_domains"]:
        status = "NOT_APPLICABLE"
        reason = (
            f"Javierre 2016 PCHi-C covers 17 primary hematopoietic cell types only; "
            f"the locked {task['domain']} domain is absent from this source"
        )
    elif task["domain"] not in policy["pchic_2016"]["available_domains"]:
        fail("PCHi-C task domain is absent from the frozen availability partition")
    else:
        cache_path = root / (args.cache or policy["pchic_2016"]["cache_path"])
        cache_provenance_path = root / (
            args.cache_provenance or policy["pchic_2016"]["cache_provenance_path"]
        )
        if not cache_path.is_file() or not cache_provenance_path.is_file():
            fail("PCHi-C overlap cache and provenance must be complete before an immune task runs")
        cache_provenance = json.loads(cache_provenance_path.read_text(encoding="utf-8"))
        variants_path, genes_path = root / "results/atlas/variants.tsv", root / "results/atlas/genes.tsv"
        if (
            cache_provenance.get("policy_sha256") != sha256(policy_path)
            or cache_provenance.get("source_id") != task["source_id"]
            or cache_provenance.get("source_release") != task["source_release"]
            or cache_provenance.get("variants_sha256") != sha256(variants_path)
            or cache_provenance.get("genes_sha256") != sha256(genes_path)
            or cache_provenance.get("cache_sha256") != sha256(cache_path)
        ):
            fail("PCHi-C overlap cache differs from its checksum-bound inputs or task source")
        cache_fields, cached = read_tsv(cache_path)
        if cache_fields != fields:
            fail("PCHi-C overlap cache differs from the locked regulatory schema")
        rows = [
            row for row in cached
            if row["locus_id"] == task["locus_id"] and row["context_domain"] == task["domain"]
            and row["source_dataset"] == task["source_id"]
            and row["source_version"] == task["source_release"]
        ]
        rows.sort(key=lambda row: (
            row["variant_id"], row["regulatory_element_id"], row["biosample"],
            row["target_gene_id"], row["regulatory_evidence_id"],
        ))
        status = "COMPLETED" if rows else "NO_EVIDENCE_FOUND"
        reason = (
            f"Pinned PCHi-C cache produced {len(rows)} qualified promoter-contact evidence rows"
            if rows else
            "No eligible fine-mapped variant overlapped a CHiCAGO>=5 contact to an exactly matched supported same-locus gene in the immune source"
        )
        cache_details = {
            "cache_path": str(cache_path), "cache_sha256": sha256(cache_path),
            "cache_provenance_path": str(cache_provenance_path),
            "cache_provenance_sha256": sha256(cache_provenance_path),
        }
    payload = table_text(fields, rows)
    out_path = root / args.out
    atomic_text(out_path, payload)
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"], "terminal_status": status, "terminal_reason": reason,
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"], "source_id": task["source_id"],
        "source_release": task["source_release"], "source_manifest_sha256": sha256(source_manifest_path),
        "source_component_sha256": component_hashes, "cache": cache_details,
        "result_rows": len(rows), "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "claim_limit": policy["pchic_2016"]["claim_limit"],
    }
    provenance_path = root / args.provenance_out
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"PCHIC_TASK_OK task={task['task_id']} status={status} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
