#!/usr/bin/env python3
"""Validate and index every checksum-locked phi-enabled GTEx-v8 model context."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import sqlite3
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import molecular_contract

REGISTRY_FIELDS = molecular_contract.MODEL_REGISTRY_FIELDS
PHI_EXCLUSION_FIELDS = molecular_contract.PHI_EXCLUSION_FIELDS
VARIANT_FIELDS = molecular_contract.MODEL_VARIANT_FIELDS
VARIANT_ID = re.compile(r"^(?:chr)?([1-9]|1[0-9]|2[0-2])[_:]([0-9]+)[_:]([ACGT])[_:]([ACGT])(?:[_:]b38)?$", re.I)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def gzip_text(path: Path, text: str) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def inspect_model(path: Path) -> tuple[str, int, int, int, list[str], float, float, dict[str, tuple[int, int, str, str]]]:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        columns = [row[1] for row in connection.execute("PRAGMA table_info(weights)")]
        required = {"gene", "weight", "ref_allele", "eff_allele"}
        if not required.issubset(columns):
            fail(f"PredictDB weights schema is incomplete: {path}")
        candidate_keys = [key for key in ("varID", "varid", "rsid") if key in columns]
        if not candidate_keys:
            fail(f"PredictDB weights table has no supported variant key: {path}")
        chosen = None
        for candidate in candidate_keys:
            sample = [row[0] for row in connection.execute(f'SELECT DISTINCT "{candidate}" FROM weights LIMIT 1000')]
            if sample and sum(VARIANT_ID.fullmatch(str(value)) is not None for value in sample) / len(sample) >= 0.95:
                chosen = candidate; break
        if chosen is None:
            fail(f"PredictDB model lacks parseable GRCh38 variant identifiers: {path}")
        gene_count = int(connection.execute("SELECT COUNT(DISTINCT gene) FROM weights").fetchone()[0])
        weight_count = int(connection.execute("SELECT COUNT(*) FROM weights").fetchone()[0])
        extra_columns = [row[1] for row in connection.execute("PRAGMA table_info(extra)")]
        if not {"gene", "phi"}.issubset(extra_columns):
            fail(f"PredictDB model lacks gene-specific phi variance-control inputs: {path}")
        missing_phi = [str(row[0]) for row in connection.execute(
            "SELECT DISTINCT w.gene FROM weights w "
            "LEFT JOIN extra e ON w.gene=e.gene AND e.phi IS NOT NULL WHERE e.gene IS NULL"
            " ORDER BY w.gene"
        )]
        phi_values = [float(row[0]) for row in connection.execute(
            "SELECT DISTINCT e.phi FROM extra e JOIN weights w ON e.gene=w.gene"
        )]
        if not phi_values or any(not math.isfinite(value) for value in phi_values):
            fail(f"PredictDB model contains invalid phi values: {path}")
        phi_gene_count = int(connection.execute(
            "SELECT COUNT(DISTINCT e.gene) FROM extra e JOIN weights w ON e.gene=w.gene WHERE e.phi IS NOT NULL"
        ).fetchone()[0])
        variants: dict[str, tuple[int, int, str, str]] = {}
        for variant_id, non_effect, effect in connection.execute(
            f'SELECT DISTINCT "{chosen}", ref_allele, eff_allele FROM weights'
        ):
            match = VARIANT_ID.fullmatch(str(variant_id))
            if match is None:
                fail(f"unparseable selected model variant ID: {path}/{variant_id}")
            chromosome, position = int(match.group(1)), int(match.group(2))
            ref, alt = match.group(3).upper(), match.group(4).upper()
            if ref == alt or {str(non_effect).upper(), str(effect).upper()} != {ref, alt}:
                fail(f"model alleles disagree with GRCh38 variant ID: {path}/{variant_id}")
            value = (chromosome, position, ref, alt)
            if str(variant_id) in variants and variants[str(variant_id)] != value:
                fail(f"model variant identity is inconsistent: {path}/{variant_id}")
            variants[str(variant_id)] = value
        if phi_gene_count + len(missing_phi) != gene_count:
            fail(f"phi coverage accounting differs from modeled-gene count: {path}")
        return chosen, gene_count, weight_count, phi_gene_count, missing_phi, min(phi_values), max(phi_values), variants
    finally:
        connection.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--registry-out", default="results/tables/twas_model_registry.tsv")
    parser.add_argument("--variant-dir", default="results/tables/twas_model_variants")
    parser.add_argument("--phi-exclusions-out", default="results/tables/twas_model_phi_exclusions.tsv")
    parser.add_argument("--lock-out", default="results/tables/twas_model_registry.lock.json")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    if not args.materialize:
        fail("model validation/indexing requires explicit --materialize")
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    twas, source = policy["twas"], policy["twas"]["phi_model_source"]
    inventory_path, inventory_lock_path = root / source["inventory_path"], root / source["inventory_lock_path"]
    download_lock_path = root / source["download_lock_path"]
    inventory_fields, inventory = read_tsv(inventory_path)
    expected_inventory_fields = [
        "file_id", "filename", "context", "file_role", "bytes", "download_url",
        "results_accessed_before_lock",
    ]
    if inventory_fields != expected_inventory_fields:
        fail("phi-model inventory schema drifted")
    inventory_lock = json.loads(inventory_lock_path.read_text(encoding="utf-8"))
    download_lock = json.loads(download_lock_path.read_text(encoding="utf-8"))
    if (
        inventory_lock.get("inventory_sha256") != sha256(inventory_path)
        or inventory_lock.get("policy_sha256") != sha256(policy_path)
        or inventory_lock.get("file_ids_in_locked_order") != [row["file_id"] for row in inventory]
        or download_lock.get("inventory_sha256") != sha256(inventory_path)
        or download_lock.get("inventory_lock_sha256") != sha256(inventory_lock_path)
        or download_lock.get("policy_sha256") != sha256(policy_path)
    ):
        fail("phi-model inventory/download evidence differs from lock")
    model_dir = root / source["install_dir"]
    installed: dict[str, Path] = {}
    for row in inventory:
        path = model_dir / row["filename"]
        evidence = download_lock.get("files", {}).get(row["filename"], {})
        if (
            not path.is_file() or path.stat().st_size != int(row["bytes"])
            or evidence.get("file_id") != row["file_id"] or evidence.get("bytes") != int(row["bytes"])
            or evidence.get("sha256") != sha256(path)
        ):
            fail(f"phi-model source file differs from download lock: {row['filename']}")
        installed[row["filename"]] = path
    by_context: dict[str, dict[str, dict[str, str]]] = {}
    for row in inventory:
        by_context.setdefault(row["context"], {})[row["file_role"]] = row
    if len(by_context) != int(source["expected_context_count"]) or any(set(files) != {"MODEL_DB", "COVARIANCE"} for files in by_context.values()):
        fail("installed phi-model family is incomplete")
    registry_out, variant_dir = root / args.registry_out, root / args.variant_dir
    phi_exclusions_out, lock_out = root / args.phi_exclusions_out, root / args.lock_out
    if registry_out.exists() or variant_dir.exists() or phi_exclusions_out.exists() or lock_out.exists():
        fail("TWAS model-index outputs already exist; refusing overwrite")
    staging_parent = root / "work/twas_model_index"; staging_parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="phi_models.", dir=staging_parent))
    try:
        family = twas["model_families"][0]
        registry: list[dict[str, object]] = []
        phi_exclusions: list[dict[str, object]] = []
        union: dict[str, tuple[int, int, str, str]] = {}
        for context in sorted(by_context):
            files = by_context[context]
            model_row, covariance_row = files["MODEL_DB"], files["COVARIANCE"]
            model, covariance = installed[model_row["filename"]], installed[covariance_row["filename"]]
            snp_key, gene_count, weight_count, phi_gene_count, missing_phi, phi_min, phi_max, variants = inspect_model(model)
            for variant_id, value in variants.items():
                if variant_id in union and union[variant_id] != value:
                    fail(f"PredictDB family has inconsistent variant identity: {variant_id}")
                union[variant_id] = value
            model_id = f"{family}__{context}"
            registry.append({
                "model_id": model_id, "model_family": family, "modality": "TWAS",
                "context": context, "model_db_path": str(model.relative_to(root)),
                "model_db_sha256": sha256(model), "covariance_path": str(covariance.relative_to(root)),
                "covariance_sha256": sha256(covariance), "model_snp_key": snp_key,
                "gene_count": gene_count, "weight_count": weight_count, "phi_gene_count": phi_gene_count,
                "phi_missing_gene_count": len(missing_phi),
                "phi_min": format(phi_min, ".15g"), "phi_max": format(phi_max, ".15g"),
                "source_release": source["release_name"], "model_file_id": model_row["file_id"],
                "covariance_file_id": covariance_row["file_id"],
            })
            phi_exclusions.extend({
                "model_id": model_id, "model_family": family, "context": context, "gene_id": gene,
                "exclusion_reason": "MODEL_FEATURE_MISSING_PHI",
            } for gene in missing_phi)
        if not union or len(registry) != int(source["expected_context_count"]):
            fail("phi-model index is empty or context-incomplete")
        variant_rows = [
            {"model_variant_id": variant_id, "chromosome_grch38": value[0], "position_grch38": value[1], "ref": value[2], "alt": value[3]}
            for variant_id, value in sorted(union.items(), key=lambda item: (item[1][0], item[1][1], item[0]))
        ]
        staged_registry, staged_variants = staging / "registry.tsv", staging / "variants.tsv.gz"
        staged_exclusions = staging / "phi_exclusions.tsv"
        staged_registry.write_text(table_text(REGISTRY_FIELDS, registry), encoding="utf-8")
        gzip_text(staged_variants, table_text(VARIANT_FIELDS, variant_rows))
        staged_exclusions.write_text(table_text(PHI_EXCLUSION_FIELDS, phi_exclusions), encoding="utf-8")
        registry_out.parent.mkdir(parents=True, exist_ok=True); variant_dir.mkdir(parents=True)
        os.replace(staged_registry, registry_out)
        final_variants = variant_dir / f"{family}.tsv.gz"; os.replace(staged_variants, final_variants)
        os.replace(staged_exclusions, phi_exclusions_out)
        lock = {
            "schema_version": "atlas-v1.0-twas-models.2",
            "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "results_accessed": False, "model_count": len(registry),
            "model_ids_in_locked_order": [row["model_id"] for row in registry],
            "variant_count_by_family": {family: len(variant_rows)}, "registry_sha256": sha256(registry_out),
            "variant_registry_sha256": {family: sha256(final_variants)},
            "phi_exclusion_count": len(phi_exclusions), "phi_exclusions_sha256": sha256(phi_exclusions_out),
            "inventory_sha256": sha256(inventory_path), "inventory_lock_sha256": sha256(inventory_lock_path),
            "download_lock_sha256": sha256(download_lock_path), "policy_sha256": sha256(policy_path),
            "script_sha256": molecular_contract.script_hashes(
                Path(__file__).resolve().parent.parent, "twas",
            ),
            "claim_limit": policy["claim_limit"],
        }
        lock_out.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        staging.rmdir()
        print(f"TWAS_MODEL_INDEX_OK models={len(registry)} variants={len(variant_rows)} phi_exclusions={len(phi_exclusions)}")
        return 0
    finally:
        if staging.exists(): shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
