#!/usr/bin/env python3
"""Freeze every molecular-source query before any interval result is accessed."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from liftover_chain import load_chain


FIELDS = [
    "search_task_id", "locus_id", "pair_id", "sleep_trait", "non_sleep_trait",
    "chromosome_grch37", "start_grch37", "end_grch37", "chromosome_grch38",
    "start_grch38", "end_grch38", "source_family_id", "study_id", "dataset_id",
    "modality", "quant_method", "context", "sample_size", "source_build", "query_mode",
    "source_url", "source_index_url", "exact_release", "planned_outcome",
    "results_accessed_before_lock",
]
SAFE_ID = re.compile(r"[^A-Za-z0-9_.-]+")


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
        fail(f"missing real non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def text_tsv(rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def atomically_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def https_ebi(url: str) -> str:
    if url.startswith("ftp://ftp.ebi.ac.uk/"):
        return "https://ftp.ebi.ac.uk/" + url.removeprefix("ftp://ftp.ebi.ac.uk/")
    if not url.startswith("https://"):
        fail(f"source path cannot be converted to HTTPS: {url}")
    return url


def clean(value: str) -> str:
    return SAFE_ID.sub("_", value).strip("_")


def dataset_family(root: Path, policy: dict) -> list[dict[str, str]]:
    metadata_path = root / policy["metadata_assets"][0]["path"]
    uniform_path = root / policy["metadata_assets"][1]["path"]
    imported_path = root / policy["metadata_assets"][2]["path"]
    _, metadata = read_tsv(metadata_path)
    _, uniform = read_tsv(uniform_path)
    meta_by_id = {row["dataset_id"]: row for row in metadata}
    if len(meta_by_id) != len(metadata) or len(uniform) != len(metadata):
        fail("release-7 metadata/path registries are not one row per dataset")
    method_map = policy["search_family"]["uniform_quant_methods"]
    family_by_method = {
        "ge": "EQTL_CATALOGUE_R7_GE",
        "leafcutter": "EQTL_CATALOGUE_R7_LEAFCUTTER",
        "aptamer": "EQTL_CATALOGUE_R7_SUN2018",
    }
    selected: list[dict[str, str]] = []
    for row in uniform:
        meta = meta_by_id.get(row["dataset_id"])
        if meta is None:
            fail(f"release-7 path has no metadata row: {row['dataset_id']}")
        for field in ("study_id", "study_label", "sample_group", "sample_size", "quant_method"):
            if row[field] != meta[field]:
                fail(f"release-7 metadata/path drift for {row['dataset_id']}/{field}")
        method = row["quant_method"]
        if method not in method_map:
            continue
        sample_size = int(row["sample_size"])
        if sample_size < policy["search_family"]["minimum_sample_size"]:
            continue
        selected.append({
            "source_family_id": family_by_method[method], "study_id": row["study_id"],
            "dataset_id": row["dataset_id"], "modality": method_map[method],
            "quant_method": method,
            "context": "|".join((row["study_label"], row["sample_group"], meta["tissue_label"], meta["condition_label"])),
            "sample_size": row["sample_size"], "source_build": "GRCh38",
            "query_mode": "TABIX_GRCH38_INTERVAL", "source_url": https_ebi(row["ftp_path"]),
            "source_index_url": https_ebi(row["ftp_path"]) + ".tbi", "exact_release": "eQTL_Catalogue_r7_June_2024",
        })
    _, imported = read_tsv(imported_path)
    for row in imported:
        if row["study"] != "GTEx_V8" or row["quant_method"] != "ge":
            fail("imported registry contains a non-GTEx-v8 gene-expression row")
        selected.append({
            "source_family_id": "EQTL_CATALOGUE_GTEX_V8_IMPORTED", "study_id": "GTEx_V8",
            "dataset_id": f"GTEx_V8__{row['qtl_group']}", "modality": "eQTL", "quant_method": "ge",
            "context": "|".join((row["qtl_group"], row["tissue_ontology_id"], row["tissue_label"])),
            "sample_size": "LOCK_FROM_GTEX_V8_TISSUE_METADATA", "source_build": "GRCh38",
            "query_mode": "TABIX_GRCH38_INTERVAL", "source_url": https_ebi(row["ftp_path"]),
            "source_index_url": https_ebi(row["ftp_path"]) + ".tbi", "exact_release": "GTEx_Analysis_v8_imported",
        })
    selected.extend([
        {
            "source_family_id": "FENLAND_CIS_PQTL_V1", "study_id": "Pietzner_2021_Fenland",
            "dataset_id": "Zenodo_6787142_Fenland_cis_pQTL_v1", "modality": "pQTL", "quant_method": "Olink",
            "context": "plasma|Fenland|N=485", "sample_size": "485", "source_build": "GRCh37",
            "query_mode": "LOCAL_ARCHIVE_GRCH37_INTERVAL", "source_url": "https://zenodo.org/api/records/6787142/files/Summary.statistics.cis.Olink.proteins.Fenland.20220107.txt.gz/content",
            "source_index_url": "NA", "exact_release": "Zenodo_6787142_v1",
        },
        {
            "source_family_id": "PSYCHENCODE_R3_BRAIN", "study_id": "PsychENCODE",
            "dataset_id": "Synapse_syn4921369_release3_discovery", "modality": "eQTL_or_sQTL", "quant_method": "per_exact_file",
            "context": "brain|PsychENCODE_release3", "sample_size": "LOCK_PER_FILE", "source_build": "MIXED_LOCK_PER_FILE",
            "query_mode": "SYNAPSE_EXACT_FILE_DISCOVERY", "source_url": "https://www.synapse.org/Synapse:syn4921369",
            "source_index_url": "NA", "exact_release": "PsychENCODE_Release_3",
        },
    ])
    keys = [(row["source_family_id"], row["dataset_id"]) for row in selected]
    if len(keys) != len(set(keys)):
        fail("duplicate molecular dataset identity after source-family construction")
    counts = {method: sum(row["quant_method"] == method for row in selected) for method in method_map}
    if counts != {"ge": 280, "leafcutter": 113, "aptamer": 1}:
        fail(f"release-7 selected dataset counts drifted: {counts}")
    if sum(row["source_family_id"] == "EQTL_CATALOGUE_GTEX_V8_IMPORTED" for row in selected) != 49:
        fail("GTEx-v8 imported search family is not exactly 49 tissues")
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--sources", default="config/molecular_source_registry.tsv")
    parser.add_argument("--preflight", default="results/tables/molecular_preflight.json")
    parser.add_argument("--loci", default="results/atlas/loci.tsv")
    parser.add_argument("--variants", default="results/atlas/variants.tsv")
    parser.add_argument("--fine-mapping-provenance", default="results/atlas/fine_mapping.provenance.json")
    parser.add_argument("--out", default="results/tables/molecular_search_plan.tsv")
    parser.add_argument("--lock-out", default="results/tables/molecular_search_plan.lock.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, sources_path = root / args.policy, root / args.sources
    preflight_path, loci_path, variants_path = root / args.preflight, root / args.loci, root / args.variants
    provenance_path = root / args.fine_mapping_provenance
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    if not preflight.get("code_ready") or preflight.get("policy_sha256") != sha256(policy_path) or preflight.get("sources_sha256") != sha256(sources_path):
        fail("molecular preflight is not code-ready or differs from current contracts")
    _, loci = read_tsv(loci_path)
    _, variants = read_tsv(variants_path)
    fm = json.loads(provenance_path.read_text(encoding="utf-8"))
    if fm.get("outputs", {}).get("results/atlas/loci.tsv") != sha256(loci_path) or fm.get("outputs", {}).get("results/atlas/variants.tsv") != sha256(variants_path):
        fail("canonical loci or variants differ from fine-mapping provenance")
    if not loci or any(row["analysis_tier"] != "PRIMARY_PHASE1" for row in loci):
        fail("molecular search entry family is not the nonempty primary fine-mapped locus family")
    variants_by_locus = {row["locus_id"] for row in variants}
    if not {row["locus_id"] for row in loci}.issubset(variants_by_locus):
        fail("not every molecular-search locus has fine-mapped variants")
    chain_spec = policy["reference_build"]
    chain, chain_provenance = load_chain(
        root / chain_spec["chain_path"], chain_spec["chain_sha256"], chain_spec["chain_bytes"],
    )
    datasets = dataset_family(root, policy)
    rows: list[dict[str, object]] = []
    for locus in loci:
        chromosome = int(locus["chromosome"])
        start, end = int(locus["start_bp"]), int(locus["end_bp"])
        start_status, mapped_start = chain.map_point(chromosome, start)
        end_status, mapped_end = chain.map_point(chromosome, end)
        if start_status != "mapped" or end_status != "mapped" or mapped_start is None or mapped_end is None:
            fail(f"locus endpoints do not both map to GRCh38: {locus['locus_id']}")
        if mapped_start[0] != mapped_end[0] or mapped_start[1] > mapped_end[1] or mapped_start[2] != mapped_end[2]:
            fail(f"locus endpoint mapping is not a single forward interval: {locus['locus_id']}")
        for dataset in datasets:
            use_grch38 = dataset["source_build"] == "GRCh38"
            task_id = "MQS__" + clean(locus["locus_id"]) + "__" + clean(dataset["source_family_id"]) + "__" + clean(dataset["dataset_id"])
            record = {
                "search_task_id": task_id, "locus_id": locus["locus_id"],
                "pair_id": f"{locus['sleep_trait']}__{locus['non_sleep_trait']}",
                "sleep_trait": locus["sleep_trait"], "non_sleep_trait": locus["non_sleep_trait"],
                "chromosome_grch37": chromosome, "start_grch37": start, "end_grch37": end,
                "chromosome_grch38": mapped_start[0], "start_grch38": mapped_start[1], "end_grch38": mapped_end[1],
                **dataset, "planned_outcome": "SOURCE_QUERY_REQUIRED",
                "results_accessed_before_lock": "NO",
            }
            if not use_grch38 and dataset["source_build"] == "GRCh37":
                record["query_mode"] = "LOCAL_ARCHIVE_GRCH37_INTERVAL"
            rows.append({field: record[field] for field in FIELDS})
    task_ids = [str(row["search_task_id"]) for row in rows]
    if len(task_ids) != len(set(task_ids)):
        fail("duplicate molecular search task identity")
    plan_text = text_tsv(rows)
    plan_sha = hashlib.sha256(plan_text.encode()).hexdigest()
    out, lock_out = root / args.out, root / args.lock_out
    locked_at = "VALIDATION_RECOMPUTE"
    if lock_out.is_file():
        locked_at = json.loads(lock_out.read_text(encoding="utf-8")).get("locked_utc", locked_at)
    lock = {
        "schema_version": "atlas-v1.0-molecular-search.1", "locked_utc": locked_at,
        "results_accessed_before_lock": False, "locus_count": len(loci), "dataset_count_per_locus": len(datasets),
        "search_task_count": len(rows), "search_task_ids_in_locked_order": task_ids,
        "plan_sha256": plan_sha, "policy_sha256": sha256(policy_path), "sources_sha256": sha256(sources_path),
        "preflight_sha256": sha256(preflight_path), "loci_sha256": sha256(loci_path),
        "variants_sha256": sha256(variants_path), "fine_mapping_provenance_sha256": sha256(provenance_path),
        "chain": chain_provenance,
        "metadata_sha256": {asset["id"]: sha256(root / asset["path"]) for asset in policy["metadata_assets"]},
        "claim_limit": policy["claim_limit"],
    }
    if not args.validate_only:
        if out.exists() or lock_out.exists():
            fail("molecular search plan/lock already exists; refusing overwrite")
        lock["locked_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        atomically_write(out, plan_text)
        lock["plan_sha256"] = sha256(out)
        atomically_write(lock_out, json.dumps(lock, indent=2, sort_keys=True) + "\n")
    else:
        if not out.is_file() or out.read_text(encoding="utf-8") != plan_text:
            fail("molecular search plan differs from deterministic recomputation")
        observed_lock = json.loads(lock_out.read_text(encoding="utf-8"))
        lock["locked_utc"] = observed_lock.get("locked_utc")
        if observed_lock != lock:
            fail("molecular search plan lock differs from deterministic recomputation")
    if not args.quiet:
        print(f"MOLECULAR_SEARCH_PLAN_OK loci={len(loci)} datasets={len(datasets)} tasks={len(rows)} result_free=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
