#!/usr/bin/env python3
"""Execute one locked eQTL Catalogue/GTEx interval query and normalize to GRCh37."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from liftover_chain import load_chain, reverse_complement


SUPPORTED_FAMILIES = {
    "EQTL_CATALOGUE_R7_GE", "EQTL_CATALOGUE_R7_LEAFCUTTER",
    "EQTL_CATALOGUE_R7_SUN2018", "EQTL_CATALOGUE_GTEX_V8_IMPORTED",
}
BASES = {"A", "C", "G", "T"}


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


def gzip_text(path: Path, text: str) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def finite(value: str, field: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field}") from exc
    if not math.isfinite(number):
        raise ValueError(f"nonfinite {field}")
    return number


def orient_to_locked_alleles(
    effect: str, other: str, beta: float, locked_effect: str, locked_other: str,
) -> tuple[str, str, float]:
    if (effect, other) == (locked_effect, locked_other):
        return locked_effect, locked_other, beta
    if (effect, other) == (locked_other, locked_effect):
        return locked_effect, locked_other, -beta
    raise ValueError("source and locked allele orientations disagree")


def header_snapshot(url: str) -> dict[str, object]:
    result = subprocess.run(["curl", "-fsSI", url], capture_output=True, text=True)
    if result.returncode:
        return {"url": url, "status": "HEAD_FAILED", "stderr": result.stderr.strip()}
    headers: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            headers[key.strip().lower()] = value.strip()
    return {"url": url, "status": "HEAD_OK", "headers": headers}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("search_task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--plan", default="results/tables/molecular_search_plan.tsv")
    parser.add_argument("--plan-lock", default="results/tables/molecular_search_plan.lock.json")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--variants", default="results/atlas/variants.tsv")
    parser.add_argument("--out-dir", default="results/molecular/search")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        fail("source access requires explicit --execute after reviewing the locked search plan")
    root = Path(args.root).resolve()
    plan_path, lock_path = root / args.plan, root / args.plan_lock
    policy_path, variants_path = root / args.policy, root / args.variants
    plan_fields, plan = read_tsv(plan_path)
    del plan_fields
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("plan_sha256") != sha256(plan_path) or lock.get("results_accessed_before_lock") is not False:
        fail("molecular search plan differs from its result-free lock")
    selected = [row for row in plan if row["search_task_id"] == args.search_task_id]
    if len(selected) != 1:
        fail("search_task_id must identify exactly one locked query")
    task = selected[0]
    if task["source_family_id"] not in SUPPORTED_FAMILIES or task["query_mode"] != "TABIX_GRCH38_INTERVAL":
        fail("this executor only accepts locked eQTL Catalogue/GTEx tabix tasks")
    if task["results_accessed_before_lock"] != "NO" or task["planned_outcome"] != "SOURCE_QUERY_REQUIRED":
        fail("query task is not a result-free source-access row")
    tabix = shutil.which("tabix")
    if not tabix:
        fail("tabix is not installed")
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    output_dir = root / args.out_dir / args.search_task_id
    if output_dir.exists():
        fail(f"query output already exists; refusing overwrite: {output_dir}")
    _, variants = read_tsv(variants_path)
    locus_variants = [row for row in variants if row["locus_id"] == task["locus_id"]]
    if not locus_variants:
        fail("locked query locus has no fine-mapped variant universe")
    reference: dict[tuple[int, int, frozenset[str]], dict[str, str]] = {}
    for row in locus_variants:
        key = (int(row["chromosome"]), int(row["position_bp"]), frozenset((row["effect_allele"], row["other_allele"])))
        if key in reference and reference[key]["variant_id"] != row["variant_id"]:
            fail(f"ambiguous fine-mapped coordinate/allele key at {task['locus_id']}")
        reference[key] = row
    chain_spec = policy["source_to_analysis_chain"]
    chain, chain_provenance = load_chain(root / chain_spec["path"], chain_spec["sha256"], chain_spec["bytes"])
    region = f"{task['chromosome_grch38']}:{task['start_grch38']}-{task['end_grch38']}"
    work_root = root / "work/molecular"
    work_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=args.search_task_id + ".", dir=work_root))
    try:
        query = subprocess.run([tabix, "-h", task["source_url"], region], capture_output=True, text=True)
        if query.returncode:
            fail(f"tabix query failed for {args.search_task_id}: {query.stderr.strip()}")
        raw_text = query.stdout
        raw_path = staging / "source_query.tsv.gz"
        gzip_text(raw_path, raw_text)
        lines = [line for line in raw_text.splitlines() if line and not line.startswith("##")]
        header_line = next((line for line in lines if line.startswith("#")), None)
        if header_line is not None:
            header = header_line.lstrip("#").split("\t")
            data_lines = lines[lines.index(header_line) + 1:]
        elif lines:
            header = lines[0].split("\t")
            data_lines = lines[1:]
        else:
            header, data_lines = [], []
        required = {"chromosome", "position", "ref", "alt", "beta", "se", "molecular_trait_id"}
        if data_lines and not required.issubset(header):
            fail(f"source response lacks required eQTL Catalogue fields: {sorted(required - set(header))}")
        source_rows = list(csv.DictReader(data_lines, fieldnames=header, delimiter="\t")) if data_lines else []
        excluded: Counter[str] = Counter()
        normalized: dict[tuple[str, str], dict[str, object]] = {}
        for row in source_rows:
            try:
                source_chr, source_pos = int(row["chromosome"].removeprefix("chr")), int(row["position"])
                ref, alt = row["ref"].upper(), row["alt"].upper()
                beta, se = finite(row["beta"], "beta"), finite(row["se"], "se")
            except (KeyError, ValueError):
                excluded["MALFORMED_CORE_FIELDS"] += 1; continue
            if source_chr > 22 or ref not in BASES or alt not in BASES or ref == alt:
                excluded["NON_AUTOSOMAL_OR_NON_BIALLELIC_SNP"] += 1; continue
            status, mapped = chain.map_point(source_chr, source_pos)
            if status != "mapped" or mapped is None:
                excluded[f"LIFTOVER_{status.upper()}"] += 1; continue
            chromosome, position, strand = mapped
            effect, other = (alt, ref) if strand == "+" else (reverse_complement(alt), reverse_complement(ref))
            if effect + other in {"AT", "TA", "CG", "GC"}:
                excluded["AMBIGUOUS_STRAND"] += 1; continue
            key = (chromosome, position, frozenset((effect, other)))
            reference_row = reference.get(key)
            if reference_row is None:
                excluded["NOT_IN_LOCKED_LD_VARIANT_UNIVERSE"] += 1; continue
            variant_id = reference_row["variant_id"]
            locked_effect, locked_other = reference_row["effect_allele"].upper(), reference_row["other_allele"].upper()
            try:
                locked_effect, locked_other, beta = orient_to_locked_alleles(
                    effect, other, beta, locked_effect, locked_other,
                )
            except ValueError:
                fail(f"allele-set lookup produced an orientation mismatch: {variant_id}")
            try:
                if row.get("maf") not in {None, "", "NA", "."}:
                    maf = finite(row["maf"], "maf")
                else:
                    ac, an = finite(row.get("ac", ""), "ac"), finite(row.get("an", ""), "an")
                    maf = min(ac / an, 1 - ac / an)
                if row.get("an") not in {None, "", "NA", "."}:
                    n = finite(row["an"], "an") / 2
                else:
                    n = finite(task["sample_size"], "sample_size")
            except ValueError:
                excluded["INVALID_MAF_OR_N"] += 1; continue
            if not 0 < maf <= 0.5 or maf < policy["variant_qc"]["minimum_maf"] or se <= 0 or n <= 0:
                excluded["MAF_SE_OR_N_QC"] += 1; continue
            feature_id = row["molecular_trait_id"]
            gene_id = row.get("gene_id", "NA") or "NA"
            record = {
                "source_variant_id": row.get("variant", f"chr{source_chr}_{source_pos}_{ref}_{alt}"),
                "rsid": variant_id,
                "source_chromosome": source_chr, "source_position": source_pos,
                "source_ref": ref, "source_alt": alt, "chromosome": chromosome,
                "position": position, "effect_allele": locked_effect, "other_allele": locked_other,
                "beta": format(beta, ".15g"), "se": format(se, ".15g"),
                "maf": format(maf, ".15g"), "n": format(n, ".15g"),
                "feature_id": feature_id, "feature_name": row.get("molecular_trait_object_id", feature_id) or feature_id,
                "gene_id": gene_id, "gene_symbol": "NA", "context": task["context"],
                "source_dataset_id": task["dataset_id"], "modality": task["modality"],
                "liftover_status": "MAPPED_FORWARD" if strand == "+" else "MAPPED_REVERSE",
            }
            identity = (feature_id, variant_id)
            if identity in normalized and normalized[identity] != record:
                fail(f"disagreeing duplicate source rows for {feature_id}/{variant_id}")
            normalized[identity] = record
        normalized_rows = sorted(normalized.values(), key=lambda row: (str(row["feature_id"]), int(row["position"]), str(row["rsid"])))
        normalized_path = staging / "normalized_qtl.tsv.gz"
        gzip_text(normalized_path, table_text(policy["normalized_qtl_schema"], normalized_rows))
        outcome = "DATA_READY" if normalized_rows else "NO_FEATURES_IN_LOCUS" if not source_rows else "NO_ANALYZABLE_FEATURES"
        remote = {
            "source": header_snapshot(task["source_url"]),
            "index": header_snapshot(task["source_index_url"]),
            "queried_region": region,
        }
        headers_path = staging / "remote_headers.json"
        headers_path.write_text(json.dumps(remote, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        provenance = {
            "schema_version": "atlas-v1.0-molecular-query.1",
            "queried_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "search_task_id": args.search_task_id, "search_plan_sha256": sha256(plan_path),
            "search_plan_lock_sha256": sha256(lock_path), "policy_sha256": sha256(policy_path),
            "source_family_id": task["source_family_id"], "study_id": task["study_id"],
            "dataset_id": task["dataset_id"], "exact_release": task["exact_release"],
            "source_url": task["source_url"], "source_index_url": task["source_index_url"],
            "query_region": region, "search_outcome": outcome, "raw_source_row_count": len(source_rows),
            "normalized_row_count": len(normalized_rows),
            "normalized_feature_count": len({row["feature_id"] for row in normalized_rows}),
            "exclusion_counts": dict(sorted(excluded.items())), "chain": chain_provenance,
            "outputs": {
                "source_query.tsv.gz": {"bytes": raw_path.stat().st_size, "sha256": sha256(raw_path)},
                "normalized_qtl.tsv.gz": {"bytes": normalized_path.stat().st_size, "sha256": sha256(normalized_path)},
                "remote_headers.json": {"bytes": headers_path.stat().st_size, "sha256": sha256(headers_path)},
            },
            "results_accessed_before_lock": False, "claim_limit": policy["claim_limit"],
        }
        provenance_path = staging / "provenance.json"
        provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staging, output_dir)
        print(f"MOLECULAR_QUERY_OK task={args.search_task_id} outcome={outcome} rows={len(normalized_rows)} features={provenance['normalized_feature_count']}")
        return 0
    finally:
        if staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
