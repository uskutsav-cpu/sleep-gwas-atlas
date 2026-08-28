#!/usr/bin/env python3
"""Map one full GRCh37 GWAS to an exact locked GTEx-v8 model-variant family."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import tempfile
from collections import Counter
from pathlib import Path

from liftover_chain import load_chain, reverse_complement


OUTPUT_FIELDS = ["SNP", "A1", "A2", "BETA", "SE"]
AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def gzip_text(path: Path, text: str) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(text.encode("utf-8"))


def table_text(rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trait_id")
    parser.add_argument("model_family")
    parser.add_argument("--root", default=".")
    parser.add_argument("--mapping-manifest", default="results/tables/twas_gwas_mapping_manifest.tsv")
    parser.add_argument("--run-lock", default="results/tables/twas_run_manifest.lock.json")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--materialize", action="store_true")
    args = parser.parse_args()
    if not args.materialize:
        fail("TWAS GWAS mapping requires explicit --materialize")
    root = Path(args.root).resolve()
    mapping_path, run_lock_path, policy_path = root / args.mapping_manifest, root / args.run_lock, root / args.policy
    _, mappings = read_tsv(mapping_path)
    run_lock = json.loads(run_lock_path.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if run_lock.get("mapping_sha256") != sha256(mapping_path) or run_lock.get("policy_sha256") != sha256(policy_path):
        fail("TWAS mapping manifest differs from its lock")
    selected = [row for row in mappings if row["trait_id"] == args.trait_id and row["model_family"] == args.model_family]
    if len(selected) != 1:
        fail("trait/model family must identify exactly one mapping task")
    task = selected[0]
    source, variant_path = root / task["full_gwas_path"], root / task["model_variant_path"]
    if sha256(source) != task["full_gwas_sha256"] or sha256(variant_path) != task["model_variant_sha256"]:
        fail("TWAS mapping source differs from lock")
    variant_fields, variants = read_tsv(variant_path)
    if variant_fields != ["model_variant_id", "chromosome_grch38", "position_grch38", "ref", "alt"]:
        fail("TWAS model-variant registry schema drifted")
    model_by_key: dict[tuple[int, int, frozenset[str]], list[str]] = {}
    for row in variants:
        key = (int(row["chromosome_grch38"]), int(row["position_grch38"]), frozenset((row["ref"], row["alt"])))
        model_by_key.setdefault(key, []).append(row["model_variant_id"])
    chain_spec = policy["reference_build"]
    chain, chain_provenance = load_chain(root / chain_spec["chain_path"], chain_spec["chain_sha256"], chain_spec["chain_bytes"])
    output: dict[str, dict[str, object]] = {}
    excluded: Counter[str] = Counter()
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE"}
        if not required.issubset(reader.fieldnames or []):
            fail("full GWAS lacks canonical fields for TWAS mapping")
        for row in reader:
            try:
                chromosome, position = int(row["CHR"]), int(row["BP"])
                a1, a2 = row["A1"].upper(), row["A2"].upper()
                beta, se = float(row["BETA"]), float(row["SE"])
            except ValueError:
                excluded["MALFORMED"] += 1; continue
            if not math.isfinite(beta) or not math.isfinite(se) or se <= 0 or frozenset((a1, a2)) in AMBIGUOUS:
                excluded["STATISTIC_OR_AMBIGUOUS_ALLELE_QC"] += 1; continue
            status, mapped = chain.map_point(chromosome, position)
            if status != "mapped" or mapped is None:
                excluded[f"LIFTOVER_{status.upper()}"] += 1; continue
            mapped_chr, mapped_position, strand = mapped
            effect, other = (a1, a2) if strand == "+" else (reverse_complement(a1), reverse_complement(a2))
            model_ids = model_by_key.get((mapped_chr, mapped_position, frozenset((effect, other))))
            if not model_ids:
                excluded["NOT_IN_MODEL_FAMILY"] += 1; continue
            for model_id in model_ids:
                record = {"SNP": model_id, "A1": effect, "A2": other, "BETA": format(beta, ".15g"), "SE": format(se, ".15g")}
                if model_id in output and output[model_id] != record:
                    fail(f"multiple disagreeing GWAS rows map to model variant: {model_id}")
                output[model_id] = record
    if not output:
        fail("TWAS mapping retained no model variants")
    ordered = [output[row["model_variant_id"]] for row in variants if row["model_variant_id"] in output]
    final, final_lock = root / task["mapped_gwas_path"], root / task["mapped_gwas_lock_path"]
    if final.exists() or final_lock.exists():
        fail("mapped TWAS GWAS already exists; refusing overwrite")
    final.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=final.name + ".", dir=final.parent)
    os.close(descriptor); temporary = Path(temporary_name)
    try:
        gzip_text(temporary, table_text(ordered))
        os.replace(temporary, final)
    finally:
        temporary.unlink(missing_ok=True)
    lock_payload = {
        "schema_version": "atlas-v1.0-twas-gwas-map.1", "mapping_id": task["mapping_id"],
        "mapped_variant_count": len(ordered), "output_sha256": sha256(final),
        "source_gwas_sha256": task["full_gwas_sha256"], "model_variant_sha256": task["model_variant_sha256"],
        "mapping_manifest_sha256": sha256(mapping_path), "run_manifest_lock_sha256": sha256(run_lock_path),
        "policy_sha256": sha256(policy_path), "chain": chain_provenance,
        "exclusion_counts": dict(sorted(excluded.items())), "results_accessed": False,
        "claim_limit": policy["claim_limit"],
    }
    final_lock.write_text(json.dumps(lock_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"TWAS_GWAS_MAPPING_OK trait={args.trait_id} family={args.model_family} variants={len(ordered)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
