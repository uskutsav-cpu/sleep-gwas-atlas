#!/usr/bin/env python3
"""Stream the pinned ABC 2021 atlas once and cache eligible same-locus links."""
from __future__ import annotations

import argparse
from bisect import bisect_right
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path


MISSING = {"", "NA"}
ABC_FIELDS = [
    "chr", "start", "end", "name", "class", "activity_base", "TargetGene",
    "TargetGeneTSS", "TargetGeneExpression", "TargetGenePromoterActivityQuantile",
    "TargetGeneIsExpressed", "distance", "isSelfPromoter", "hic_contact",
    "powerlaw_contact", "powerlaw_contact_reference", "hic_contact_pl_scaled",
    "hic_pseudocount", "hic_contact_pl_scaled_adj", "ABC.Score.Numerator",
    "ABC.Score", "powerlaw.Score.Numerator", "powerlaw.Score", "CellType",
]


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


def eligible_variant(row: dict[str, str]) -> bool:
    return row["qc_status"] == "PASS" and any(
        row[field] not in MISSING
        for field in ("credible_set_sleep", "credible_set_non_sleep", "shared_signal_posterior")
    )


def fmt(value: float) -> str:
    return format(value, ".12g")


def evidence_id(*values: str) -> str:
    return "AUTOABC::" + hashlib.sha256("\x1f".join(values).encode()).hexdigest()[:24]


def abc_components(root: Path, policy: dict[str, object]) -> tuple[Path, dict[str, object], Path]:
    spec = policy["abc_2021"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("ABC component manifest is absent or differs from the policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components = manifest.get("components", [])
    if len(components) != 1 or components[0].get("component_id") != "ABC_ALL_PREDICTIONS":
        fail("ABC component manifest is not the exact required one-file bundle")
    component = components[0]
    source_path = root / component["path"]
    if (
        not source_path.is_file() or source_path.stat().st_size != component["bytes"]
        or sha256(source_path) != component["sha256"]
    ):
        fail("ABC all-predictions source differs from its exact pin")
    return source_path, manifest, manifest_path


def prepare_rows(
    policy: dict[str, object], manifest: dict[str, object], source_path: Path,
    variants: list[dict[str, str]], genes: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, object]]:
    selected_variants = [row for row in variants if eligible_variant(row)]
    positions: dict[int, list[tuple[int, dict[str, str]]]] = {}
    for row in selected_variants:
        try:
            chromosome, position = int(row["chromosome"].removeprefix("chr")), int(row["position_bp"])
        except (KeyError, ValueError) as exc:
            fail(f"invalid eligible atlas variant coordinate: {row.get('variant_id', 'UNKNOWN')}")
            raise AssertionError from exc
        if not 1 <= chromosome <= 22 or position < 1:
            fail(f"eligible atlas variant is outside autosomal GRCh37: {row['variant_id']}")
        positions.setdefault(chromosome, []).append((position, row))
    position_only: dict[int, list[int]] = {}
    for chromosome, values in positions.items():
        values.sort(key=lambda item: (item[0], item[1]["locus_id"], item[1]["variant_id"]))
        position_only[chromosome] = [value[0] for value in values]

    genes_by_locus_symbol: dict[tuple[str, str], list[str]] = {}
    for row in genes:
        symbol = row.get("gene_symbol", "NA")
        if symbol in MISSING:
            continue
        genes_by_locus_symbol.setdefault((row["locus_id"], symbol), []).append(row["gene_id"])

    domain_map = manifest["domain_biosamples"]
    biosample_domain = {
        biosample: domain for domain, biosamples in domain_map.items() for biosample in biosamples
    }
    if len(biosample_domain) != sum(len(values) for values in domain_map.values()):
        fail("ABC domain map assigns a biosample more than once")
    released_threshold = float(policy["abc_2021"]["released_score_threshold"])
    primary_threshold = float(policy["abc_2021"]["primary_score_threshold"])
    counters = {
        "source_rows": 0, "mapped_domain_rows": 0, "self_promoter_rows_excluded": 0,
        "interval_variant_overlaps": 0, "unsupported_target_overlaps_excluded": 0,
        "ambiguous_target_overlaps_excluded": 0,
    }
    source_biosamples: set[str] = set()
    rows: list[dict[str, str]] = []
    with gzip.open(source_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ABC_FIELDS:
            fail("ABC source header differs from the frozen 24-column schema")
        for source in reader:
            counters["source_rows"] += 1
            biosample = source["CellType"]
            source_biosamples.add(biosample)
            domain = biosample_domain.get(biosample)
            if domain is None:
                continue
            counters["mapped_domain_rows"] += 1
            if policy["abc_2021"]["exclude_self_promoters"] and source["isSelfPromoter"] == "True":
                counters["self_promoter_rows_excluded"] += 1
                continue
            try:
                chromosome = int(source["chr"].removeprefix("chr"))
                start, end = int(source["start"]), int(source["end"])
                score = float(source["ABC.Score"])
            except ValueError:
                continue
            if chromosome not in positions or start < 0 or end <= start:
                continue
            if not math.isfinite(score) or score < released_threshold - 1e-12 or score > 1:
                fail("ABC source contains an out-of-policy score")
            left = bisect_right(position_only[chromosome], start)
            right = bisect_right(position_only[chromosome], end)
            for _position, variant in positions[chromosome][left:right]:
                counters["interval_variant_overlaps"] += 1
                candidates = sorted(set(genes_by_locus_symbol.get(
                    (variant["locus_id"], source["TargetGene"]), []
                )))
                if not candidates:
                    counters["unsupported_target_overlaps_excluded"] += 1
                    continue
                if len(candidates) != 1:
                    counters["ambiguous_target_overlaps_excluded"] += 1
                    continue
                gene_id = candidates[0]
                element = f"{source['chr']}:{start}-{end}:{source['name']}"
                tier = "PRIMARY_GE_0.02" if score >= primary_threshold else "RELEASE_GE_0.015"
                rows.append({
                    "regulatory_evidence_id": evidence_id(
                        variant["locus_id"], variant["variant_id"], element,
                        biosample, gene_id, fmt(score),
                    ),
                    "regulatory_element_id": element, "variant_id": variant["variant_id"],
                    "locus_id": variant["locus_id"], "element_type": "enhancer_promoter_link",
                    "annotation": (
                        f"ABC_CLASS::{source['class']}::TIER::{tier}::DISTANCE_BP::{source['distance']}"
                    ),
                    "effect": fmt(score), "p_value": "NA", "biosample": biosample,
                    "tissue": "NA", "cell_type": biosample, "context_domain": domain,
                    "target_gene_id": gene_id, "link_method": "ACTIVITY_BY_CONTACT",
                    "source_dataset": policy["abc_2021"]["source_id"],
                    "source_version": manifest["release"],
                    "evidence_level": f"ABC_{tier}", "provenance_id": "ABC_CACHE_PENDING",
                })
    if counters["source_rows"] != manifest["row_count"] or len(source_biosamples) != manifest["biosample_count"]:
        fail("ABC streamed source shape differs from the release manifest")
    rows.sort(key=lambda row: (
        row["locus_id"], row["context_domain"], row["variant_id"],
        row["regulatory_element_id"], row["biosample"], row["target_gene_id"],
    ))
    identities = [row["regulatory_evidence_id"] for row in rows]
    if len(identities) != len(set(identities)):
        fail("ABC cache contains duplicate regulatory evidence identities")
    counters.update({
        "source_biosamples": len(source_biosamples),
        "eligible_variant_rows": len(selected_variants),
        "eligible_variant_keys": len({(row["locus_id"], row["variant_id"]) for row in selected_variants}),
        "supported_gene_rows": len(genes), "result_rows": len(rows),
    })
    return rows, counters


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--variants", default="results/atlas/variants.tsv")
    parser.add_argument("--genes", default="results/atlas/genes.tsv")
    parser.add_argument("--out")
    parser.add_argument("--provenance-out")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    source_path, manifest, manifest_path = abc_components(root, policy)
    _, registry = read_tsv(root / policy["source_registry"])
    source = next((row for row in registry if row["source_id"] == policy["abc_2021"]["source_id"]), None)
    if (
        source is None or source["source_status"] != "SOURCE_VERIFIED"
        or source["exact_release"] != manifest["release"]
        or root / source["local_path"] != source_path
    ):
        fail("ABC source registry is absent, unverified, or release/path mismatched")
    variants_path, genes_path = root / args.variants, root / args.genes
    _, variants = read_tsv(variants_path)
    _, genes = read_tsv(genes_path)
    rows, counters = prepare_rows(policy, manifest, source_path, variants, genes)
    fields = policy["regulatory_mapping"]["canonical_fields"]
    if any(list(row) != fields for row in rows):
        fail("ABC cache row differs from the locked regulatory schema")
    payload = table_text(fields, rows)
    out_path = root / (args.out or policy["abc_2021"]["cache_path"])
    provenance_path = root / (args.provenance_out or policy["abc_2021"]["cache_provenance_path"])
    atomic_text(out_path, payload)
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "source_id": source["source_id"], "source_release": source["exact_release"],
        "policy_sha256": sha256(policy_path), "component_manifest_sha256": sha256(manifest_path),
        "source_path": str(source_path), "source_bytes": source_path.stat().st_size,
        "source_sha256": sha256(source_path), "variants_sha256": sha256(variants_path),
        "genes_sha256": sha256(genes_path), "cache_path": str(out_path),
        "cache_sha256": hashlib.sha256(payload.encode()).hexdigest(), "counts": counters,
        "domain_biosamples": manifest["domain_biosamples"],
        "claim_limit": policy["abc_2021"]["claim_limit"],
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(
        f"ABC_OVERLAP_CACHE_OK source_rows={counters['source_rows']} "
        f"eligible_variants={counters['eligible_variant_rows']} rows={len(rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
