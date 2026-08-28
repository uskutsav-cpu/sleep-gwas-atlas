#!/usr/bin/env python3
"""Cache eligible variant-to-promoter contacts from pinned Javierre PCHi-C data."""
from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path


MISSING = {"", ".", "NA"}
BASE_FIELDS = [
    "baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart",
    "oeEnd", "oeID", "oeName", "dist",
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
    return "AUTOPCHIC::" + hashlib.sha256("\x1f".join(values).encode()).hexdigest()[:24]


def pchic_components(root: Path, policy: dict[str, object]) -> tuple[dict[str, Path], dict[str, object], Path]:
    spec = policy["pchic_2016"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("PCHi-C component manifest is absent or differs from the policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components: dict[str, Path] = {}
    for component in manifest.get("components", []):
        path = root / component["path"]
        if (
            not path.is_file() or path.stat().st_size != component["bytes"]
            or sha256(path) != component["sha256"]
        ):
            fail(f"PCHi-C component differs from its exact pin: {component['component_id']}")
        components[component["component_id"]] = path
    if set(components) != {"PCHIC_PEAK_MATRIX_CUTOFF5", "PCHIC_PEAK_MATRIX_README"}:
        fail("PCHi-C manifest is not the exact matrix-plus-README bundle")
    return components, manifest, manifest_path


def promoter_symbols(value: str) -> list[str]:
    return sorted({symbol for symbol in value.split(";") if symbol not in MISSING})


def prepare_rows(
    policy: dict[str, object], manifest: dict[str, object], matrix_path: Path,
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
        if symbol not in MISSING:
            genes_by_locus_symbol.setdefault((row["locus_id"], symbol), []).append(row["gene_id"])

    def overlapping(chromosome_text: str, start_text: str, end_text: str) -> list[dict[str, str]]:
        try:
            chromosome, start, end = int(chromosome_text), int(start_text), int(end_text)
        except ValueError:
            return []
        if chromosome not in positions or start < 1 or end < start:
            return []
        left = bisect_left(position_only[chromosome], start)
        right = bisect_right(position_only[chromosome], end)
        return [row for _position, row in positions[chromosome][left:right]]

    cells = manifest["cell_types"]
    threshold = float(policy["pchic_2016"]["chicago_score_threshold"])
    counters = {
        "source_rows": 0, "qualifying_cell_contacts": 0, "interval_variant_overlaps": 0,
        "unannotated_opposite_promoter_overlaps_excluded": 0,
        "unsupported_target_overlaps_excluded": 0, "ambiguous_target_overlaps_excluded": 0,
        "reciprocal_or_duplicate_rows_collapsed": 0,
    }
    retained: dict[tuple[str, ...], dict[str, str]] = {}
    expected_header = BASE_FIELDS + list(cells) + ["clusterID", "clusterPostProb"]
    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected_header:
            fail("PCHi-C source header differs from the frozen matrix schema")
        for source in reader:
            counters["source_rows"] += 1
            try:
                scores = {cell: float(source[cell]) for cell in cells}
            except ValueError as exc:
                fail("PCHi-C source contains a nonnumeric CHiCAGO score")
                raise AssertionError from exc
            if any(not math.isfinite(value) or value < 0 for value in scores.values()):
                fail("PCHi-C source contains an invalid CHiCAGO score")
            qualifying = [(cell, value) for cell, value in scores.items() if value >= threshold]
            if not qualifying:
                fail("PCHi-C cutoff-5 source row has no qualifying cell type")
            counters["qualifying_cell_contacts"] += len(qualifying)
            pair = "-".join(sorted((source["baitID"], source["oeID"]), key=lambda value: int(value)))
            directions = [
                (
                    overlapping(source["oeChr"], source["oeStart"], source["oeEnd"]),
                    source["baitName"], source["oeChr"], source["oeStart"], source["oeEnd"],
                    source["oeID"], "OTHER_END_TO_BAIT_PROMOTER",
                ),
                (
                    overlapping(source["baitChr"], source["baitStart"], source["baitEnd"]),
                    source["oeName"], source["baitChr"], source["baitStart"], source["baitEnd"],
                    source["baitID"], "BAIT_TO_OTHER_BAIT_PROMOTER",
                ),
            ]
            for variant_rows, target_names, fragment_chr, fragment_start, fragment_end, fragment_id, direction in directions:
                if not variant_rows:
                    continue
                counters["interval_variant_overlaps"] += len(variant_rows)
                symbols = promoter_symbols(target_names)
                if not symbols:
                    counters["unannotated_opposite_promoter_overlaps_excluded"] += len(variant_rows)
                    continue
                for variant in variant_rows:
                    for symbol in symbols:
                        candidates = sorted(set(genes_by_locus_symbol.get((variant["locus_id"], symbol), [])))
                        if not candidates:
                            counters["unsupported_target_overlaps_excluded"] += 1
                            continue
                        if len(candidates) != 1:
                            counters["ambiguous_target_overlaps_excluded"] += 1
                            continue
                        gene_id = candidates[0]
                        element = f"chr{fragment_chr}:{fragment_start}-{fragment_end}:RF{fragment_id}"
                        for cell, score in qualifying:
                            key = (
                                variant["locus_id"], variant["variant_id"], element,
                                gene_id, cell,
                            )
                            row = {
                                "regulatory_evidence_id": evidence_id(*key),
                                "regulatory_element_id": element, "variant_id": variant["variant_id"],
                                "locus_id": variant["locus_id"],
                                "element_type": "3D_contact_where_available",
                                "annotation": f"CHICAGO_FRAGMENT_PAIR::{pair}::DIRECTION::{direction}",
                                "effect": fmt(score), "p_value": "NA", "biosample": cell,
                                "tissue": "blood", "cell_type": cells[cell],
                                "context_domain": "immune", "target_gene_id": gene_id,
                                "link_method": "PROMOTER_CAPTURE_HIC_CHICAGO_GE5",
                                "source_dataset": policy["pchic_2016"]["source_id"],
                                "source_version": manifest["release"],
                                "evidence_level": "PHYSICAL_PROMOTER_CONTACT_CHICAGO_GE5",
                                "provenance_id": "PCHIC_CACHE_PENDING",
                            }
                            previous = retained.get(key)
                            if previous is None or float(previous["effect"]) < score:
                                if previous is not None:
                                    counters["reciprocal_or_duplicate_rows_collapsed"] += 1
                                retained[key] = row
                            else:
                                counters["reciprocal_or_duplicate_rows_collapsed"] += 1
    if counters["source_rows"] != manifest["row_count"]:
        fail("PCHi-C streamed row count differs from the source manifest")
    rows = sorted(retained.values(), key=lambda row: (
        row["locus_id"], row["variant_id"], row["regulatory_element_id"],
        row["target_gene_id"], row["biosample"], row["regulatory_evidence_id"],
    ))
    counters.update({
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
    components, manifest, manifest_path = pchic_components(root, policy)
    _, registry = read_tsv(root / policy["source_registry"])
    source = next((row for row in registry if row["source_id"] == policy["pchic_2016"]["source_id"]), None)
    matrix_path = components["PCHIC_PEAK_MATRIX_CUTOFF5"]
    if (
        source is None or source["source_status"] != "SOURCE_VERIFIED"
        or source["exact_release"] != manifest["release"]
        or root / source["local_path"] != matrix_path
    ):
        fail("PCHi-C source registry is absent, unverified, or release/path mismatched")
    variants_path, genes_path = root / args.variants, root / args.genes
    _, variants = read_tsv(variants_path)
    _, genes = read_tsv(genes_path)
    rows, counters = prepare_rows(policy, manifest, matrix_path, variants, genes)
    fields = policy["regulatory_mapping"]["canonical_fields"]
    if any(list(row) != fields for row in rows):
        fail("PCHi-C cache row differs from the locked regulatory schema")
    payload = table_text(fields, rows)
    out_path = root / (args.out or policy["pchic_2016"]["cache_path"])
    provenance_path = root / (args.provenance_out or policy["pchic_2016"]["cache_provenance_path"])
    atomic_text(out_path, payload)
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "source_id": source["source_id"], "source_release": source["exact_release"],
        "policy_sha256": sha256(policy_path), "component_manifest_sha256": sha256(manifest_path),
        "source_path": str(matrix_path), "source_bytes": matrix_path.stat().st_size,
        "source_sha256": sha256(matrix_path), "readme_sha256": sha256(components["PCHIC_PEAK_MATRIX_README"]),
        "variants_sha256": sha256(variants_path), "genes_sha256": sha256(genes_path),
        "cache_path": str(out_path), "cache_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "counts": counters, "cell_types": manifest["cell_types"],
        "claim_limit": policy["pchic_2016"]["claim_limit"],
    }
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(
        f"PCHIC_OVERLAP_CACHE_OK source_rows={counters['source_rows']} "
        f"eligible_variants={counters['eligible_variant_rows']} rows={len(rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
