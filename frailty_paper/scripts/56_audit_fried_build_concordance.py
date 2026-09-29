#!/usr/bin/env python3
"""Compare Fried FFS coordinates with the pinned GRCh37 HapMap3 map."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable, Iterator

Variant = tuple[int, int, str, str, str]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def group_by_coordinate(rows: Iterable[Variant], sort_errors: list[int]) -> Iterator[tuple[tuple[int, int], list[Variant]]]:
    iterator = iter(rows)
    first = next(iterator, None)
    if first is None:
        return
    coordinate = first[:2]
    group = [first]
    previous = coordinate
    for row in iterator:
        key = row[:2]
        if key < previous:
            sort_errors[0] += 1
        previous = key
        if key == coordinate:
            group.append(row)
        else:
            yield coordinate, group
            coordinate, group = key, [row]
    yield coordinate, group


def merge_concordance(gwas_rows: Iterable[Variant], reference_rows: Iterable[Variant]) -> dict[str, int]:
    """Streaming coordinate-group intersection preserving multi-allelic positions."""
    counts = {
        "summary_variants_compared": 0,
        "coordinate_intersection": 0,
        "same_coordinate_same_rsid": 0,
        "unmatched_variant_rows_at_shared_coordinates": 0,
        "same_rsid_same_allele_order": 0,
        "same_rsid_swapped_alleles": 0,
        "same_rsid_allele_mismatch": 0,
        "summary_sort_violations": 0,
        "reference_sort_violations": 0,
    }
    summary_sort = [0]
    reference_sort = [0]
    a, b = iter(group_by_coordinate(gwas_rows, summary_sort)), iter(group_by_coordinate(reference_rows, reference_sort))
    x, y = next(a, None), next(b, None)
    while x is not None and y is not None:
        if x[0] < y[0]:
            counts["summary_variants_compared"] += len(x[1])
            x = next(a, None)
        elif y[0] < x[0]:
            y = next(b, None)
        else:
            counts["coordinate_intersection"] += 1
            counts["summary_variants_compared"] += len(x[1])
            ref_by_id: dict[str, list[Variant]] = {}
            for variant in y[1]:
                ref_by_id.setdefault(variant[2], []).append(variant)
            matched_ids = 0
            for variant in x[1]:
                candidates = ref_by_id.get(variant[2], [])
                if not candidates:
                    continue
                reference = candidates.pop()
                matched_ids += 1
                counts["same_coordinate_same_rsid"] += 1
                if variant[3:] == reference[3:]:
                    counts["same_rsid_same_allele_order"] += 1
                elif variant[3:] == (reference[4], reference[3]):
                    counts["same_rsid_swapped_alleles"] += 1
                else:
                    counts["same_rsid_allele_mismatch"] += 1
            counts["unmatched_variant_rows_at_shared_coordinates"] += len(x[1]) + len(y[1]) - 2 * matched_ids
            x, y = next(a, None), next(b, None)
    while x is not None:
        counts["summary_variants_compared"] += len(x[1])
        x = next(a, None)
    while y is not None:
        y = next(b, None)
    counts["summary_sort_violations"] = summary_sort[0]
    counts["reference_sort_violations"] = reference_sort[0]
    return counts


def read_gwas(path: Path) -> Iterator[Variant]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"CHROM", "POS", "A1", "A2", "SNP"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"missing required summary columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            if not re.fullmatch(r"rs\d+", row["SNP"]):
                continue
            if row["A1"] not in {"A", "C", "G", "T"} or row["A2"] not in {"A", "C", "G", "T"}:
                continue
            yield int(row["CHROM"]), int(row["POS"]), row["SNP"], row["A1"], row["A2"]


def read_reference(path: Path) -> Iterator[Variant]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"CHR", "BP", "SNP", "A1", "A2"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"missing required reference columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            yield int(row["CHR"]), int(row["BP"]), row["SNP"], row["A1"], row["A2"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-provenance", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    provenance = json.loads(args.reference_provenance.read_text(encoding="utf-8"))
    if provenance.get("genome_build") != "GRCh37/hg19":
        raise ValueError("pinned coordinate map provenance does not identify GRCh37/hg19")
    counts = merge_concordance(read_gwas(args.summary), read_reference(args.reference))
    counts["reference_rows"] = int(provenance["mapped_rows"])
    result = {
        "schema_version": "frailty.fried-grch37-coordinate-concordance.v1",
        "summary_file": str(args.summary.resolve()),
        "summary_sha256": sha256(args.summary),
        "reference_file": str(args.reference.resolve()),
        "reference_sha256": sha256(args.reference),
        "reference_genome_build": provenance["genome_build"],
        "reference_provenance_sha256": sha256(args.reference_provenance),
        "matching_key": "chromosome + 1-based position grouped before comparing rsID and unordered allele pairs",
        "counts": counts,
        "interpretation": "Observed coordinate and rsID concordance supports GRCh37/hg19 coordinates for the acquired Fried FFS file; this does not validate effect-allele direction or participant overlap.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
