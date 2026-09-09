#!/usr/bin/env python3
"""Build ranked and domain-specific tables from the isolated extension family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rg", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
    )
    parser.add_argument(
        "--positive-out", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_ranked_positive.tsv"),
    )
    parser.add_argument(
        "--negative-out", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_ranked_negative.tsv"),
    )
    parser.add_argument(
        "--domain-pairs-out", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_domain_pair_views.tsv"),
    )
    parser.add_argument(
        "--domain-summary-out", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_domain_summary.tsv"),
    )
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/extension_views.json"),
    )
    args = parser.parse_args()
    rows = read_tsv(args.rg)
    required = {
        "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
        "rg", "se", "p", "extension_fdr", "analysis_status",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("ERROR: extension rg matrix is empty or lacks publication fields")
    pairs = [(row["sleep_trait"], row["extension_trait_id"]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise SystemExit("ERROR: duplicate pair in extension rg matrix")
    if any(row["analysis_status"] != "PRIMARY_EXTENSION_RG_COMPLETE" for row in rows):
        raise SystemExit("ERROR: ranked views may contain only completed primary extension pairs")

    base_fields = list(rows[0])
    positive_source = sorted(
        (row for row in rows if float(row["rg"]) > 0),
        key=lambda row: (-float(row["rg"]), float(row["extension_fdr"]), row["sleep_trait"], row["extension_trait_id"]),
    )
    negative_source = sorted(
        (row for row in rows if float(row["rg"]) < 0),
        key=lambda row: (float(row["rg"]), float(row["extension_fdr"]), row["sleep_trait"], row["extension_trait_id"]),
    )
    positive = [{"direction_rank": index, **row} for index, row in enumerate(positive_source, 1)]
    negative = [{"direction_rank": index, **row} for index, row in enumerate(negative_source, 1)]
    ranked_fields = ["direction_rank", *base_fields]
    write_tsv(args.positive_out, ranked_fields, positive)
    write_tsv(args.negative_out, ranked_fields, negative)

    by_domain: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_domain[row["phenotype_domain"]].append(row)
    domain_pairs: list[dict[str, object]] = []
    domain_summary: list[dict[str, object]] = []
    for domain in sorted(by_domain):
        domain_rows = by_domain[domain]
        ranked = sorted(
            domain_rows,
            key=lambda row: (float(row["extension_fdr"]), -abs(float(row["rg"])), row["sleep_trait"], row["extension_trait_id"]),
        )
        domain_pairs.extend({"within_domain_rank": index, **row} for index, row in enumerate(ranked, 1))
        rg_values = [float(row["rg"]) for row in domain_rows]
        domain_summary.append({
            "phenotype_domain": domain,
            "extension_trait_count": len({row["extension_trait_id"] for row in domain_rows}),
            "primary_pair_count": len(domain_rows),
            "extension_fdr_significant_pair_count": sum(float(row["extension_fdr"]) < 0.05 for row in domain_rows),
            "fdr_and_abs_rg_ge_0.15_pair_count": sum(float(row["extension_fdr"]) < 0.05 and abs(float(row["rg"])) >= 0.15 for row in domain_rows),
            "minimum_extension_fdr": min(float(row["extension_fdr"]) for row in domain_rows),
            "maximum_positive_rg": max(rg_values),
            "minimum_negative_rg": min(rg_values),
            "median_rg": statistics.median(rg_values),
            "analysis_status": "DOMAIN_VIEW_COMPLETE",
        })
    write_tsv(args.domain_pairs_out, ["within_domain_rank", *base_fields], domain_pairs)
    write_tsv(args.domain_summary_out, list(domain_summary[0]), domain_summary)

    provenance = {
        "schema_version": "1.0.0", "source_rg": str(args.rg),
        "source_rg_sha256": sha256(args.rg), "source_pair_count": len(rows),
        "positive_pair_count": len(positive), "negative_pair_count": len(negative),
        "zero_rg_pair_count": len(rows) - len(positive) - len(negative),
        "domain_count": len(domain_summary),
        "sorting": {
            "positive": "rg descending, then extension FDR ascending",
            "negative": "rg ascending, then extension FDR ascending",
            "domain": "extension FDR ascending, then abs(rg) descending",
        },
        "outputs": {
            str(args.positive_out): sha256(args.positive_out),
            str(args.negative_out): sha256(args.negative_out),
            str(args.domain_pairs_out): sha256(args.domain_pairs_out),
            str(args.domain_summary_out): sha256(args.domain_summary_out),
        },
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(
        f"EXTENSION_VIEWS_OK pairs={len(rows)} positive={len(positive)} "
        f"negative={len(negative)} domains={len(domain_summary)}"
    )


if __name__ == "__main__":
    main()
