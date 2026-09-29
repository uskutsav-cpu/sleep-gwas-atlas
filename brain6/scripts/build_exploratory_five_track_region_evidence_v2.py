#!/usr/bin/env python3
"""Add cached Track B lead context to the unchanged 20-region v1 ledger."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / "brain6/results/exploratory_five_track_v1"
B = V1 / "b_annotation"
OUT = ROOT / "brain6/results/exploratory_five_track_v2"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, data: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned exploratory v2 output already exists")
    v1_receipt = json.loads((V1 / "geographic_region_evidence.provenance.json").read_text())
    b_receipt = json.loads((B / "provenance.json").read_text())
    for name, digest in v1_receipt["output_sha256"].items():
        if sha(V1 / name) != digest:
            raise ValueError(f"v1 output hash mismatch: {name}")
    for name, digest in b_receipt["output_sha256"].items():
        if sha(B / name) != digest:
            raise ValueError(f"Track B annotation hash mismatch: {name}")
    old = read(V1 / "geographic_region_evidence.tsv")
    b_leads = read(B / "lead_annotation_summary.tsv")
    b_catalog = read(B / "catalog_exact_rsid_associations.tsv")
    if len(old) != 20 or len(b_leads) != 6:
        raise ValueError("20 geographic groups or six B leads missing")
    studies = defaultdict(set)
    for record in b_catalog:
        studies[record["locus_id"]].add(record["study_accession"])
    updated = []
    encountered = set()
    for record in old:
        ids = set(record["candidate_locus_ids"].split(";"))
        matches = [x for x in b_leads if x["locus_id"] in ids]
        encountered.update(x["locus_id"] for x in matches)
        entry = dict(record)
        entry["b_nearest_gene_positional_context"] = ";".join(
            sorted({x["nearest_gene_symbol"] for x in matches})) or "NOT_APPLICABLE"
        entry["b_catalog_exact_rsid_records_unfiltered"] = str(sum(
            int(x["catalog_exact_rsid_associations"]) for x in matches))
        entry["b_catalog_studies_unfiltered"] = ";".join(sorted({
            study for x in matches for study in studies[x["locus_id"]]
        })) or "NONE_IN_OFFICIAL_SNAPSHOT"
        entry["b_ensembl_regulatory_overlaps_at_lead"] = str(sum(
            int(x["ensembl_regulatory_overlaps_at_lead"]) for x in matches))
        entry["b_annotation_interpretation"] = (
            "EXPLORATORY_COORDINATE_AND_CATALOG_CONTEXT_ONLY" if matches else "NOT_APPLICABLE"
        )
        updated.append(entry)
    if encountered != {x["locus_id"] for x in b_leads}:
        raise ValueError("B lead not assigned to a geographic region")
    OUT.mkdir(parents=True, exist_ok=False)
    target = OUT / "geographic_region_evidence.tsv"
    write(target, updated)
    report = OUT / "region_evidence_report.md"
    lines = [
        "# Five-track PLACO exploratory 20-region evidence, with admitted B annotation",
        "",
        "**EXPLORATORY. No region is LAVA-confirmed or formally tiered.**",
        "",
        "The detailed TSV retains all 25 pair-specific candidates within 20",
        "geographic groups. It incorporates the official-source, cached",
        "GRCh37 positional annotation of six admitted insomnia-ADHD leads.",
        "GWAS Catalog counts for B are unfiltered exact-rsID records; they are",
        "separate from the older four-pair trait-relevant snapshot and are not",
        "independent replication. Ensembl regulatory overlap at a lead is",
        "coordinate context, not functional validation.",
        "",
        "| Geographic region | Pairs | Leads | Candidate rows | Existing positional genes | B positional gene | B catalog records |",
        "|---|---|---|---:|---|---|---:|",
    ]
    for x in updated:
        lines.append(
            f"| {x['geographic_region']} | {x['pairs']} | {x['lead_variants']} | "
            f"{x['n_candidate_variant_rows']} | {x['nearest_gene_context']} | "
            f"{x['b_nearest_gene_positional_context']} | {x['b_catalog_exact_rsid_records_unfiltered']} |"
        )
    lines.extend([
        "",
        "The only available replication is **global pair-level** FinnGen",
        "directional evidence for insomnia-ADHD. No individual region has",
        "locus-specific replication. No region has validated fine-mapping,",
        "trait-trait or QTL colocalization, regulatory, tissue/cell-type, or",
        "pathway support. The strongest biological support category is therefore",
        "**none established**. Four multi-pair geographic overlaps on",
        "chromosomes 5, 7, 15, and 17 remain descriptive hypotheses, with no",
        "formal ordering or tier promotion.",
        "",
        "Canonical LAVA v3 remains `FAILED_QC_NOT_PROMOTED`: 3,720 of 17,465",
        "cells were NOT_RUN versus the frozen maximum of 873. Local-rg claims,",
        "confirmatory region tiers, and dependent cross-layer conclusions stay",
        "`BLOCKED_LAVA`. The 24 reference-unmatched PLACO candidate rows remain",
        "preserved in the v1 unassigned table.",
        "",
    ])
    report.write_text("\n".join(lines))
    provenance = {
        "analysis_id": "brain6_exploratory_five_track_region_evidence_v2",
        "analysis_label": "EXPLORATORY",
        "geographic_regions": 20,
        "pair_specific_candidates": 25,
        "source_sha256": {
            str((V1 / "geographic_region_evidence.tsv").relative_to(ROOT)): sha(V1 / "geographic_region_evidence.tsv"),
            str((V1 / "geographic_region_evidence.provenance.json").relative_to(ROOT)): sha(V1 / "geographic_region_evidence.provenance.json"),
            str((B / "lead_annotation_summary.tsv").relative_to(ROOT)): sha(B / "lead_annotation_summary.tsv"),
            str((B / "catalog_exact_rsid_associations.tsv").relative_to(ROOT)): sha(B / "catalog_exact_rsid_associations.tsv"),
            str((B / "provenance.json").relative_to(ROOT)): sha(B / "provenance.json"),
        },
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {target.name: sha(target), report.name: sha(report)},
        "limitations": [
            "All gene and regulatory records are descriptive coordinates, not molecular support.",
            "GWAS Catalog records are not independent replication and may reuse discovery cohorts.",
            "Frozen LAVA and candidate promotion gates remain unchanged.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"regions": len(updated), "B_leads": len(encountered),
                      "B_catalog_records": len(b_catalog)}, sort_keys=True))


if __name__ == "__main__":
    main()
