#!/usr/bin/env python3
"""Assemble the complete five-track exploratory region evidence ledger."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
V5 = ROOT / "brain6/results/exploratory_five_track_v5"
REG = ROOT / "brain6/results/exploratory_five_track_v6/lead_regulatory_context"
OUT = REG.parent


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check(folder: Path) -> None:
    receipt = json.loads((folder / "provenance.json").read_text())
    for name, digest in receipt["output_sha256"].items():
        if sha(folder / name) != digest:
            raise ValueError(f"Receipt hash mismatch: {folder / name}")
    for name, digest in receipt.get("raw_query_sha256", {}).items():
        if sha(folder / name) != digest:
            raise ValueError(f"Raw regulatory response hash mismatch: {folder / name}")


def write(path: Path, records: list[dict[str, str]]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    if (OUT / "provenance.json").exists():
        raise FileExistsError("Versioned v6 evidence report already exists")
    check(V5)
    check(V5 / "gtex_brain_qtl_credible_sets")
    check(REG)
    previous = read(V5 / "geographic_region_evidence.tsv")
    annotations = read(REG / "lead_regulatory_summary.tsv")
    features = read(REG / "lead_regulatory_features.tsv")
    if len(previous) != 20 or len(annotations) != 27 or len(features) != 6:
        raise ValueError("Twenty-region/27-lead regulatory accounting changed")
    by_region_annotations = defaultdict(list)
    by_region_features = defaultdict(list)
    for item in annotations:
        by_region_annotations[item["region_group"]].append(item)
    for item in features:
        by_region_features[item["region_group"]].append(item)
    if set(by_region_annotations) != {r["geographic_region"] for r in previous}:
        raise ValueError("Regulatory summary does not cover all regions")
    updated = []
    for prior in previous:
        row = dict(prior)
        region = row["geographic_region"]
        leads = by_region_annotations[region]
        observations = by_region_features[region]
        row["ensembl_exact_lead_regulatory_feature_count"] = str(len({
            (x["lead_snp"], x["feature_id"]) for x in observations
        }))
        row["ensembl_exact_leads_with_regulatory_overlap"] = str(len({x["lead_snp"] for x in observations}))
        row["ensembl_exact_lead_regulatory_feature_ids"] = ";".join(sorted({x["feature_id"] for x in observations})) or "NONE_AT_EXACT_LEADS"
        row["ensembl_exact_lead_regulatory_interpretation"] = "EXPLORATORY_COORDINATE_OVERLAP_ONLY"
        if len({x["lead_snp"] for x in leads}) != len(set(row["lead_variants"].split(";"))):
            raise ValueError(f"Region lead identity changed: {region}")
        updated.append(row)
    if sum(int(r["n_pair_specific_candidates"]) for r in updated) != 25:
        raise ValueError("25 pair-specific candidates not preserved")
    if any(r["lava_local_rg"] != "BLOCKED_LAVA" or r["final_region_tier"] != "BLOCKED_LAVA" for r in updated):
        raise ValueError("Protected LAVA/tier gate changed")
    table = OUT / "geographic_region_evidence.tsv"
    write(table, updated)

    lines = [
        "# Brain6 five-track exploratory evidence across 20 geographic regions",
        "",
        "**EXPLORATORY. Canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED; no final shared-locus tier is assigned.**",
        "",
        "The source-bound table retains all 25 pair-specific PLACO candidate loci, 20 coordinate",
        "groups, 2,686 candidate rows, 27 lead rows, and all prior evidence columns. Pair-level",
        "FinnGen directional replication applies only to insomnia–ADHD globally; it is not",
        "locus-specific replication for any of its six rows below. The GTEx columns distinguish",
        "nominal lead eQTL associations from exact PLACO-candidate overlap with published",
        "molecular credible sets. Neither is GWAS fine-mapping or colocalization.", "",
        "| Geographic region | Pair(s) | Lead(s) | Valid signed LD / candidates | Global pair replication | Exact lead eQTL queries with rows | Candidates in molecular QTL sets | Highest molecular eQTL / sQTL PIP | Exact lead regulatory features | Prior trait-relevant Catalog rows |",
        "|---|---|---|---:|---|---:|---:|---|---:|---:|",
    ]
    for r in updated:
        eqtl_pip = r["gtex_brain_eqtl_cs_max_molecular_pip"]
        sqtl_pip = r["gtex_brain_sqtl_cs_max_molecular_pip"]
        eqtl_pip = "none" if eqtl_pip == "NO_EXACT_CANDIDATE_CS_OVERLAP" else eqtl_pip
        sqtl_pip = "none" if sqtl_pip == "NO_EXACT_CANDIDATE_CS_OVERLAP" else sqtl_pip
        lines.append(
            f"| {r['geographic_region']} | {r['pairs']} | {r['lead_variants']} | "
            f"{r['n_valid_signed_ld_rows']}/{r['n_candidate_variant_rows']} | "
            f"{r['pair_level_global_replication']} | "
            f"{r['gtex_brain_lead_eqtl_queries_with_exact_rows']}/{r['gtex_brain_lead_eqtl_unique_queries']} | "
            f"{r['gtex_brain_qtl_cs_candidate_variants']} | {eqtl_pip} / {sqtl_pip} | "
            f"{r['ensembl_exact_lead_regulatory_feature_count']} | {r['prior_catalog_trait_relevant_records']} |"
        )
    lines.extend([
        "",
        "The separate TSV carries positional genes, source-specific Catalog studies, five",
        "independent-genotype cross-pair LD edges, all signed-LD QC exceptions, and every",
        "method state. Six candidate rows have out-of-range raw signed correlations and 24",
        "have no exact frozen-reference allele; all 27 lead rows pass signed LD. In the",
        "regulatory lookup, six feature rows overlap five exact leads; these are coordinate",
        "annotations only. An empty exact-lead query cannot exclude a regulatory feature",
        "or QTL elsewhere in the region.", "",
        "## Current interpretation", "",
        "The chr11:61.1–62.1 Mb long-sleep–bipolar region has five candidate variants in",
        "one published brain eQTL credible set, a molecular maximum PIP of 0.199, and an",
        "enhancer coordinate overlap at its lead. The chr17:43.0–45.3 Mb long-sleep–PD/SCZ",
        "group has extensive candidate overlap with published brain expression and splicing",
        "credible sets; it also contains 1,845 PLACO candidate rows, so raw overlap count",
        "is particularly sensitive to candidate density. These are descriptive molecular",
        "hypotheses, not ranked or validated shared loci. The chr15 group also has candidate",
        "eQTL-set overlap, with low molecular PIP in this panel. Nearby gene labels and",
        "Catalog rows remain positional/prior-study context; they do not establish causality.", "",
        "**No region has locus-specific independent replication, validated GWAS SuSiE",
        "fine-mapping, trait–trait coloc, trait–eQTL/sQTL coloc, tissue/cell-type",
        "enrichment, or pathway support.** The available nominal eQTL and molecular",
        "credible-set observations do not fill those missing method states. LAVA local-rg,",
        "final evidence tiers, and confirmatory cross-layer claims remain BLOCKED_LAVA.", "",
    ])
    report = OUT / "region_evidence_report.md"
    report.write_text("\n".join(lines))
    receipt = {"analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA",
               "source_v5_receipt_sha256": sha(V5 / "provenance.json"),
               "source_v6_regulatory_receipt_sha256": sha(REG / "provenance.json"),
               "output_sha256": {p.name: sha(p) for p in (table, report)}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"regions": len(updated), "pair_candidates": 25,
                      "regulatory_feature_rows": len(features),
                      "regions_with_molecular_cs_overlap": sum(
                          int(r["gtex_brain_qtl_cs_candidate_variants"]) > 0 for r in updated)}, sort_keys=True))


if __name__ == "__main__":
    main()
