#!/usr/bin/env python3
"""Attach source-bound lead eQTL context to the 20-region exploratory ledger."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
V3 = ROOT / "brain6/results/exploratory_five_track_v3"
QTL = ROOT / "brain6/results/exploratory_five_track_v4/lead_eqtl_context"
OUT = QTL.parent


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check_receipt(folder: Path) -> None:
    receipt = json.loads((folder / "provenance.json").read_text())
    for name, digest in receipt["output_sha256"].items():
        if sha(folder / name) != digest:
            raise ValueError(f"Output hash mismatch: {folder / name}")
    for name, digest in receipt.get("query_sha256", {}).items():
        if sha(folder / name) != digest:
            raise ValueError(f"Raw QTL response hash mismatch: {folder / name}")


def write(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if (OUT / "geographic_region_evidence.tsv").exists():
        raise FileExistsError("Versioned v4 region table already exists")
    check_receipt(V3)
    check_receipt(QTL)
    previous = read(V3 / "geographic_region_evidence.tsv")
    summary = read(QTL / "lead_eqtl_query_summary.tsv")
    associations = read(QTL / "lead_eqtl_associations.tsv")
    if len(previous) != 20 or len(summary) != 54 or len(associations) != 826:
        raise ValueError("20-region/27-lead/two-dataset coverage changed")
    by_region_summary = defaultdict(list)
    by_region_associations = defaultdict(list)
    for row in summary:
        by_region_summary[row["region_group"]].append(row)
    for row in associations:
        by_region_associations[row["region_group"]].append(row)
    if set(by_region_summary) != {r["geographic_region"] for r in previous}:
        raise ValueError("QTL summary failed to cover every geographic region")
    updated = []
    for prior in previous:
        row = dict(prior)
        region = row["geographic_region"]
        summaries = by_region_summary[region]
        context = by_region_associations[region]
        distinct_queries = {(q["lead_snp"], q["chr37"], q["bp37"], q["dataset_id"]) for q in summaries}
        distinct_hits = {(a["lead_snp"], a["dataset_id"], a["molecular_trait_id"], a["qtl_variant"]) for a in context}
        best = min(context, key=lambda a: float(a["nominal_pvalue"])) if context else None
        row["gtex_brain_lead_eqtl_unique_queries"] = str(len(distinct_queries))
        row["gtex_brain_lead_eqtl_queries_with_exact_rows"] = str(len({
            (q["lead_snp"], q["chr37"], q["bp37"], q["dataset_id"]) for q in summaries
            if q["status"] == "EXACT_LEAD_QTL_ROWS_PRESENT"
        }))
        row["gtex_brain_lead_eqtl_distinct_nominal_rows"] = str(len(distinct_hits))
        row["gtex_brain_lead_eqtl_min_nominal_p"] = best["nominal_pvalue"] if best else "NOT_OBSERVED_AT_EXACT_LEADS"
        row["gtex_brain_lead_eqtl_min_p_gene_id"] = best["gene_id"] if best else "NOT_OBSERVED_AT_EXACT_LEADS"
        row["gtex_brain_lead_eqtl_interpretation"] = "EXPLORATORY_LEAD_LEVEL_NOMINAL_ONLY"
        updated.append(row)
    if sum(int(r["gtex_brain_lead_eqtl_unique_queries"]) for r in updated) != 52:
        raise ValueError("Distinct source queries failed to reconcile")
    table = OUT / "geographic_region_evidence.tsv"
    write(table, updated)
    lines = [
        "# Five-track PLACO exploratory 20-region evidence — lead eQTL context added",
        "",
        "**EXPLORATORY. Zero LAVA-confirmed regions and zero final shared-locus tiers.**",
        "",
        "All 25 pair-specific candidates, 27 lead rows, and 20 coordinate-union groups are retained.",
        "Exact lead GRCh37 coordinates were lifted uniquely to GRCh38 and queried against two",
        "prespecified GTEx v8 brain gene-expression datasets (cortex and DLPFC). All 52 distinct",
        "query responses are cached with hashes. The table reports nominal context only; it does",
        "not assign an eQTL-supported gene or a shared causal signal.",
        "",
        "| Region | Pair candidates | Leads | Exact lead eQTL queries / queried | Distinct nominal rows | Lowest nominal p | Gene ID at lowest p |",
        "|---|---:|---|---:|---:|---:|---|",
    ]
    for r in updated:
        lines.append(
            f"| {r['geographic_region']} | {r['n_pair_specific_candidates']} | {r['lead_variants']} | "
            f"{r['gtex_brain_lead_eqtl_queries_with_exact_rows']}/{r['gtex_brain_lead_eqtl_unique_queries']} | "
            f"{r['gtex_brain_lead_eqtl_distinct_nominal_rows']} | "
            f"{r['gtex_brain_lead_eqtl_min_nominal_p']} | {r['gtex_brain_lead_eqtl_min_p_gene_id']} |"
        )
    lines.extend([
        "",
        "The lowest p values are unadjusted minima over varying numbers of genes, lead variants,",
        "and two tissues. They cannot rank regions or justify gene/tissue claims. An absent exact",
        "lead row does not imply that no QTL exists elsewhere in that region. The full v3 evidence",
        "columns remain in the TSV, including replication scope, signed LD, independent LD, and",
        "explicit unrun fine-mapping, coloc, splicing QTL, regulatory, cell-type, and pathway states.",
        "No biologically supported shared region is formally validated or ranked. LAVA local-rg",
        "and final evidence-tier conclusions remain BLOCKED_LAVA under the unchanged canonical",
        "FAILED_QC_NOT_PROMOTED decision.", "",
    ])
    report = OUT / "region_evidence_report.md"
    report.write_text("\n".join(lines))
    receipt = {
        "analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA", "source_v3_receipt_sha256": sha(V3 / "provenance.json"),
        "source_eqtl_receipt_sha256": sha(QTL / "provenance.json"),
        "output_sha256": {p.name: sha(p) for p in (table, report)},
    }
    with (OUT / "provenance.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"regions": len(updated), "unique_queries": 52,
                      "regions_with_exact_lead_qtl_rows": sum(int(r["gtex_brain_lead_eqtl_distinct_nominal_rows"]) > 0 for r in updated)}, sort_keys=True))


if __name__ == "__main__":
    main()
