#!/usr/bin/env python3
"""Cache official descriptive annotation for the six admitted Track B leads."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1"
OUT = ROOT / "brain6/results/exploratory_five_track_v1/b_annotation"
ENSEMBL = "https://grch37.rest.ensembl.org"
CATALOG = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_rows(path: Path, data: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def get_json(url: str) -> object:
    result = subprocess.run(
        ["curl", "--fail", "--silent", "--show-error", "--location", "--retry", "3",
         "--retry-delay", "2", "--max-time", "30", "--header", "Accept: application/json",
         "--header", "Content-type: application/json", url],
        check=True, capture_output=True, text=True,
    )
    return json.loads(result.stdout)


def distance(gene: dict, pos: int) -> tuple[int, int, str]:
    start, end = int(gene["start"]), int(gene["end"])
    body = max(start - pos, pos - end, 0)
    tss = start if int(gene["strand"]) == 1 else end
    return body, abs(tss - pos), str(gene.get("id", gene.get("gene_id", "")))


def main() -> None:
    if OUT.exists():
        raise FileExistsError("Versioned B annotation output exists; retain it and use a new version")
    loci = [x for x in rows(SOURCE / "candidate_loci.tsv") if x["pair_id"] == "insomnia__adhd"]
    candidate_rows = rows(SOURCE / "candidate_variants.tsv")
    leads = [x for x in candidate_rows if x["pair_id"] == "insomnia__adhd" and
             x["candidate_status"] == "LEAD"]
    if len(loci) != 6 or len(leads) != 6 or {x["locus_id"] for x in leads} != {x["locus_id"] for x in loci}:
        raise ValueError("Track B six-lead source invariant failed")
    receipt_path = ROOT / "brain6/results/placo/insomnia__adhd/admission_receipt.json"
    if sha(receipt_path) != "25c304d7107cc68ed7d8825ad4403ff031b17910d1ac3627c4751a621a985738":
        raise ValueError("Protected Track B admission receipt changed")

    raw = []
    summary = []
    catalog_rows = []
    gene_rows = []
    regulatory_rows = []
    for lead in sorted(leads, key=lambda x: (int(x["CHR"]), int(x["BP"]))):
        chrom, pos, snp = lead["CHR"], int(lead["BP"]), lead["SNP"]
        start, stop = max(1, pos - 500_000), pos + 500_000
        gene_url = f"{ENSEMBL}/overlap/region/human/{chrom}:{start}-{stop}?feature=gene"
        reg_url = f"{ENSEMBL}/overlap/region/human/{chrom}:{pos}-{pos}?feature=regulatory"
        catalog_url = f"{CATALOG}?rs_id={snp}&size=500"
        genes = get_json(gene_url)
        regulatory = get_json(reg_url)
        catalog = get_json(catalog_url)
        if not isinstance(genes, list) or not isinstance(regulatory, list) or not isinstance(catalog, dict):
            raise ValueError(f"Unexpected official API response for {snp}")
        pages = [catalog]
        association_list = list(catalog.get("_embedded", {}).get("associations", []))
        page = catalog.get("page", {})
        total_pages = int(page.get("total_pages", 1))
        if total_pages > 20:
            raise ValueError(f"Unexpectedly large catalog pagination for {snp}")
        for n in range(1, total_pages):
            next_page = get_json(f"{catalog_url}&page={n}")
            if not isinstance(next_page, dict):
                raise ValueError("Invalid GWAS Catalog page")
            pages.append(next_page)
            association_list.extend(next_page.get("_embedded", {}).get("associations", []))
        if len(association_list) != int(page.get("total_elements", len(association_list))):
            raise ValueError(f"Incomplete GWAS Catalog pagination for {snp}")
        raw.append({"snp": snp, "chromosome": chrom, "bp_grch37": pos,
                    "gene_url": gene_url, "gene_response": genes,
                    "regulatory_url": reg_url, "regulatory_response": regulatory,
                    "catalog_url": catalog_url, "catalog_pages": pages})
        nearest = min(genes, key=lambda g: distance(g, pos)) if genes else None
        if nearest:
            gene_rows.append({"locus_id": lead["locus_id"], "lead_snp": snp,
                              "gene_id": nearest.get("id", nearest.get("gene_id", "")),
                              "gene_symbol": nearest.get("external_name", ""),
                              "biotype": nearest.get("biotype", ""),
                              "distance_to_gene_body_bp": distance(nearest, pos)[0],
                              "distance_to_tss_bp": distance(nearest, pos)[1],
                              "status": "EXPLORATORY_POSITIONAL_CONTEXT_ONLY"})
        for feature in regulatory:
            regulatory_rows.append({"locus_id": lead["locus_id"], "lead_snp": snp,
                                    "feature_id": feature.get("id", ""),
                                    "feature_type": feature.get("feature_type", ""),
                                    "description": feature.get("description", ""),
                                    "start": feature.get("start", ""), "end": feature.get("end", ""),
                                    "status": "EXPLORATORY_COORDINATE_OVERLAP_ONLY"})
        for association in association_list:
            catalog_rows.append({"locus_id": lead["locus_id"], "lead_snp": snp,
                                 "association_id": association.get("association_id", ""),
                                 "study_accession": association.get("accession_id", ""),
                                 "pubmed_id": association.get("pubmed_id", ""),
                                 "first_author": association.get("first_author", ""),
                                 "reported_trait": ";".join(map(str, association.get("reported_trait", []))),
                                 "efo_trait": ";".join(str(x.get("efo_trait", "")) for x in association.get("efo_traits", [])),
                                 "p_value": association.get("p_value", ""),
                                 "status": "EXPLORATORY_CATALOG_CONTEXT_ONLY"})
        summary.append({"locus_id": lead["locus_id"], "lead_snp": snp,
                        "chr": chrom, "bp_grch37": pos,
                        "nearest_gene_symbol": nearest.get("external_name", "") if nearest else "NO_GENE_IN_500KB",
                        "nearest_gene_id": nearest.get("id", nearest.get("gene_id", "")) if nearest else "",
                        "catalog_exact_rsid_associations": len(association_list),
                        "ensembl_regulatory_overlaps_at_lead": len(regulatory),
                        "analysis_label": "EXPLORATORY",
                        "functional_support": "NOT_ESTABLISHED",
                        "lava_local_rg": "BLOCKED_LAVA"})
        print(f"{snp}: genes {len(genes)}, catalog {len(association_list)}, regulatory overlaps {len(regulatory)}", flush=True)

    OUT.mkdir(parents=True, exist_ok=False)
    raw_path = OUT / "raw_official_api_responses.json"
    raw_path.write_text(json.dumps(raw, indent=2, sort_keys=True) + "\n")
    outputs = {
        "lead_annotation_summary.tsv": (summary, list(summary[0])),
        "nearest_gene.tsv": (gene_rows, list(gene_rows[0])),
        "catalog_exact_rsid_associations.tsv": (catalog_rows,
            list(catalog_rows[0]) if catalog_rows else ["locus_id", "lead_snp", "association_id", "study_accession", "pubmed_id", "first_author", "reported_trait", "efo_trait", "p_value", "status"]),
        "regulatory_lead_overlaps.tsv": (regulatory_rows,
            list(regulatory_rows[0]) if regulatory_rows else ["locus_id", "lead_snp", "feature_id", "feature_type", "description", "start", "end", "status"]),
    }
    for name, (data, fields) in outputs.items():
        write_rows(OUT / name, data, fields)
    provenance = {
        "analysis_id": "brain6_admitted_b_exploratory_annotation_v1",
        "analysis_label": "EXPLORATORY", "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "assembly": "GRCh37", "lead_count": len(leads),
        "official_api_docs": {
            "ensembl": f"{ENSEMBL}/documentation/info/overlap_region",
            "gwas_catalog": "https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/",
        },
        "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                         (SOURCE / "candidate_loci.tsv", SOURCE / "candidate_variants.tsv",
                          SOURCE / "provenance.json", receipt_path)},
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {p.name: sha(p) for p in [raw_path] + [OUT / x for x in outputs]},
        "limitations": [
            "Nearest gene and regulatory coordinate overlap are not causal-gene or functional support.",
            "Catalog exact-rsID associations are unfiltered context; cohort overlap and phenotype equivalence are not established.",
            "No LAVA local-rg or confirmatory shared-locus tier is promoted.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
