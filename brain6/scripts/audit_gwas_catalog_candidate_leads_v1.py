#!/usr/bin/env python3
"""Collect current GWAS Catalog v2 association context for 21 PLACO leads.

Curated top associations are contextual only: they do not prove independent
replication, bivariate sharing, or non-overlap with the discovery cohorts.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv"
OUT = ROOT / "brain6/results/loci/gwas_catalog_candidate_leads_v1"
API = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
KEYWORDS = {
    "sleep": ("sleep", "insomnia"),
    "mdd": ("depress", "major depressive"),
    "bipolar": ("bipolar",),
    "parkinson": ("parkinson",),
    "scz": ("schizophren",),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Brain6/1.0 research-audit"})
    with urllib.request.urlopen(request, timeout=45) as response:
        return json.load(response)


def main() -> None:
    with SOURCE.open(newline="") as handle:
        leads = [x for x in csv.DictReader(handle, delimiter="\t") if x["candidate_status"] == "LEAD"]
    if len(leads) != 21:
        raise ValueError(f"Expected 21 leads, got {len(leads)}")
    by_snp: dict[str, list[dict]] = {}
    for row in leads:
        by_snp.setdefault(row["SNP"], []).append(row)

    records = []
    receipts = []
    for snp in sorted(by_snp):
        url = API + "?" + urllib.parse.urlencode({"rs_id": snp, "size": 500})
        n_pages = 0
        while url:
            page = fetch_json(url)
            n_pages += 1
            for item in page.get("_embedded", {}).get("associations", []):
                trait = "; ".join(x.get("efo_trait", "") for x in item.get("efo_traits", []))
                reported = "; ".join(item.get("reported_trait", []))
                records.append({"lead_snp": snp, "association_id": item.get("association_id", ""),
                                "study_accession": item.get("accession_id", ""),
                                "pubmed_id": item.get("pubmed_id", ""),
                                "first_author": item.get("first_author", ""),
                                "efo_trait": trait, "reported_trait": reported,
                                "p_value": item.get("p_value", ""),
                                "beta": item.get("beta", ""),
                                "risk_frequency": item.get("risk_frequency", ""),
                                "snp_effect_allele": "; ".join(item.get("snp_effect_allele", [])),
                                "catalog_locations": "; ".join(item.get("locations", []))})
            url = page.get("_links", {}).get("next", {}).get("href")
            if n_pages > 100:
                raise RuntimeError(f"Unexpected pagination for {snp}")
            time.sleep(0.08)
        receipts.append({"snp": snp, "pages": n_pages,
                         "associations": sum(x["lead_snp"] == snp for x in records),
                         "query": API + "?" + urllib.parse.urlencode({"rs_id": snp, "size": 500})})
        print(f"{snp}: {receipts[-1]['associations']} curated associations", flush=True)

    relevant = []
    for lead in leads:
        disorder = lead["pair_id"].split("__")[1]
        for item in records:
            if item["lead_snp"] != lead["SNP"]:
                continue
            haystack = (item["efo_trait"] + " " + item["reported_trait"]).lower()
            tags = [label for label, words in (("sleep", KEYWORDS["sleep"]),
                                                (disorder, KEYWORDS[disorder]))
                    if any(word in haystack for word in words)]
            if tags:
                relevant.append({"locus_id": lead["locus_id"], "pair_id": lead["pair_id"],
                                 "evidence_trait_class": ";".join(tags), **item})

    OUT.mkdir(parents=True, exist_ok=False)
    def write(name: str, values: list[dict], fields: list[str]) -> None:
        with (OUT / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(values)
    write("catalog_associations.tsv", records, list(records[0]))
    write("trait_relevant_context.tsv", relevant,
          ["locus_id", "pair_id", "evidence_trait_class"] + list(records[0]))
    receipt = {"analysis_id": "brain6_gwas_catalog_candidate_leads_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "api_documentation": "https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/",
               "limitations": ["Curated top associations only; absence is not evidence of no association",
                               "No study sample overlap or independent replication established",
                               "Trait tags use broad string matching and require scientific review"],
               "input_sha256": {str(SOURCE): sha256(SOURCE), str(Path(__file__)): sha256(Path(__file__))},
               "queries": receipts,
               "output_sha256": {name: sha256(OUT / name) for name in
                                 ("catalog_associations.tsv", "trait_relevant_context.tsv")}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"{len(records)} total associations, {len(relevant)} tagged trait-context rows", flush=True)


if __name__ == "__main__":
    main()
