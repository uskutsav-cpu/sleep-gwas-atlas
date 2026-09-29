#!/usr/bin/env python3
"""Make a finite, conservative decision for each of the 19 partial PLACO regions."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/loci"
LOCI = BASE / "placo_factor_normalized_sensitivity_v2/placo_candidate_loci.tsv"
PHASE3 = BASE / "independent_1000g_phase3_eur_ld_v1/region_summary.tsv"
HIGHCOV = BASE / "independent_1000g_highcov_targeted_ld_v1/targeted_highcov_ld.tsv"
RAW = BASE / "independent_1000g_ld_range_exceptions_v1/provenance.json"
CATALOG = BASE / "gwas_catalog_candidate_leads_v1/trait_relevant_context.tsv"
STABILITY = BASE / "placo_candidate_locus_lead_stability_v1.tsv"
SIGNED = BASE / "placo_candidate_signed_effect_ld_partial_v1.tsv"
TIER = ROOT / "brain6/config/shared_locus_evidence_tiers_v1.json"
OUT = BASE / "candidate_region_decisions_v1"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    loci = rows(LOCI)
    phase3 = {x["locus_id"]: x for x in rows(PHASE3)}
    highcov = rows(HIGHCOV)
    catalog = rows(CATALOG)
    stability = {(x["pair_id"], x["lead_variant"]): x for x in rows(STABILITY)}
    directions = defaultdict(list)
    for x in rows(SIGNED):
        if x["locus_id"] != "NA":
            directions[x["locus_id"]].append(x)
    tier = json.loads(TIER.read_text())
    if len(loci) != 19 or tier["status"] != "FROZEN_PRE_ANNOTATION":
        raise ValueError("Unexpected locus count or tier lock state")
    if tier["prerequisites"]["placo"].find("five-track") < 0:
        raise ValueError("Frozen five-track prerequisite was altered")

    output = []
    for locus in loci:
        lid, pair = locus["locus_id"], locus["pair_id"]
        check = phase3[lid]
        lead_ids = locus["lead_variants"].split(";")
        if not all((pair, snp) in stability and stability[(pair, snp)]["stability_class"] == "STABLE_ALL_THREE"
                   for snp in lead_ids):
            raise ValueError(f"Unexpected unstable candidate lead: {lid}")
        local_highcov = [x for x in highcov if x["locus_id"] == lid]
        highcov_matches = sum(int(x["highcov_exact_alleles"]) for x in local_highcov)
        highcov_leads = sum(int(x["highcov_exact_alleles"]) for x in local_highcov if x["candidate_status"] == "LEAD")
        catalog_rows = [x for x in catalog if x["locus_id"] == lid]
        relation = Counter(x["direction_relation"] for x in directions[lid])
        ld_note = (f"Phase3 EUR exact alleles {check['n_exact_vcf']}/{check['n_candidates']} candidates, "
                   f"{check['n_leads_exact_vcf']}/{check['n_leads']} leads; "
                   f"{check['n_edges_measured']} lead edges measured")
        if local_highcov:
            ld_note += f"; 30x EUR exact alleles {highcov_matches}/{len(local_highcov)} candidates, {highcov_leads}/{check['n_leads']} leads"
        special = []
        if pair == "insomnia__mdd" and locus["CHR"] == "5":
            special.append("Phase3 lead absent, recovered by unique hg19-to-hg38 liftOver in 30x genotypes; all 36 lead LD memberships remain above r2=0.1")
        if pair == "longsleep__parkinson" and locus["CHR"] == "17":
            special.append("two Parkinson chr17 region leads have EUR r2=0.650, so independent signals are unproven")
        if pair == "longsleep__scz" and locus["CHR"] == "17":
            special.append("three schizophrenia chr17 leads have pairwise EUR r2=0.825-0.973; the sole Phase3 lead-edge cutoff conflict is resolved by 30x r2=0.995")
        if locus["CHR"] == "15":
            special.append("rs1051168 is the same lead in bipolar and schizophrenia pair regions, not separate replication")
        if relation["OPPOSING"]:
            special.append("signed sleep and disorder effects oppose at the lead; this is compatible with shared association")
        if pair == "insomnia__mdd" and locus["CHR"] == "7":
            special.append("GWAS Catalog has both insomnia and MDD lead associations, but source overlap and independent bivariate replication are unresolved")
        if pair == "longsleep__scz" and locus["CHR"] == "2":
            special.append("GWAS Catalog sleep hit is continuous sleep duration, not the locked long-sleep tail")
        reason = ("Terminal-valid protected insomnia–ADHD Track B is unavailable, so the frozen complete five-track PLACO prerequisite cannot be satisfied; "
                  "canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED and no eligible independent exact-pair replication is established. "
                  + " ".join(special))
        output.append({"Region": f"chr{locus['CHR']}:{locus['START']}-{locus['STOP']}",
                       "Pair": pair, "Decision": "BLOCKED_EXTERNAL", "Reason": reason,
                       "Fine-map": "NOT_RUN_GATE", "Coloc": "NOT_RUN_GATE", "Gene": "NOT_PRIORITIZED",
                       "Tissue": "NOT_PRIORITIZED", "Cell type": "NOT_PRIORITIZED",
                       "Therapeutic target": "NONE_DEFENSIBLE",
                       "locus_id": lid, "lead_variants": locus["lead_variants"],
                       "PLACO_lead_P": locus["lead_P_PLACO"], "LD_check": ld_note,
                       "stable_partial_leads": "YES", "signed_effect_relation": ";".join(f"{k}:{v}" for k, v in sorted(relation.items())),
                       "catalog_trait_context_rows": len(catalog_rows),
                       "external_replication_adjudication": "NOT_ESTABLISHED"})

    if len(output) != 19 or Counter(x["Decision"] for x in output) != {"BLOCKED_EXTERNAL": 19}:
        raise ValueError("Incomplete finite decision table")
    OUT.mkdir(parents=True, exist_ok=False)
    fields = list(output[0])
    with (OUT / "decisions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)

    header = ["# Brain6: finite decision for 19 current PLACO candidate regions", "",
              "The decision overlay does not change the frozen evidence tier or canonical LAVA v3. "
              "All 19 partial-family regions are `BLOCKED_EXTERNAL`; none is promoted or rejected as an association. "
              "This is a negative promotion result, not evidence that these trait pairs have no shared genetic architecture.", "",
              "| Region | Pair | Decision | Reason | Fine-map | Coloc | Gene | Tissue | Cell type | Therapeutic target |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for x in output:
        detail = x["Reason"].replace("Terminal-valid protected insomnia–ADHD Track B is unavailable, so the frozen complete five-track PLACO prerequisite cannot be satisfied; canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED and no eligible independent exact-pair replication is established. ", "")
        concise = "Missing five-track PLACO; LAVA v3 failed QC; no eligible exact-pair replication."
        if detail:
            concise += " " + detail
        header.append("| " + " | ".join((x["Region"], x["Pair"], x["Decision"], concise,
                                       x["Fine-map"], x["Coloc"], x["Gene"], x["Tissue"],
                                       x["Cell type"], x["Therapeutic target"])) + " |")
    footer = ["", "## New evidence and limits", "",
              "- The official 1000 Genomes Phase 3 GRCh37 EUR panel (503 participants) gives exact SNP allele matches for 2,262/2,297 candidate rows and 20/21 leads. It measures 2,210 candidate-to-lead edges. One Phase 3 edge crosses the r²=0.1 threshold relative to LAVA-derived LD.",
              "- An official 30× GRCh38 call set using the same 503 EUR sample IDs and a uniquely mapped UCSC chain recovers all 36 chr5 candidate SNPs, including the missing lead rs40465. Its 36/36 lead LD memberships agree with the LAVA-derived r² threshold. It matches 1,824/1,832 schizophrenia chr17 candidate SNPs. For rs55938136 versus rs62054437, its r²=0.995 agrees with LAVA-derived r²=0.997, whereas Phase 3 reports r²=0.002. The Phase 3 call is discordant at this SNP and should not disqualify the full region.",
              "- Independent Phase 3 EUR genotypes measure 12,151 of 12,321 listed raw LAVA LD-range-exception rows; all have r²≥0.922. This supports very high LD at those variant pairs but does not make the original >1 factor products valid correlations or repair the failed canonical LAVA family.",
              "- The two Parkinson chr17 region leads have EUR r²=0.650. Three schizophrenia chr17 leads have pairwise r²=0.825–0.973. Separate causal signals or locus independence cannot be asserted from the partial clumping alone.",
              "- Current GWAS Catalog v2 yields 228 curated associations across the 20 unique lead rsIDs and 24 broadly tagged relevant rows. Catalog top-hit coverage is incomplete; the relevant hits do not document a non-overlapping, exact-phenotype bivariate replication. Insomnia/MDD rs2894699 has hits for both traits, but Watanabe's insomnia meta-analysis includes UK Biobank, which is also used in the locked insomnia discovery. A long-sleep/SCZ lead hit for sleep duration is a continuous trait and does not reproduce the locked ≥9 h long-sleep phenotype.",
              "", "## Exact unblockers", "",
              "1. A terminal-valid, source-complete protected insomnia–ADHD Track B package with aligned input/source/lock identities, followed by a complete five-track PLACO family correction under the frozen rule. This is required before any current region can receive even exploratory TIER_3.",
              "2. For TIER_1 or TIER_2 local sharing, a separately predeclared, complete, QC-passing local-rg family with defensible LD and the required robustness gates. Frozen canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED.",
              "3. For TIER_1 independent replication, exact sleep and disorder phenotypes in a source-verified non-overlapping cohort with direction and P-value checks. The current GWAS Catalog context cannot substitute.",
              "", "No fine-mapping, colocalization, QTL, regulatory, gene, tissue, cell-type, mechanism, or therapeutic claim is scientifically eligible under the frozen annotation gate. No intervention shortlist is emitted.", ""]
    (OUT / "decision_report.md").write_text("\n".join(header + footer))
    receipt = {"analysis_id": "brain6_candidate_region_decisions_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "status_counts": dict(Counter(x["Decision"] for x in output)),
               "immutability": "Decision overlay only; canonical LAVA v3 and prior PLACO outputs untouched",
               "inputs_sha256": {str(x): sha256(x) for x in (LOCI, PHASE3, HIGHCOV, RAW, CATALOG, STABILITY, SIGNED, TIER, Path(__file__))},
               "outputs_sha256": {name: sha256(OUT / name) for name in ("decisions.tsv", "decision_report.md")}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt["status_counts"]))


if __name__ == "__main__":
    main()
