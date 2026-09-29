#!/usr/bin/env python3
"""Independently check admitted B candidate LD in 503 1000G EUR samples."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import pysam
from pyliftover import LiftOver

from audit_1000g_phase3_eur_ld_v1 import BASE_URL, VCF_NAME, dosage, pearson_ld
from rebuild_partial_candidates_1000g_ld_v1 import HIGHCOV_URL

ROOT = Path(__file__).resolve().parents[2]
LOCI = ROOT / "brain6/results/loci/placo_five_track_candidate_loci_v1"
LEDGER = ROOT / "brain6/results/placo/insomnia__adhd/B.full.tsv.gz"
PANEL = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/integrated_call_samples_v3.20130502.ALL.panel")
CHAIN = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz")
OUT = ROOT / "brain6/results/loci/admitted_B_1000g_ld_v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def table(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream,delimiter="\t"))


def emit(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    receipt = json.loads((LOCI / "provenance.json").read_text())
    if receipt["status"] != "PASS_COMPLETE_PAIR_FAMILY_CANDIDATE_LOCI_ONLY" or receipt["candidate_variants_sha256"] != sha(LOCI / "candidate_variants.tsv"):
        raise ValueError("five-track candidate input is not validated")
    if OUT.exists():
        raise FileExistsError("B 1000G audit directory already exists")
    candidates = [x for x in table(LOCI / "candidate_variants.tsv") if x["pair_id"] == "insomnia__adhd"]
    loci = [x for x in table(LOCI / "candidate_loci.tsv") if x["pair_id"] == "insomnia__adhd"]
    if len(candidates) != 389 or len(loci) != 6 or any(x["reference_status"] != "EXACT_MATCH" for x in candidates):
        raise ValueError("admitted B candidate family differs from frozen UKB clumping")
    by_snp = {x["SNP"]:x for x in candidates}
    alleles = {}
    with gzip.open(LEDGER,"rt",newline="") as stream:
        for row in csv.DictReader(stream,delimiter="\t"):
            if row["SNP"] in by_snp:
                if row["SNP"] in alleles or (row["CHR"],row["BP"]) != (by_snp[row["SNP"]]["CHR"],by_snp[row["SNP"]]["BP"]):
                    raise ValueError("B variant identity mismatch")
                alleles[row["SNP"]] = {row["A1"],row["A2"]}
    if set(alleles) != set(by_snp):
        raise ValueError("missing B source alleles")
    samples = [x["sample"] for x in table(PANEL) if x["super_pop"] == "EUR"]
    if len(samples) != 503 or len(set(samples)) != 503:
        raise ValueError("1000G EUR sample panel changed")
    lifter = LiftOver(str(CHAIN))
    by_chr_pos: dict[str,dict[int,list[dict]]] = defaultdict(lambda:defaultdict(list))
    mapping = {}
    for row in candidates:
        chrom,bp = row["CHR"],int(row["BP"])
        if chrom == "5":
            options = lifter.convert_coordinate("chr5",bp-1)
            plus = [(c,p+1) for c,p,strand,_ in options if c=="chr5" and strand=="+"]
            if len(plus) != 1:
                mapping[row["SNP"]] = {"status":"LIFTOVER_UNRESOLVED","target_bp":""}
                continue
            target_bp = plus[0][1]
            mapping[row["SNP"]] = {"status":"HIGHCOV_GRCH38","target_bp":target_bp}
        else:
            target_bp = bp
            mapping[row["SNP"]] = {"status":"PHASE3_GRCH37","target_bp":target_bp}
        by_chr_pos[chrom][target_bp].append(row)
    genotypes = {}
    scans = []
    for chrom in sorted(by_chr_pos,key=int):
        highcov = chrom == "5"
        url = HIGHCOV_URL.format(chrom=chrom) if highcov else f"{BASE_URL}/{VCF_NAME.format(chrom=chrom)}"
        vcf = pysam.VariantFile(url)
        if set(samples) - set(vcf.header.samples):
            raise ValueError(f"EUR panel absent from chr{chrom}")
        vcf.subset_samples(samples)
        positions = by_chr_pos[chrom]
        intervals = []
        for bp in sorted(positions):
            if intervals and bp - intervals[-1][1] <= 1_000_000:
                intervals[-1][1] = bp
            else:
                intervals.append([bp,bp])
        for start,stop in intervals:
            n_records = 0
            for record in vcf.fetch(f"chr{chrom}" if highcov else chrom,start-1,stop):
                n_records += 1
                if record.pos not in positions or len(record.alts or ()) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                    continue
                observed = {record.ref,record.alts[0]}
                for row in positions[record.pos]:
                    snp = row["SNP"]
                    if observed != alleles[snp]:
                        continue
                    if snp in genotypes:
                        raise ValueError(f"duplicate exact-allele 1000G genotype for {snp}")
                    genotypes[snp] = dosage(record,samples)
            scans.append({"CHR":chrom,"start":start,"stop":stop,"records_scanned":n_records,
                          "reference":"1000G_30x_GRCh38" if highcov else "1000G_Phase3_GRCh37",
                          "url":url})
            print(f"chr{chrom}:{start}-{stop} records={n_records}",flush=True)
        vcf.close()
    matches = []
    for row in candidates:
        snp = row["SNP"]
        matches.append({"SNP":snp,"CHR":row["CHR"],"BP":row["BP"],
                        "status":"EXACT_ALLELE_GENOTYPE" if snp in genotypes else "NO_EXACT_ALLELE_GENOTYPE",
                        "reference":mapping[snp]["status"],"target_BP":mapping[snp]["target_bp"]})
    leads = []
    assignments = []
    for chrom in sorted({r["CHR"] for r in candidates},key=int):
        ordered = sorted((r for r in candidates if r["CHR"]==chrom),key=lambda r:(float(r["P_PLACO"]),int(r["BP"]),r["SNP"]))
        remaining = {r["SNP"] for r in ordered if r["SNP"] in genotypes}
        for lead in ordered:
            snp = lead["SNP"]
            if snp not in remaining:
                continue
            remaining.remove(snp)
            leads.append({"CHR":chrom,"SNP":snp,"BP":lead["BP"],"P_PLACO":lead["P_PLACO"],"reference":mapping[snp]["status"]})
            assignments.append({"CHR":chrom,"SNP":snp,"BP":lead["BP"],"lead_SNP":snp,"status":"LEAD","r2_to_lead":1.0})
            for other in ordered:
                second = other["SNP"]
                if second not in remaining or abs(int(other["BP"])-int(lead["BP"])) > 500_000:
                    continue
                _,r2,n,_,_ = pearson_ld(genotypes[snp],genotypes[second])
                if math.isfinite(r2) and n>=30 and r2>=0.1:
                    remaining.remove(second)
                    assignments.append({"CHR":chrom,"SNP":second,"BP":other["BP"],"lead_SNP":snp,"status":"LD_CLUMPED","r2_to_lead":r2})
    for row in candidates:
        if row["SNP"] not in genotypes:
            assignments.append({"CHR":row["CHR"],"SNP":row["SNP"],"BP":row["BP"],"lead_SNP":"","status":"NO_EXACT_ALLELE_GENOTYPE","r2_to_lead":""})
    if len(assignments) != len(candidates):
        raise ValueError("1000G B assignment accounting incomplete")
    ukb_leads = {s for locus in loci for s in locus["lead_variants"].split(";")}
    observed_leads = {x["SNP"] for x in leads}
    OUT.mkdir(parents=True,exist_ok=False)
    emit(OUT / "variant_matches.tsv",["SNP","CHR","BP","status","reference","target_BP"],matches)
    emit(OUT / "variant_assignments.tsv",["CHR","SNP","BP","lead_SNP","status","r2_to_lead"],assignments)
    emit(OUT / "independent_leads.tsv",["CHR","SNP","BP","P_PLACO","reference"],leads)
    result = {
        "analysis_id":"brain6-admitted-B-independent-1000G-LD-v1",
        "status":"PASS_DIAGNOSTIC_1000G_LD_AUDIT",
        "five_track_loci_receipt_sha256":sha(LOCI / "provenance.json"),
        "B_ledger_sha256":sha(LEDGER),"panel_sha256":sha(PANEL),"chain_sha256":sha(CHAIN),
        "auditor_sha256":sha(Path(__file__)),"sample_count_EUR":len(samples),
        "n_candidates":len(candidates),"n_exact_allele_genotypes":len(genotypes),
        "n_no_exact_allele_genotypes":len(candidates)-len(genotypes),
        "UKB_frozen_rule_leads":sorted(ukb_leads),"independent_1000G_leads":sorted(observed_leads),
        "lead_sets_identical":ukb_leads==observed_leads,
        "scans":scans,
        "outputs_sha256":{p.name:sha(p) for p in (OUT / "variant_matches.tsv",OUT / "variant_assignments.tsv",OUT / "independent_leads.tsv")},
        "interpretation":"Diagnostic independent 1000G EUR LD check only. The frozen locus rule uses the pinned UKB/LAVA reference; any 1000G difference is reported, not substituted into the frozen rule."
    }
    (OUT / "provenance.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status","n_candidates","n_exact_allele_genotypes","n_no_exact_allele_genotypes","lead_sets_identical")}))


if __name__ == "__main__":
    main()
