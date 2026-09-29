#!/usr/bin/env python3
"""Independent, exact-allele LD audit of the 19 partial PLACO regions.

Reads only requested GRCh37 intervals from the official 1000 Genomes Phase 3
indexed VCFs. The 503 European samples are selected from the official panel.
No canonical LAVA, PLACO, or clumping output is changed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pysam


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE_DIR = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1")
BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF_NAME = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, values: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merged_intervals(regions: list[dict[str, str]]) -> list[tuple[int, int]]:
    ordered = sorted((int(x["START"]), int(x["STOP"])) for x in regions)
    merged: list[list[int]] = []
    for start, stop in ordered:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(stop, merged[-1][1])
        else:
            merged.append([start, stop])
    return [(a, b) for a, b in merged]


def dosage(record: pysam.VariantRecord, sample_ids: list[str]) -> np.ndarray:
    result = np.empty(len(sample_ids), dtype=np.float64)
    for index, sample in enumerate(sample_ids):
        gt = record.samples[sample]["GT"]
        result[index] = float(sum(gt)) if gt is not None and None not in gt else math.nan
    return result


def pearson_ld(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int, float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    n = int(valid.sum())
    if n < 30:
        return math.nan, math.nan, n, math.nan, math.nan
    a, b = x[valid], y[valid]
    maf_a, maf_b = min(a.mean() / 2, 1 - a.mean() / 2), min(b.mean() / 2, 1 - b.mean() / 2)
    if a.std() == 0 or b.std() == 0:
        return math.nan, math.nan, n, maf_a, maf_b
    r = float(np.corrcoef(a, b)[0, 1])
    return r, r * r, n, float(maf_a), float(maf_b)


def analyze(reference_dir: Path, out: Path) -> None:
    source_loci = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_loci.tsv"
    source_variants = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv"
    source_signed = ROOT / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
    source_panel = reference_dir / "integrated_call_samples_v3.20130502.ALL.panel"
    loci, candidates = rows(source_loci), rows(source_variants)
    signed = {(x["pair_id"], x["SNP"]): x for x in rows(source_signed)}
    sample_ids = [x["sample"] for x in rows(source_panel) if x["super_pop"] == "EUR"]
    if len(loci) != 19 or len(sample_ids) != 503:
        raise ValueError(f"Unexpected region/sample count: {len(loci)}, {len(sample_ids)}")
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Duplicate sample IDs in official panel")

    by_chr: dict[str, list[dict[str, str]]] = defaultdict(list)
    for locus in loci:
        by_chr[locus["CHR"]].append(locus)
    wanted: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    for candidate in candidates:
        wanted[(candidate["CHR"], int(candidate["BP"]))].append(candidate)

    matches: dict[tuple[str, str], dict] = {}
    genotypes: dict[tuple[str, str], np.ndarray] = {}
    query_receipts: list[dict] = []
    for chrom in sorted(by_chr, key=int):
        url = f"{BASE_URL}/{VCF_NAME.format(chrom=chrom)}"
        vcf = pysam.VariantFile(url)
        vcf.subset_samples(sample_ids)
        for start, stop in merged_intervals(by_chr[chrom]):
            n_records, n_target_positions = 0, 0
            for record in vcf.fetch(chrom, start - 1, stop):
                n_records += 1
                target_rows = wanted.get((chrom, record.pos), [])
                if not target_rows:
                    continue
                n_target_positions += 1
                for candidate in target_rows:
                    key = candidate["pair_id"], candidate["SNP"]
                    ref = signed.get(key)
                    if ref is None:
                        continue  # unmatched candidates never define an exact-allele LD edge
                    expected = {ref["reference_A1"], ref["reference_A2"]}
                    if len(record.alts or ()) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                        continue
                    observed = {record.ref, record.alts[0]}
                    if expected != observed:
                        continue
                    if key in matches:
                        raise ValueError(f"More than one exact allele record for {key}")
                    oriented = dosage(record, sample_ids)
                    if ref["reference_A1"] == record.ref:
                        oriented = 2 - oriented
                    matches[key] = {"vcf_id": record.id or "", "vcf_ref": record.ref,
                                    "vcf_alt": record.alts[0], "allele_match": "EXACT_UNORDERED",
                                    "n_nonmissing": int(np.isfinite(oriented).sum())}
                    genotypes[key] = oriented
            query_receipts.append({"chrom": chrom, "start": start, "stop": stop,
                                   "vcf_url": url, "records_scanned": n_records,
                                   "records_at_candidate_positions": n_target_positions})
            print(f"chr{chrom}:{start}-{stop} scanned {n_records}; target-position records {n_target_positions}", flush=True)
        vcf.close()

    variant_rows: list[dict] = []
    edge_rows: list[dict] = []
    for candidate in candidates:
        key = candidate["pair_id"], candidate["SNP"]
        match = matches.get(key)
        lead_key = candidate["pair_id"], candidate["lead_SNP"]
        record = {"locus_id": candidate["locus_id"], "pair_id": candidate["pair_id"],
                  "SNP": candidate["SNP"], "CHR": candidate["CHR"], "BP": candidate["BP"],
                  "candidate_status": candidate["candidate_status"],
                  "locked_reference_status": candidate["reference_status"],
                  "vcf_match_status": "EXACT_ALLELES" if match else "NO_EXACT_ALLELE_MATCH",
                  "vcf_ref": match["vcf_ref"] if match else "",
                  "vcf_alt": match["vcf_alt"] if match else "",
                  "vcf_id": match["vcf_id"] if match else "",
                  "n_nonmissing": match["n_nonmissing"] if match else ""}
        variant_rows.append(record)
        if candidate["candidate_status"] == "LEAD":
            continue
        ld = (math.nan,) * 5
        if key in genotypes and lead_key in genotypes:
            ld = pearson_ld(genotypes[key], genotypes[lead_key])
        edge_rows.append({"locus_id": candidate["locus_id"], "pair_id": candidate["pair_id"],
                          "SNP": candidate["SNP"], "lead_SNP": candidate["lead_SNP"],
                          "CHR": candidate["CHR"], "BP": candidate["BP"],
                          "candidate_vcf_match": int(key in genotypes),
                          "lead_vcf_match": int(lead_key in genotypes),
                          "lava_r2": candidate["r2_to_lead"],
                          "eur_r": ld[0], "eur_r2": ld[1], "eur_n": ld[2],
                          "eur_maf_candidate": ld[3], "eur_maf_lead": ld[4]})

    lead_rows: list[dict] = []
    leads = [x for x in candidates if x["candidate_status"] == "LEAD"]
    for i, a in enumerate(leads):
        for b in leads[i + 1:]:
            if a["CHR"] != b["CHR"] or a["locus_id"] == b["locus_id"]:
                continue
            if abs(int(a["BP"]) - int(b["BP"])) > 2_000_000:
                continue
            key_a, key_b = (a["pair_id"], a["SNP"]), (b["pair_id"], b["SNP"])
            ld = (math.nan,) * 5
            if key_a in genotypes and key_b in genotypes:
                ld = pearson_ld(genotypes[key_a], genotypes[key_b])
            lead_rows.append({"locus_a": a["locus_id"], "locus_b": b["locus_id"],
                              "pair_a": a["pair_id"], "pair_b": b["pair_id"],
                              "lead_a": a["SNP"], "lead_b": b["SNP"],
                              "distance_bp": abs(int(a["BP"]) - int(b["BP"])),
                              "both_vcf_match": int(key_a in genotypes and key_b in genotypes),
                              "eur_r": ld[0], "eur_r2": ld[1], "eur_n": ld[2]})

    summary = []
    for locus in loci:
        lid = locus["locus_id"]
        selected = [x for x in variant_rows if x["locus_id"] == lid]
        edges = [x for x in edge_rows if x["locus_id"] == lid]
        measured = [x for x in edges if math.isfinite(x["eur_r2"])]
        leads_here = [x for x in selected if x["candidate_status"] == "LEAD"]
        discordant = [x for x in measured if math.isfinite(float(x["lava_r2"])) and
                      (float(x["lava_r2"]) >= 0.1) != (x["eur_r2"] >= 0.1)]
        summary.append({"locus_id": lid, "pair_id": locus["pair_id"],
                        "CHR": locus["CHR"], "START": locus["START"], "STOP": locus["STOP"],
                        "lead_variants": locus["lead_variants"],
                        "n_candidates": len(selected),
                        "n_exact_vcf": sum(x["vcf_match_status"] == "EXACT_ALLELES" for x in selected),
                        "n_leads": len(leads_here),
                        "n_leads_exact_vcf": sum(x["vcf_match_status"] == "EXACT_ALLELES" for x in leads_here),
                        "n_edges_measured": len(measured),
                        "n_r2_threshold_discordant": len(discordant),
                        "min_eur_r2": min((x["eur_r2"] for x in measured), default=math.nan),
                        "max_eur_r2": max((x["eur_r2"] for x in measured), default=math.nan)})

    out.mkdir(parents=True, exist_ok=False)
    outputs = {
        "variant_matches.tsv": (variant_rows, list(variant_rows[0])),
        "candidate_lead_ld.tsv": (edge_rows, list(edge_rows[0])),
        "cross_region_lead_ld.tsv": (lead_rows, list(lead_rows[0]) if lead_rows else
                                     ["locus_a", "locus_b", "pair_a", "pair_b", "lead_a", "lead_b",
                                      "distance_bp", "both_vcf_match", "eur_r", "eur_r2", "eur_n"]),
        "region_summary.tsv": (summary, list(summary[0])),
    }
    for name, (values, fields) in outputs.items():
        write_rows(out / name, values, fields)
    receipt = {"analysis_id": "brain6_1000g_phase3_eur_ld_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "reference": "1000 Genomes Phase 3, GRCh37, 20130502 integrated VCF",
               "panel_url": f"{BASE_URL}/{source_panel.name}",
               "sample_selection": {"super_pop": "EUR", "n": len(sample_ids),
                                    "populations": sorted({x["pop"] for x in rows(source_panel) if x["super_pop"] == "EUR"})},
               "limitations": ["Different population/sample size from UK Biobank LAVA reference",
                               "An absent exact-allele lead cannot have its LD validated",
                               "No external trait association or replication is supplied by this genotype panel"],
               "query_receipts": query_receipts,
               "inputs_sha256": {str(x): sha256(x) for x in (source_loci, source_variants, source_signed, source_panel, Path(__file__))},
               "outputs_sha256": {name: sha256(out / name) for name in outputs}}
    (out / "provenance.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"regions": len(summary), "lead_matches": sum(x["n_leads_exact_vcf"] for x in summary),
                      "candidate_matches": sum(x["n_exact_vcf"] for x in summary),
                      "measured_edges": sum(x["n_edges_measured"] for x in summary)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--out", type=Path, default=ROOT / "brain6/results/loci/independent_1000g_phase3_eur_ld_v1")
    args = parser.parse_args()
    analyze(args.reference_dir, args.out)
