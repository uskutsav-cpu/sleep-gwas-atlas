#!/usr/bin/env python3
"""Snapshot exact lead-level GTEx brain eQTL context, without coloc claims.

Run with the pinned independent-reference virtual environment, which provides
pyliftover and pysam. Query files are durable checkpoints and are never
overwritten; a completed provenance receipt makes the whole output immutable.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pysam
from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/exploratory_five_track_v1"
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
OUT = ROOT / "brain6/results/exploratory_five_track_v4/lead_eqtl_context"
CHAIN = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz")
METADATA_URL = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r7.tsv"
DATASETS = {
    "QTD000171": ("brain_cortex", "ge"),
    "QTD000176": ("brain_frontal_cortex", "ge"),
}
SOURCE_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"
COLUMNS = ("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
MIN_QUERY_GAP_SECONDS = 2.5


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def save_json_once(path: Path, data: dict) -> None:
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write("\n")


def download_metadata(path: Path) -> None:
    if path.exists():
        return
    response = subprocess.run(
        ["curl", "-L", "--fail", "--silent", "--show-error", "--max-time", "40", METADATA_URL],
        check=True, capture_output=True,
    ).stdout
    if not response.startswith(b"study_id\tdataset_id\t"):
        raise ValueError("Official metadata header mismatch")
    with path.open("xb") as stream:
        stream.write(response)


def allele_status(gwas_a1: str, gwas_a2: str, ref: str, alt: str) -> str:
    if len({gwas_a1, gwas_a2, ref, alt} - {"A", "C", "G", "T"}):
        return "NON_SNP_ALLELES"
    if {gwas_a1, gwas_a2} == {ref, alt}:
        return "EXACT_UNORDERED"
    complement = str.maketrans("ACGT", "TGCA")
    if {gwas_a1, gwas_a2} == {ref.translate(complement), alt.translate(complement)}:
        return "COMPLEMENT_UNORDERED"
    return "ALLELE_MISMATCH"


def main() -> None:
    if (OUT / "provenance.json").exists():
        raise FileExistsError("Completed lead-eQTL snapshot already exists")
    if not CHAIN.is_file():
        raise FileNotFoundError(f"Required pinned liftOver chain unavailable: {CHAIN}")
    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = OUT / "raw_queries"
    raw_dir.mkdir(exist_ok=True)
    metadata = OUT / "official_dataset_metadata_r7.tsv"
    download_metadata(metadata)
    meta_rows = {r["dataset_id"]: r for r in read_tsv(metadata)}
    for dataset_id, (sample_group, method) in DATASETS.items():
        row = meta_rows[dataset_id]
        if (row["study_id"], row["study_label"], row["sample_group"], row["quant_method"]) != (
            "QTS000015", "GTEx", sample_group, method
        ):
            raise ValueError(f"Official dataset identity changed for {dataset_id}")

    region = {r["candidate_locus_id"]: r["region_group"] for r in read_tsv(MEMBERS)}
    if len(region) != 25:
        raise ValueError("Expected 25 protected candidate loci")
    signed_paths = [BASE / "b_signed_ld/candidate_signed_ld.tsv", BASE / "fourpair_signed_ld/candidate_signed_ld.tsv"]
    leads = [r for path in signed_paths for r in read_tsv(path) if r["candidate_status"] == "LEAD"]
    if len(leads) != 27 or len({(r["CHR"], r["BP"], r["SNP"]) for r in leads}) != 26:
        raise ValueError("Expected 27 lead rows and 26 unique variants")
    if {r["locus_id"] for r in leads} != set(region):
        raise ValueError("Protected 25-locus membership changed")

    lift = LiftOver(str(CHAIN))
    coords = {}
    for lead in leads:
        key = (lead["CHR"], lead["BP"], lead["SNP"])
        mapped = lift.convert_coordinate("chr" + lead["CHR"], int(lead["BP"]) - 1)
        valid = [x for x in mapped if x[2] == "+"]
        coords[key] = (valid[0][0][3:], valid[0][1] + 1) if len(valid) == 1 else None
    last_fetch = 0.0
    query_hashes = {}
    datasets = {}
    for dataset_id in DATASETS:
        source_url = f"{SOURCE_BASE}/{dataset_id}/{dataset_id}.all.tsv.gz"
        datasets[dataset_id] = {"url": source_url, "index_url": source_url + ".tbi", "metadata": meta_rows[dataset_id]}
        with pysam.TabixFile(source_url) as tabix:
            for chrom37, bp37, snp in sorted(coords, key=lambda x: (int(x[0]), int(x[1]), x[2])):
                mapped = coords[(chrom37, bp37, snp)]
                if mapped is None:
                    continue
                chrom38, bp38 = mapped
                name = f"{dataset_id}_{snp}_{chrom37}_{bp37}.json"
                path = raw_dir / name
                query = {"dataset_id": dataset_id, "source_url": source_url, "chromosome38": chrom38,
                         "position38": bp38, "snp": snp, "chromosome37": chrom37, "position37": int(bp37)}
                if path.exists():
                    previous = json.loads(path.read_text())
                    if previous["query"] != query:
                        raise ValueError(f"Saved query identity changed: {path}")
                else:
                    gap = MIN_QUERY_GAP_SECONDS - (time.monotonic() - last_fetch)
                    if gap > 0:
                        time.sleep(gap)
                    raw_lines = list(tabix.fetch(chrom38, bp38 - 1, bp38))
                    last_fetch = time.monotonic()
                    parsed = [dict(zip(COLUMNS, line.split("\t"), strict=True)) for line in raw_lines]
                    if any(r["chromosome"] != chrom38 or int(r["position"]) != bp38 for r in parsed):
                        raise ValueError(f"Tabix returned a different coordinate: {query}")
                    save_json_once(path, {"query": query, "records": parsed})
                query_hashes[str(path.relative_to(OUT))] = sha(path)
                print(f"{dataset_id} {snp} chr{chrom38}:{bp38} cached", flush=True)

    associations = []
    summary = []
    for lead in leads:
        key = (lead["CHR"], lead["BP"], lead["SNP"])
        for dataset_id in DATASETS:
            mapped = coords[key]
            base = {"region_group": region[lead["locus_id"]], "pair_id": lead["pair_id"],
                    "locus_id": lead["locus_id"], "lead_snp": lead["SNP"],
                    "chr37": lead["CHR"], "bp37": lead["BP"],
                    "gwas_reference_a1": lead["reference_A1"], "gwas_reference_a2": lead["reference_A2"],
                    "dataset_id": dataset_id, "sample_group": DATASETS[dataset_id][0],
                    "sample_size_metadata": meta_rows[dataset_id]["sample_size"],
                    "chr38": mapped[0] if mapped else "NA", "bp38": str(mapped[1]) if mapped else "NA"}
            if mapped is None:
                summary.append(base | {"all_coordinate_rows": "0", "exact_rsid_allele_rows": "0", "status": "LIFTOVER_UNRESOLVED"})
                continue
            name = f"{dataset_id}_{lead['SNP']}_{lead['CHR']}_{lead['BP']}.json"
            raw = json.loads((raw_dir / name).read_text())
            matching = []
            for item in raw["records"]:
                status = allele_status(lead["reference_A1"], lead["reference_A2"], item["ref"], item["alt"])
                if item["rsid"] == lead["SNP"] and status in ("EXACT_UNORDERED", "COMPLEMENT_UNORDERED"):
                    matching.append(item)
                    associations.append(base | {"gene_id": item["gene_id"], "molecular_trait_id": item["molecular_trait_id"],
                                                "qtl_variant": item["variant"], "qtl_ref": item["ref"], "qtl_alt": item["alt"],
                                                "allele_status": status, "nominal_pvalue": item["pvalue"],
                                                "beta_alt": item["beta"], "se": item["se"], "maf": item["maf"],
                                                "analysis_label": "EXPLORATORY_LEAD_LEVEL_ONLY"})
            summary.append(base | {"all_coordinate_rows": str(len(raw["records"])),
                                   "exact_rsid_allele_rows": str(len(matching)),
                                   "status": "EXACT_LEAD_QTL_ROWS_PRESENT" if matching else "NO_EXACT_LEAD_QTL_ROW"})
    if len(summary) != 54:
        raise ValueError("Lead-by-dataset accounting incomplete")
    association_columns = list(associations[0]) if associations else list(summary[0]) + ["gene_id", "molecular_trait_id", "qtl_variant", "qtl_ref", "qtl_alt", "allele_status", "nominal_pvalue", "beta_alt", "se", "maf", "analysis_label"]
    write_tsv(OUT / "lead_eqtl_associations.tsv", associations, association_columns)
    write_tsv(OUT / "lead_eqtl_query_summary.tsv", summary, list(summary[0]))
    counts = Counter(r["status"] for r in summary)
    report = ["# GTEx v8 brain lead-eQTL lookup — exploratory", "",
              "This is an exact lead-coordinate lookup for the protected 27 lead rows (26 unique variants) in two",
              "prespecified eQTL Catalogue GTEx v8 datasets: brain cortex and brain frontal cortex (DLPFC).",
              "The official release-7 dataset metadata, GRCh37-to-GRCh38 chain, exact rsIDs, and unordered SNP alleles",
              "are bound in the provenance receipt. All nominal associations at matching lead coordinates are retained.", "",
              f"Queries: {len(query_hashes)}; lead-by-tissue rows: {len(summary)}; exact matched association rows: {len(associations)}.",
              f"Statuses: {dict(counts)}.", "",
              "**EXPLORATORY_LEAD_LEVEL_ONLY.** These rows do not establish locus-wide eQTL support, a significant",
              "QTL after gene/tissue multiplicity correction, colocalization, a causal gene, or tissue specificity.",
              "A missing exact lead row does not exclude an eQTL elsewhere in the locus. Splicing QTL and full-region",
              "molecular colocalization remain unrun; the Catalogue leafcutter folders inspected here expose selected",
              "conditional association files but not the full nominal all-association file required for unbiased coloc.",
              "The canonical LAVA decision and final tiers remain BLOCKED_LAVA.", ""]
    (OUT / "report.md").write_text("\n".join(report))
    inputs = {str(path.relative_to(ROOT)): sha(path) for path in signed_paths + [MEMBERS]}
    receipt = {"created_utc": datetime.now(timezone.utc).isoformat(), "analysis_label": "EXPLORATORY",
               "lava_local_rg": "BLOCKED_LAVA", "source_metadata_url": METADATA_URL,
               "source_metadata_sha256": sha(metadata), "chain_path": str(CHAIN), "chain_sha256": sha(CHAIN),
               "datasets": datasets, "input_sha256": inputs, "query_sha256": query_hashes,
               "pysam_version": pysam.__version__, "min_query_gap_seconds": MIN_QUERY_GAP_SECONDS,
               "output_sha256": {p.name: sha(p) for p in (OUT / "lead_eqtl_associations.tsv", OUT / "lead_eqtl_query_summary.tsv", OUT / "report.md")}}
    save_json_once(OUT / "provenance.json", receipt)
    print(json.dumps({"queries": len(query_hashes), "lead_by_tissue": len(summary),
                      "associations": len(associations), "statuses": dict(counts)}, sort_keys=True))


if __name__ == "__main__":
    main()
