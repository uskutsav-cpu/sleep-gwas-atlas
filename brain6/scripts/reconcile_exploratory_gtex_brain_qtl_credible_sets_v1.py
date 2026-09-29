#!/usr/bin/env python3
"""Match frozen PLACO candidates to published GTEx brain QTL credible sets.

This reports molecular credible-set membership, never GWAS fine-mapping or
colocalization. All four datasets are fixed before inspecting matched values.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[2]
SIGNED = ROOT / "brain6/results/exploratory_five_track_v1"
V4 = ROOT / "brain6/results/exploratory_five_track_v4"
OUT = ROOT / "brain6/results/exploratory_five_track_v5/gtex_brain_qtl_credible_sets"
REGION_OUT = OUT.parent
CHAIN = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz")
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
DATASETS = {
    "QTD000171": ("brain_cortex", "gene_expression"),
    "QTD000176": ("brain_frontal_cortex", "gene_expression"),
    "QTD000175": ("brain_cortex", "leafcutter_splicing"),
    "QTD000180": ("brain_frontal_cortex", "leafcutter_splicing"),
}
CS_COLUMNS = ("molecular_trait_id", "gene_id", "cs_id", "variant", "rsid", "cs_size", "pip", "pvalue", "beta", "se", "z", "cs_min_r2", "region")


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
            raise ValueError(f"Receipt hash mismatch: {folder / name}")


def write(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if (OUT / "provenance.json").exists() or (REGION_OUT / "provenance.json").exists():
        raise FileExistsError("Completed v5 credible-set snapshot already exists")
    check_receipt(V4)
    check_receipt(V4 / "lead_eqtl_context")
    if not CHAIN.is_file():
        raise FileNotFoundError("Pinned GRCh37-to-GRCh38 chain is unavailable")
    metadata = {r["dataset_id"]: r for r in read(V4 / "lead_eqtl_context/official_dataset_metadata_r7.tsv")}
    for dataset_id, (sample_group, method) in DATASETS.items():
        m = metadata[dataset_id]
        if (m["study_id"], m["study_label"], m["sample_group"], m["quant_method"]) != (
            "QTS000015", "GTEx", sample_group, "ge" if method == "gene_expression" else "leafcutter"
        ):
            raise ValueError(f"Unexpected source identity: {dataset_id}")
    region = {r["candidate_locus_id"]: r["region_group"] for r in read(MEMBERS)}
    prior = read(V4 / "geographic_region_evidence.tsv")
    if len(prior) != 20 or len(region) != 25:
        raise ValueError("Protected region/candidate accounting changed")
    signed_paths = [SIGNED / "b_signed_ld/candidate_signed_ld.tsv", SIGNED / "fourpair_signed_ld/candidate_signed_ld.tsv"]
    candidates = [r for path in signed_paths for r in read(path)]
    if len(candidates) != 2686 or len({(r["pair_id"], r["SNP"]) for r in candidates}) != 2686:
        raise ValueError("Five-track candidate identity changed")
    lift = LiftOver(str(CHAIN))
    by_coordinate = defaultdict(list)
    coverage = []
    for candidate in candidates:
        item = {"region_group": region.get(candidate["locus_id"], "UNASSIGNED_NO_REFERENCE_VARIANT"),
                "pair_id": candidate["pair_id"], "locus_id": candidate["locus_id"], "snp": candidate["SNP"],
                "chr37": candidate["CHR"], "bp37": candidate["BP"],
                "candidate_status": candidate["candidate_status"], "signed_ld_status": candidate["ld_qc_status"],
                "reference_a1": candidate["reference_A1"], "reference_a2": candidate["reference_A2"],
                "chr38": "NA", "bp38": "NA", "liftover_status": "UNRESOLVED",
                "n_molecular_cs_memberships": "0", "max_molecular_pip": "NA",
                "analysis_label": "EXPLORATORY_QTL_ONLY", "lava_local_rg": "BLOCKED_LAVA"}
        if candidate["reference_A1"] == "NA" or candidate["reference_A2"] == "NA":
            item["liftover_status"] = "NO_EXACT_REFERENCE_ALLELES"
        else:
            mapped = [x for x in lift.convert_coordinate("chr" + candidate["CHR"], int(candidate["BP"]) - 1) if x[2] == "+"]
            if len(mapped) == 1:
                item["chr38"], item["bp38"] = mapped[0][0][3:], str(mapped[0][1] + 1)
                item["liftover_status"] = "UNIQUE_PLUS"
                by_coordinate[(mapped[0][0], mapped[0][1] + 1)].append(len(coverage))
        coverage.append(item)
    lift_counts = Counter(r["liftover_status"] for r in coverage)
    if lift_counts != {"UNIQUE_PLUS": 2662, "NO_EXACT_REFERENCE_ALLELES": 24}:
        raise ValueError(f"Unexpected build mapping/allele coverage: {lift_counts}")

    memberships = []
    file_records = {}
    for dataset_id, (sample_group, method) in DATASETS.items():
        path = OUT / "raw" / f"{dataset_id}.credible_sets.tsv.gz"
        if not path.is_file():
            raise FileNotFoundError(path)
        source_url = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/"
                      f"{dataset_id}/{dataset_id}.credible_sets.tsv.gz")
        row_count = 0
        with gzip.open(path, "rt", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if tuple(reader.fieldnames or ()) != CS_COLUMNS:
                raise ValueError(f"Official credible-set columns changed: {dataset_id}")
            for qtl in reader:
                row_count += 1
                parts = qtl["variant"].split("_", 3)
                if len(parts) != 4:
                    raise ValueError(f"Malformed molecular variant: {qtl['variant']}")
                indices = by_coordinate.get((parts[0], int(parts[1])), ())
                for index in indices:
                    candidate = coverage[index]
                    if qtl["rsid"] != candidate["snp"] or {parts[2], parts[3]} != {
                        candidate["reference_a1"], candidate["reference_a2"]
                    }:
                        continue
                    pip = float(qtl["pip"])
                    if not 0 <= pip <= 1:
                        raise ValueError(f"Invalid molecular PIP: {qtl['pip']}")
                    memberships.append({
                        "region_group": candidate["region_group"], "pair_id": candidate["pair_id"],
                        "locus_id": candidate["locus_id"], "snp": candidate["snp"],
                        "candidate_status": candidate["candidate_status"], "signed_ld_status": candidate["signed_ld_status"],
                        "chr37": candidate["chr37"], "bp37": candidate["bp37"],
                        "chr38": candidate["chr38"], "bp38": candidate["bp38"],
                        "dataset_id": dataset_id, "sample_group": sample_group, "qtl_method": method,
                        "gene_id": qtl["gene_id"], "molecular_trait_id": qtl["molecular_trait_id"],
                        "cs_id": qtl["cs_id"], "cs_size": qtl["cs_size"], "molecular_pip": qtl["pip"],
                        "molecular_cs_min_r2": qtl["cs_min_r2"], "molecular_variant": qtl["variant"],
                        "molecular_nominal_p": qtl["pvalue"],
                        "matching_rule": "EXACT_RSID_GRCH38_POSITION_UNORDERED_ALLELES",
                        "analysis_label": "EXPLORATORY_MOLECULAR_CS_MEMBERSHIP_ONLY", "lava_local_rg": "BLOCKED_LAVA",
                    })
        file_records[dataset_id] = {"source_url": source_url, "size_bytes": path.stat().st_size,
                                    "sha256": sha(path), "rows": row_count,
                                    "sample_group": sample_group, "qtl_method": method,
                                    "sample_size_metadata": metadata[dataset_id]["sample_size"]}
    if not memberships:
        raise ValueError("No matched molecular credible-set memberships")
    by_candidate = defaultdict(list)
    by_region = defaultdict(list)
    for entry in memberships:
        by_candidate[(entry["pair_id"], entry["snp"])].append(entry)
        by_region[entry["region_group"]].append(entry)
    for item in coverage:
        hits = by_candidate[(item["pair_id"], item["snp"])]
        item["n_molecular_cs_memberships"] = str(len(hits))
        item["max_molecular_pip"] = str(max(float(x["molecular_pip"]) for x in hits)) if hits else "NA"

    updated = []
    for old in prior:
        item = dict(old)
        region_id = item["geographic_region"]
        hits = by_region[region_id]
        distinct_memberships = {
            (x["dataset_id"], x["snp"], x["molecular_trait_id"], x["cs_id"], x["molecular_variant"]) for x in hits
        }
        item["gtex_brain_qtl_cs_distinct_memberships"] = str(len(distinct_memberships))
        item["gtex_brain_qtl_cs_candidate_variants"] = str(len({x["snp"] for x in hits}))
        item["gtex_brain_qtl_cs_distinct_sets"] = str(len({(x["dataset_id"], x["cs_id"]) for x in hits}))
        for method, label in (("gene_expression", "eqtl"), ("leafcutter_splicing", "sqtl")):
            subset = [x for x in hits if x["qtl_method"] == method]
            best = max(subset, key=lambda x: float(x["molecular_pip"])) if subset else None
            item[f"gtex_brain_{label}_cs_max_molecular_pip"] = best["molecular_pip"] if best else "NO_EXACT_CANDIDATE_CS_OVERLAP"
            item[f"gtex_brain_{label}_cs_max_pip_gene_id"] = best["gene_id"] if best else "NO_EXACT_CANDIDATE_CS_OVERLAP"
        item["gtex_brain_qtl_cs_interpretation"] = "EXPLORATORY_MOLECULAR_FINE_MAP_OVERLAP_ONLY"
        updated.append(item)
    if len(updated) != 20 or sum(int(r["n_pair_specific_candidates"]) for r in updated) != 25:
        raise ValueError("Twenty-region/25-candidate preservation failed")
    if sum(int(r["n_molecular_cs_memberships"]) for r in coverage) != len(memberships):
        raise ValueError("Candidate membership counts failed to reconcile")

    write(OUT / "candidate_molecular_credible_set_memberships.tsv", memberships, list(memberships[0]))
    write(OUT / "candidate_qtl_finemap_coverage.tsv", coverage, list(coverage[0]))
    write(REGION_OUT / "geographic_region_evidence.tsv", updated, list(updated[0]))
    report = [
        "# Five-track exploratory molecular credible-set overlap",
        "",
        "**EXPLORATORY. Zero LAVA-confirmed regions and zero final shared-locus tiers.**",
        "",
        "The four prespecified GTEx v8 eQTL Catalogue datasets cover brain cortex and DLPFC",
        "gene expression and Leafcutter splicing. All 2,686 pair-specific candidate rows remain",
        "in the coverage file. Of those, 2,662 have exact UKB reference alleles and unique",
        "plus-strand GRCh37-to-GRCh38 maps; 24 lack exact reference alleles. Matching requires",
        "the source rsID, lifted coordinate, and unordered allele pair to agree. All four raw",
        "credible-set files are retained and hashed. Molecular PIP is the QTL model's PIP only.",
        "",
        "| Region | Candidates in molecular credible sets | Distinct molecular sets | Highest eQTL PIP | Highest sQTL PIP |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in updated:
        report.append(f"| {r['geographic_region']} | {r['gtex_brain_qtl_cs_candidate_variants']} | "
                      f"{r['gtex_brain_qtl_cs_distinct_sets']} | {r['gtex_brain_eqtl_cs_max_molecular_pip']} | "
                      f"{r['gtex_brain_sqtl_cs_max_molecular_pip']} |")
    report.extend([
        "",
        "These are exact PLACO-candidate overlaps with pre-existing **molecular** credible sets.",
        "A PLACO candidate has no established GWAS PIP here, and an overlapping molecular",
        "credible set is not colocalization or proof of a shared causal variant. The maximum",
        "PIP values are descriptive and must not rank shared loci or genes. No absence of QTL",
        "is inferred from an unmatched candidate. Trait fine-mapping, trait-trait coloc,",
        "trait-QTL coloc, independent locus-specific replication, and LAVA local-rg remain",
        "unestablished; final region tiers stay BLOCKED_LAVA.", "",
    ])
    (REGION_OUT / "region_evidence_report.md").write_text("\n".join(report))
    receipt = {"analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA",
               "source_v4_receipt_sha256": sha(V4 / "provenance.json"),
               "source_v4_lead_eqtl_receipt_sha256": sha(V4 / "lead_eqtl_context/provenance.json"),
               "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in signed_paths + [MEMBERS]},
               "chain_path": str(CHAIN), "chain_sha256": sha(CHAIN),
               "metadata_sha256": sha(V4 / "lead_eqtl_context/official_dataset_metadata_r7.tsv"),
               "source_files": file_records,
               "candidate_rows": len(coverage), "molecular_membership_rows": len(memberships),
               "output_sha256": {p.name: sha(p) for p in (
                   OUT / "candidate_molecular_credible_set_memberships.tsv", OUT / "candidate_qtl_finemap_coverage.tsv")}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    parent_receipt = {"analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA",
                      "source_molecular_receipt_sha256": sha(OUT / "provenance.json"),
                      "source_v4_receipt_sha256": sha(V4 / "provenance.json"),
                      "output_sha256": {p.name: sha(p) for p in (
                          REGION_OUT / "geographic_region_evidence.tsv", REGION_OUT / "region_evidence_report.md")}}
    (REGION_OUT / "provenance.json").write_text(json.dumps(parent_receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"candidates": len(coverage), "molecular_membership_rows": len(memberships),
                      "regions": len(updated), "regions_with_exact_qtl_cs_overlap": sum(
                          int(r["gtex_brain_qtl_cs_candidate_variants"]) > 0 for r in updated)}, sort_keys=True))


if __name__ == "__main__":
    main()
