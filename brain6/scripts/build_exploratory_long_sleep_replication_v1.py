#!/usr/bin/env python3
"""Source-bound, phenotype-adjacent long-sleep lookups for frozen Brain6 leads."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import zipfile
from collections import Counter
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def alignment(original, other):
    if not original:
        return "ORIGINAL_VARIANT_UNAVAILABLE", "UNAVAILABLE"
    oa1, oa0 = original["ALLELE1"].upper(), original["ALLELE0"].upper()
    ea1, ea2 = other["A1"].upper(), other["A2"].upper()
    if {oa1, oa0} != {ea1, ea2}:
        return "ALLELE_MISMATCH", "UNAVAILABLE"
    if {oa1, oa0} in ({"A", "T"}, {"C", "G"}):
        return "PALINDROMIC_UNRESOLVED", "UNAVAILABLE"
    flip = 1 if oa1 == ea1 else -1
    sign = float(original["BETA_LONGSLEEP"]) * float(other["Effect"]) * flip
    return ("MATCH" if flip == 1 else "SWAPPED"), ("CONCORDANT" if sign > 0 else "DISCORDANT" if sign < 0 else "ZERO_EFFECT")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--dashti-zip", required=True)
    parser.add_argument("--yale-gzip", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    repo, out = Path(args.repo), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    base = repo / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1"
    regions = read_tsv(base / "geographic_region_groups.tsv")
    candidates = read_tsv(base / "pair_candidate_members.tsv")
    variant_rows = read_tsv(repo / "brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_variants.tsv")
    coordinates = {row["SNP"]: (row["CHR"], row["BP"]) for row in variant_rows}
    leads = sorted({lead for row in candidates for lead in row["lead_variants"].split(";")})
    assert len(regions) == 20 and len(candidates) == 25 and len(leads) == 26
    assert all(lead in coordinates for lead in leads)
    wanted = set(leads)
    original = {}
    with zipfile.ZipFile(args.dashti_zip) as archive:
        with archive.open("longsumstats.txt") as handle:
            reader = csv.DictReader((line.decode("utf-8") for line in handle), delimiter="\t")
            for row in reader:
                if row["SNP"] in wanted and (row["CHR"], row["BP"]) == coordinates[row["SNP"]]:
                    original[row["SNP"]] = row
    region_by_chr = {}
    for row in regions:
        region_by_chr.setdefault(row["chromosome"], []).append(row)
    coord_to_lead = {(chrom, pos): lead for lead, (chrom, pos) in coordinates.items() if lead in wanted}
    yale = {}
    region_min = {}
    n_values = Counter()
    rows_scanned = 0
    with gzip.open(args.yale_gzip, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == ["rsid", "CHR", "BP", "A1", "A2", "Freq1", "N", "Effect", "P"]
        for row in reader:
            rows_scanned += 1
            n_values[row["N"]] += 1
            chrom, pos = row["CHR"], int(row["BP"])
            lead = coord_to_lead.get((chrom, row["BP"]))
            if lead == row["rsid"]:
                yale[lead] = row
            try:
                p = float(row["P"])
            except ValueError:
                continue
            if not math.isfinite(p) or p < 0 or p > 1:
                continue
            for region in region_by_chr.get(chrom, []):
                if int(region["start"]) <= pos <= int(region["stop"]):
                    key = region["region_group"]
                    if key not in region_min or p < region_min[key][0]:
                        region_min[key] = (p, row["rsid"])
    lead_rows = []
    for lead in leads:
        source = original.get(lead)
        match = yale.get(lead)
        allele, direction = alignment(source, match) if match else ("YALE_VARIANT_UNAVAILABLE", "UNAVAILABLE")
        lead_rows.append({
            "lead_variant": lead, "chromosome": coordinates[lead][0], "position_grch37": coordinates[lead][1],
            "dashti_present": bool(source), "dashti_a1": source["ALLELE1"] if source else "",
            "dashti_a0": source["ALLELE0"] if source else "",
            "dashti_beta": source["BETA_LONGSLEEP"] if source else "",
            "yale_present": bool(match), "yale_a1": match["A1"] if match else "",
            "yale_a2": match["A2"] if match else "",
            "yale_effect": match["Effect"] if match else "",
            "yale_p": match["P"] if match else "", "yale_n": match["N"] if match else "",
            "allele_alignment": allele, "direction_vs_dashti": direction,
        })
    lead_lookup = {row["lead_variant"]: row for row in lead_rows}
    candidate_rows = []
    for row in candidates:
        lrows = [lead_lookup[x] for x in row["lead_variants"].split(";")]
        valid = [x for x in lrows if x["yale_present"] and x["allele_alignment"] in ("MATCH", "SWAPPED")]
        exact = min(valid, key=lambda x: float(x["yale_p"])) if valid else None
        minimum = region_min.get(row["region_group"])
        if exact and float(exact["yale_p"]) < 0.05 / 25 and exact["direction_vs_dashti"] == "CONCORDANT":
            label = "EXPLORATORY_EXACT_VARIANT_REPLICATION"
        elif minimum and minimum[0] < 5e-8:
            label = "EXPLORATORY_LOCUS_REPLICATION"
        elif exact and exact["direction_vs_dashti"] == "CONCORDANT":
            label = "EXPLORATORY_DIRECTION_ONLY"
        elif exact:
            label = "NO_REPLICATION"
        else:
            label = "UNAVAILABLE"
        candidate_rows.append({
            "region_group": row["region_group"], "pair_id": row["pair_id"],
            "candidate_locus_id": row["candidate_locus_id"], "lead_variants": row["lead_variants"],
            "best_exact_lead": exact["lead_variant"] if exact else "",
            "best_exact_p": exact["yale_p"] if exact else "",
            "best_exact_direction": exact["direction_vs_dashti"] if exact else "",
            "minimum_region_p": minimum[0] if minimum else "",
            "minimum_region_variant": minimum[1] if minimum else "",
            "classification": label,
            "phenotype_relation": "EXPLORATORY_ADJACENT_GE10_VS_7_TO_8;UKB_OVERLAP",
        })
    region_rows = []
    for row in regions:
        members = [c for c in candidate_rows if c["region_group"] == row["region_group"]]
        minimum = region_min.get(row["region_group"])
        region_rows.append({"region_group": row["region_group"], "chromosome": row["chromosome"],
                            "start": row["start"], "stop": row["stop"],
                            "n_pair_specific_candidates": len(members), "pairs": row["pairs"],
                            "minimum_yale_p": minimum[0] if minimum else "",
                            "minimum_yale_variant": minimum[1] if minimum else "",
                            "candidate_classifications": ";".join(x["classification"] for x in members),
                            "confirmatory_status": "NOT_PROMOTED"})
    write_tsv(out / "lead_variant_lookups.tsv", lead_rows, list(lead_rows[0]))
    write_tsv(out / "pair_candidate_replication.tsv", candidate_rows, list(candidate_rows[0]))
    write_tsv(out / "geographic_region_replication.tsv", region_rows, list(region_rows[0]))
    provenance = {"analysis_id": "exploratory_long_sleep_replication_v1", "source": "Austin-Zimmerman 2023 UKB+MVP EUR long sleep >=10h vs 7-8h",
                  "source_url": "https://hugesumstats.yale.edu/dl/eur_meta_long_dbgap.txt.gz",
                  "source_sha256": sha256(args.yale_gzip), "dashti_archive_sha256": sha256(args.dashti_zip),
                  "rows_scanned": rows_scanned, "n_values": dict(n_values), "regions": len(region_rows),
                  "pair_candidates": len(candidate_rows), "distinct_leads": len(lead_rows),
                  "exact_yale_leads": len(yale), "exact_dashti_leads": len(original),
                  "classification_counts": dict(Counter(x["classification"] for x in candidate_rows)),
                  "thresholds_frozen": {"exact_variant": 0.05 / 25, "locus_descriptive": 5e-8}}
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
