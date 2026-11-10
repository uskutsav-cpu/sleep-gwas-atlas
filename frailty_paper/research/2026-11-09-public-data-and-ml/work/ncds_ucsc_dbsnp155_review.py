#!/usr/bin/env python3
"""Check unresolved count-score positions against public UCSC dbSNP 155.

This checks a reference record, not the availability or imputation quality of
any variant in an individual-level validation cohort.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import subprocess
import urllib.parse
from pathlib import Path

TARGETS = {"rs28528997", "rs2815764", "rs2155292", "rs817509", "rs2838788",
           "rs12671330", "rs2604551", "rs9845387", "rs12524473"}
parser = argparse.ArgumentParser()
parser.add_argument("--all", action="store_true", help="Audit all 22 markers")
args = parser.parse_args()
with open("outputs/ncds_score_variant_request.csv", newline="") as f:
    rows = [r for r in csv.DictReader(f) if r["score"] == "22_count" and (args.all or r["variant"] in TARGETS)]
assert len(rows) == (22 if args.all else len(TARGETS))

def check(row):
    pos = int(row["pos38"])
    url = ("https://api.genome.ucsc.edu/getData/track?" + urllib.parse.urlencode({
        "genome": "hg38", "track": "dbSnp155", "chrom": "chr" + row["chr38"],
        "start": pos - 1, "end": pos}))
    result = {"variant": row["variant"], "chrom": row["chr38"], "pos38": pos,
              "requested_alleles": [row["counted_allele_38_candidate"], row["other_allele_38_candidate"]],
              "url": url}
    try:
        body = subprocess.run(["curl", "-fsSL", "--max-time", "20", url],
                              check=True, capture_output=True).stdout
        obj = json.loads(body)
        hits = [h for h in obj.get("dbSnp155", []) if h.get("name") == row["variant"]
                and h.get("chromStart") == pos - 1 and h.get("chromEnd") == pos]
        result.update(status="ok", response_sha256=hashlib.sha256(body).hexdigest(),
                      exact_records=len(hits))
        if len(hits) == 1:
            h = hits[0]
            alleles = [h["ref"]] + [x for x in h["alts"].split(",") if x]
            result.update(reference_alleles=alleles,
                          both_requested_alleles_listed=set(result["requested_alleles"]) <= set(alleles),
                          multiallelic=len(set(alleles)) > 2)
    except Exception as error:
        result.update(status="error", error=str(error))
    return result

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    out = list(executor.map(check, rows))
destination = "outputs/ncds_ucsc_dbsnp155_all22.json" if args.all else "outputs/ncds_ucsc_dbsnp155_review.json"
Path(destination).write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"requested": len(out), "ok": sum(x["status"] == "ok" for x in out),
                  "exact_pair": sum(x.get("both_requested_alleles_listed", False) for x in out),
                  "multiallelic": sum(x.get("multiallelic", False) for x in out),
                  "errors": [x["variant"] for x in out if x["status"] != "ok"]}, indent=2))
