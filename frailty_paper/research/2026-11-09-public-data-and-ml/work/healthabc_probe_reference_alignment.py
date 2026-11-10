#!/usr/bin/env python3
"""Align public probe flanks to hg18; does not inspect target genotypes."""
import concurrent.futures
import csv
import hashlib
import json
import re
import subprocess
import urllib.parse
from pathlib import Path

data = json.loads(Path("outputs/healthabc_geo_platform_coordinate_screen.json").read_text())
assert data["complete_platform_table"]
requests = list(csv.DictReader(open("outputs/ncds_score_variant_request.csv")))
model = {r["variant"]: {r["counted_allele_37"], r["other_allele_37"]} for r in requests}
comp = str.maketrans("ACGT", "TGCA")
def reverse_complement(s):
    return s.translate(comp)[::-1]
def run(row):
    result = {"variant": row["SNP_ID"], "chr36": row["Chr"], "pos36": row["MapInfo"]}
    m = re.fullmatch(r"([ACGT]+)\[([ACGT])/([ACGT])\]([ACGT]+)", (row["TopGenomicSeq"] or "").upper())
    if not m:
        return {**result, "status": "unsupported_flanking_sequence"}
    left, a, b, right = m.groups()
    # Extract a symmetric reference window so either orientation has correct extents.
    radius = max(len(left), len(right)); pos = int(row["MapInfo"])
    url = "https://api.genome.ucsc.edu/getData/sequence?" + urllib.parse.urlencode({"genome":"hg18", "chrom":"chr"+row["Chr"], "start":pos-1-radius, "end":pos+radius})
    result["url"] = url
    try:
        raw = subprocess.run(["curl", "-fsSL", "--max-time", "20", url], capture_output=True, check=True).stdout
        seq = json.loads(raw)["dna"].upper()
        if len(seq) != 2 * radius + 1:
            raise ValueError("unexpected reference window length")
        before, after = seq[:radius], seq[radius+1:]
        forward = before.endswith(left) and after.startswith(right)
        backward = before.endswith(reverse_complement(right)) and after.startswith(reverse_complement(left))
        orientation = "+" if forward and not backward else "-" if backward and not forward else "ambiguous_or_mismatch"
        aligned = {a,b} if orientation == "+" else {a.translate(comp),b.translate(comp)} if orientation == "-" else set()
        result.update(status="ok", reference_sha256=hashlib.sha256(raw).hexdigest(),
                      flank_orientation=orientation, top_alleles=sorted({a,b}),
                      reference_aligned_probe_alleles=sorted(aligned),
                      hg18_reference_base=seq[radius], reference_base_in_probe=seq[radius] in aligned,
                      project_alleles=sorted(model[row["SNP_ID"]]),
                      pair_matches_project=aligned == model[row["SNP_ID"]])
    except Exception as exc:
        result.update(status="error", error=str(exc))
    return result
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    results=list(pool.map(run,data["matches"]))
Path("outputs/healthabc_probe_reference_alignment.json").write_text(json.dumps(results,indent=2)+"\n")
print(json.dumps({"records":len(results),"successful":sum(r["status"]=="ok" for r in results),
                  "pair_matches":sum(r.get("pair_matches_project",False) for r in results),
                  "unresolved":[r["variant"] for r in results if r.get("flank_orientation") not in ("+","-")]},indent=2))
