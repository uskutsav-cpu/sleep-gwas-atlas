#!/usr/bin/env python3
"""Scan a streamed GPL6984 SOFT platform table without saving the large file."""
import csv
import json
import sys
from pathlib import Path
from pyliftover import LiftOver

requests = list(csv.DictReader(open("outputs/ncds_score_variant_request.csv")))
ids = {r["variant"] for r in requests}
chain = LiftOver("work/hg19ToHg18.over.chain.gz")
coordinates = {}
for r in requests:
    mapped = chain.convert_coordinate("chr" + r["chr37"], int(r["pos37"]) - 1) or []
    if len(mapped) == 1 and mapped[0][2] == "+":
        key = (mapped[0][0].removeprefix("chr"), str(mapped[0][1] + 1))
        coordinates.setdefault(key, set()).add(r["variant"])
started = finished = False
header = None
hits = []
for line in sys.stdin:
    if line.startswith("!platform_table_begin"):
        started = True
        continue
    if line.startswith("!platform_table_end"):
        finished = True
        break
    if not started:
        continue
    if header is None:
        header = line.rstrip("\r\n").split("\t")
        continue
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != len(header):
        continue
    row = dict(zip(header, fields))
    by_name = row.get("SNP_ID") in ids
    by_coordinate = row.get("GenomeBuild") == "36" and (row.get("Chr"), row.get("MapInfo")) in coordinates
    if by_name or by_coordinate:
        result = {k: row.get(k) for k in ("ID", "Name", "SNP_ID", "GenomeBuild", "Chr", "MapInfo", "SNP", "IlmnStrand", "SourceStrand", "TopGenomicSeq")}
        result["candidate_coordinate_ids"] = sorted(coordinates.get((row.get("Chr"), row.get("MapInfo")), ())) if by_coordinate else []
        hits.append(result)
report = {"source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6984&targ=self&form=text&view=full",
          "complete_platform_table": finished, "candidate_count": len(ids), "matches": hits,
          "limitations": "Probe inventory only, no target genotypes, QC, imputation or outcome validation. Strand and build require harmonization."}
Path("outputs/healthabc_geo_platform_coordinate_screen.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({"complete": finished, "unique_matched_rsids": len({x["SNP_ID"] for x in hits if x["SNP_ID"] in ids}), "unique_coordinate_candidates": len({v for x in hits for v in x["candidate_coordinate_ids"]}), "rows": len(hits)}))
if not finished:
    sys.exit("Incomplete platform stream; absence is unknown")
