#!/usr/bin/env python3
"""Fetch predeclared insomnia–ADHD candidate windows from official brain QTL files.

Requires pysam and pyliftover. Writes only versioned SSD slices with receipts;
completed slices are verified and never overwritten.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path

import pysam
from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[2]
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
CHAIN = Path("/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz")
OUT = Path("/Volumes/Extreme SSD/brain6-work/brain6-exploratory-functional-v1/qtl_regional")
DATASETS = {"QTD000171": "all", "QTD000175": "cc", "QTD000176": "all", "QTD000180": "cc"}
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"
HEADER = ("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf",
          "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def get_regions() -> list[tuple[str, int, int, str]]:
    with MEMBERS.open(newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r["pair_id"] == "insomnia__adhd"]
    if len(rows) != 6 or len({r["region_group"] for r in rows}) != 6:
        raise ValueError("Expected exactly six frozen insomnia–ADHD geographic regions")
    lift = LiftOver(str(CHAIN))
    out = []
    for row in rows:
        label = row["region_group"]
        chrom, span = label.removeprefix("chr").split(":")
        lo, hi = map(int, span.split("-"))
        hits = [lift.convert_coordinate("chr" + chrom, p - 1) for p in (lo, hi)]
        if any(len(x) != 1 or x[0][2] != "+" or x[0][0] != "chr" + chrom for x in hits):
            raise ValueError(f"Ambiguous endpoint liftover: {label}")
        start, stop = hits[0][0][1] + 1, hits[1][0][1] + 1
        if not 0 < start < stop:
            raise ValueError(f"Invalid lifted span: {label}")
        out.append((chrom, start, stop, label))
    return sorted(out, key=lambda r: (int(r[0]), r[1]))


def remote_metadata(url: str) -> dict[str, str]:
    completed = subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location", "--head",
                                "--max-time", "60", url], check=True, capture_output=True, text=True)
    header = {}
    for line in completed.stdout.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            header[key.strip().lower()] = value.strip()
    keys = ("ETag", "Last-Modified", "Content-Length", "Accept-Ranges")
    return {key: header.get(key.lower(), "") for key in keys}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    regions = get_regions()
    source_inputs = {str(MEMBERS): digest(MEMBERS), str(CHAIN): digest(CHAIN)}
    for dataset, suffix in DATASETS.items():
        url = f"{BASE}/{dataset}/{dataset}.{suffix}.tsv.gz"
        index_url = url + ".tbi"
        meta = remote_metadata(url)
        index_meta = remote_metadata(index_url)
        reader = pysam.TabixFile(url)
        if not all(r[0] in reader.contigs for r in regions):
            raise ValueError(f"Missing chromosome in {dataset}")
        for chrom, start, stop, label in regions:
            stem = f"{dataset}_{label.replace(':', '_').replace('-', '_')}"
            output = OUT / (stem + ".tsv.gz")
            receipt = OUT / (stem + ".receipt.json")
            if output.exists() and receipt.exists():
                old = json.loads(receipt.read_text())
                if old["sha256"] != digest(output) or old["source_metadata"]["ETag"] != meta["ETag"]:
                    raise ValueError(f"Existing slice changed: {output}")
                continue
            if output.exists() or receipt.exists():
                raise FileExistsError(f"Incomplete slice state: {output}")
            temp = OUT / (stem + ".partial")
            if temp.exists():
                raise FileExistsError(temp)
            count = 0
            molecular_traits = set()
            with gzip.open(temp, "wt", newline="") as f:
                f.write("\t".join(HEADER) + "\n")
                for line in reader.fetch(chrom, start - 1, stop):
                    fields = line.split("\t")
                    if len(fields) != len(HEADER):
                        raise ValueError(f"Unexpected QTL field count: {dataset} {label}")
                    count += 1
                    molecular_traits.add(fields[0])
                    f.write(line + "\n")
            os.replace(temp, output)
            record = {"dataset_id": dataset, "study_id": "QTS000015", "source_url": url,
                      "source_index_url": index_url, "source_metadata": meta, "index_metadata": index_meta,
                      "source_scope": "REMOTE_TABIX_REGION_NOT_FULL_FILE", "original_region_grch37": label,
                      "query_region_grch38": f"{chrom}:{start}-{stop}", "n_rows": count,
                      "n_molecular_traits": len(molecular_traits), "sha256": digest(output),
                      "size_bytes": output.stat().st_size, "input_sha256": source_inputs,
                      "pysam_version": pysam.__version__, "python_version": platform.python_version(),
                      "analysis_label": "EXPLORATORY_SOURCE_INPUT_ONLY"}
            with receipt.open("x") as f:
                json.dump(record, f, indent=2, sort_keys=True)
                f.write("\n")
            print(dataset, label, count, len(molecular_traits), flush=True)
            time.sleep(2)
        reader.close()


if __name__ == "__main__":
    main()
