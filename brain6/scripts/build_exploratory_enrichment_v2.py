#!/usr/bin/env python3
"""Frozen-protocol positional enrichment for the 20 Brain6 PLACO regions.

Run ``prepare`` before ``enrich``. This script never reads PLACO P/Q columns.
Its outputs are exploratory annotations, not locus or gene validation.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import random
import re
import statistics
import zipfile
from array import array
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/brain6_exploratory_enrichment_v2"
EXT = Path("/Volumes/Extreme SSD/brain6-work/brain6-exploratory-enrichment-v2/inputs")
REGIONS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/geographic_region_groups.tsv"
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
BLOCKS = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
GTF = EXT / "gencode.v19.annotation.gtf.gz"
GTEX = EXT / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
REACTOME = EXT / "ReactomePathways.gmt.zip"
PAIRS = ("insomnia__adhd", "insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz")
PLACOS = {p: ROOT / f"brain6/results/placo/{p}/variants.tsv.gz" for p in PAIRS}
SOURCE_URLS = {
    GTF: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
    GTEX: "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    REACTOME: "https://reactome.org/download/current/ReactomePathways.gmt.zip",
    REGIONS: "repository:frozen-five-track-region-reconciliation-v1",
    MEMBERS: "repository:frozen-25-pair-candidate-members-v1",
    BLOCKS: "repository:pinned-LAVA-GRCh37-2495-blocks",
    **{v: "repository:complete-five-track-PLACO:" + k for k, v in PLACOS.items()},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def load_regions() -> list[dict]:
    with REGIONS.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 20
    for r in rows:
        r["chrom"] = int(r["chromosome"])
        r["start"] = int(r["start"])
        r["stop"] = int(r["stop"])
        r["pairs_list"] = r["pairs"].split(";")
    assert sum(int(r["n_pair_specific_candidates"]) for r in rows) == 25
    return rows


def load_blocks() -> dict[int, list[dict]]:
    by_chr: dict[int, list[dict]] = defaultdict(list)
    with BLOCKS.open() as f:
        next(f)
        for line in f:
            loc, ch, s, e = map(int, line.split())
            if 1 <= ch <= 22:
                by_chr[ch].append({"loc": loc, "start": s, "stop": e, "mid": (s + e) // 2, "length": e - s + 1})
    assert sum(map(len, by_chr.values())) == 2495
    return by_chr


ATTR = re.compile(r'([A-Za-z_]+) "([^"]+)"')


def load_genes() -> tuple[dict[int, list[dict]], dict[str, dict]]:
    by_chr: dict[int, list[dict]] = defaultdict(list)
    by_id: dict[str, dict] = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) < 9 or col[2] != "gene":
                continue
            chrom_text = col[0].removeprefix("chr")
            if not chrom_text.isdigit():
                continue
            ch = int(chrom_text)
            if not 1 <= ch <= 22:
                continue
            at = dict(ATTR.findall(col[8]))
            if at.get("gene_type", at.get("gene_biotype")) != "protein_coding":
                continue
            gid = at["gene_id"].split(".")[0]
            gene = {"id": gid, "symbol": at.get("gene_name", ""), "chrom": ch,
                    "start": int(col[3]), "stop": int(col[4])}
            assert gid not in by_id
            by_id[gid] = gene
            by_chr[ch].append(gene)
    for genes in by_chr.values():
        genes.sort(key=lambda x: x["start"])
    return by_chr, by_id


def gene_ids(genes: list[dict], left: int, right: int) -> list[str]:
    # The GTF gene records are sorted by start; exact body overlap is inclusive.
    return [g["id"] for g in genes if g["start"] <= right and g["stop"] >= left]


def load_tested_positions() -> tuple[dict[str, dict[int, array]], dict[str, int]]:
    positions: dict[str, dict[int, array]] = {}
    counts: dict[str, int] = {}
    for pair, path in PLACOS.items():
        by_chr: dict[int, array] = defaultdict(lambda: array("I"))
        n = 0
        with gzip.open(path, "rt") as f:
            header = next(f).rstrip("\n").split("\t")
            ci, pi, si = (header.index(x) for x in ("CHR", "BP", "status"))
            for line in f:
                c = line.rstrip("\n").split("\t")
                if c[si] != "TESTED":
                    continue
                try:
                    ch, pos = int(c[ci]), int(c[pi])
                except ValueError:
                    continue
                if 1 <= ch <= 22 and pos > 0:
                    by_chr[ch].append(pos)
                    n += 1
        for ch, a in by_chr.items():
            if any(a[i] > a[i + 1] for i in range(len(a) - 1)):
                by_chr[ch] = array("I", sorted(a))
        positions[pair] = by_chr
        counts[pair] = n
    return positions, counts


def npos(pos: array, left: int, right: int) -> int:
    return bisect.bisect_right(pos, right) - bisect.bisect_left(pos, left)


def interval_tested(positions: dict, pairs: list[str], ch: int, left: int, right: int) -> dict[str, int]:
    return {p: npos(positions[p].get(ch, array("I")), left, right) for p in pairs}


def midpoint_block(blocks: list[dict], point: int) -> dict | None:
    # Blocks are non-overlapping and sorted by chromosome position.
    starts = [b["start"] for b in blocks]
    j = bisect.bisect_right(starts, point) - 1
    return blocks[j] if j >= 0 and point <= blocks[j]["stop"] else None


def intervals_overlap(a: dict, b: dict) -> bool:
    return a["chrom"] == b["chrom"] and a["start"] <= b["stop"] and b["start"] <= a["stop"]


def prepare() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    regions = load_regions()
    blocks = load_blocks()
    genes, genes_by_id = load_genes()
    positions, tested_counts = load_tested_positions()
    target_rows = []
    pool_rows = []
    audits = []
    all_blocks = [b | {"chrom": ch} for ch, bs in blocks.items() for b in bs]
    universe: set[str] = set()
    tested_block_count = 0
    for b in all_blocks:
        if sum(interval_tested(positions, list(PAIRS), b["chrom"], b["start"], b["stop"]).values()):
            tested_block_count += 1
            universe.update(gene_ids(genes.get(b["chrom"], []), b["start"], b["stop"]))
    for r in regions:
        ch, left, right = r["chrom"], r["start"], r["stop"]
        target_genes = gene_ids(genes.get(ch, []), left, right)
        target_tested = interval_tested(positions, r["pairs_list"], ch, left, right)
        anchor = midpoint_block(blocks[ch], (left + right) // 2)
        target = {"region_group": r["region_group"], "chrom": ch, "start": left, "stop": right,
                  "pairs": r["pairs_list"], "genes": target_genes, "tested": target_tested,
                  "anchor_loc": anchor["loc"] if anchor else None,
                  "anchor_length": anchor["length"] if anchor else None}
        target_rows.append(target)
        excluded = defaultdict(int)
        eligible = []
        width = right - left + 1
        chromosome_stop = max(b["stop"] for b in blocks[ch])
        for b in blocks[ch]:
            s = b["mid"] - (width - 1) // 2
            e = s + width - 1
            n = {"region_group": r["region_group"], "chrom": ch, "start": s, "stop": e,
                 "anchor_loc": b["loc"], "anchor_length": b["length"]}
            if s < 1 or e > chromosome_stop:
                excluded["chromosome_edge"] += 1
                continue
            if any(intervals_overlap(n, q) for q in regions):
                excluded["candidate_overlap"] += 1
                continue
            if anchor is None:
                excluded["target_without_anchor"] += 1
                continue
            if not 0.5 <= b["length"] / anchor["length"] <= 2.0:
                excluded["ld_block_span"] += 1
                continue
            null_tested = interval_tested(positions, r["pairs_list"], ch, s, e)
            if any((null_tested[p] != 0 if target_tested[p] == 0 else
                    not 0.5 <= null_tested[p] / target_tested[p] <= 2.0)
                   for p in r["pairs_list"]):
                excluded["tested_density"] += 1
                continue
            null_genes = gene_ids(genes.get(ch, []), s, e)
            tg, ng = len(target_genes), len(null_genes)
            if (ng != 0 if tg == 0 else abs(ng - tg) > max(2, math.ceil(0.5 * tg))):
                excluded["gene_count"] += 1
                continue
            n["genes"] = null_genes
            n["tested"] = null_tested
            eligible.append(n)
            pool_rows.append(n)
        audits.append({"region_group": r["region_group"], "chromosome": ch, "width_bp": width,
                       "pairs": ";".join(r["pairs_list"]), "target_genes": len(target_genes),
                       "target_tested": json.dumps(target_tested, sort_keys=True),
                       "target_anchor_loc": target["anchor_loc"], "target_anchor_length": target["anchor_length"],
                       "same_chrom_anchors": len(blocks[ch]), "eligible_null_anchors": len(eligible),
                       **{f"excluded_{k}": excluded[k] for k in
                          ("chromosome_edge", "candidate_overlap", "target_without_anchor", "ld_block_span", "tested_density", "gene_count")}})
    payload = {"protocol_commit": "cfb19bc1", "targets": target_rows, "pools": pool_rows,
               "universe_gene_ids": sorted(universe), "tested_block_count": tested_block_count,
               "genes_by_id": genes_by_id, "total_protein_coding_genes": len(genes_by_id),
               "tested_variant_counts": tested_counts}
    (OUT / "prepared.json").write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    fields = list(audits[0])
    write_tsv(OUT / "matching_feasibility.tsv", audits, fields)
    receipts = [{"input": str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),
                 "source_url": url, "bytes": p.stat().st_size, "sha256": sha256(p)}
                for p, url in SOURCE_URLS.items()]
    write_tsv(OUT / "source_receipts.tsv", receipts, ["input", "source_url", "bytes", "sha256"])
    summary = {"targets": len(target_rows), "pair_candidates": 25, "lava_blocks": len(all_blocks),
               "tested_blocks": tested_block_count, "protein_coding_genes": len(genes_by_id),
               "tested_universe_genes": len(universe), "tested_variant_counts": tested_counts,
               "eligible_pool_min": min(a["eligible_null_anchors"] for a in audits),
               "regions_below_10": [a["region_group"] for a in audits if a["eligible_null_anchors"] < 10]}
    (OUT / "preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    n = len(p_values)
    q = [0.0] * n
    running = 1.0
    for i in sorted(range(n), key=lambda j: p_values[j], reverse=True):
        rank = sum(v <= p_values[i] for v in p_values)
        running = min(running, p_values[i] * n / rank)
        q[i] = running
    return q


def load_gtex(genes: dict[str, dict]) -> tuple[list[str], dict[str, list[float]], int]:
    scores = {}
    matched = 0
    with gzip.open(GTEX, "rt") as f:
        next(f)
        next(f)
        header = next(f).rstrip("\n").split("\t")
        tissues = header[2:]
        assert len(tissues) == 54
        for line in f:
            cols = line.rstrip("\n").split("\t")
            gid = cols[0].split(".")[0]
            if gid not in genes:
                continue
            values = [math.log2(1.0 + float(x)) for x in cols[2:]]
            avg = statistics.fmean(values)
            scores[gid] = [x - avg for x in values]
            matched += 1
    return tissues, scores, matched


def load_reactome(genes: dict[str, dict], universe: set[str]) -> tuple[list[dict], dict]:
    symbols: dict[str, list[str]] = defaultdict(list)
    for gid, g in genes.items():
        if g["symbol"]:
            symbols[g["symbol"]].append(gid)
    pathways = []
    counts = defaultdict(int)
    with zipfile.ZipFile(REACTOME) as z:
        with z.open("ReactomePathways.gmt") as f:
            for raw in f:
                cols = raw.decode("utf-8").rstrip("\n").split("\t")
                counts["gmt_pathways"] += 1
                ids = {symbols[s][0] for s in cols[2:] if len(symbols[s]) == 1}
                ids &= universe
                if 10 <= len(ids) <= 500:
                    pathways.append({"name": cols[0], "reactome_id": cols[1], "genes": sorted(ids)})
                else:
                    counts["outside_size_bounds"] += 1
    return pathways, dict(counts)


def region_tissue(genes: list[str], scores: dict[str, list[float]], n_tissues: int) -> tuple[list[float], int]:
    present = [scores[g] for g in genes if g in scores]
    if not present:
        return [0.0] * n_tissues, 0
    return [sum(row[i] for row in present) / len(present) for i in range(n_tissues)], len(present)


def select_nulls(targets: list[dict], pools: dict[str, list[dict]], rng: random.Random) -> list[dict] | None:
    selected = []
    for target in sorted(targets, key=lambda t: len(pools[t["region_group"]])):
        options = pools[target["region_group"]]
        available = [x for x in options if not any(intervals_overlap(x, q) for q in selected)]
        if not available:
            return None
        selected.append(rng.choice(available))
    return selected


def enrich() -> None:
    data = json.loads((OUT / "prepared.json").read_text())
    targets = data["targets"]
    pools: dict[str, list[dict]] = defaultdict(list)
    for row in data["pools"]:
        pools[row["region_group"]].append(row)
    genes = data["genes_by_id"]
    tissues, scores, gtex_matched = load_gtex(genes)
    pathways, pathway_counts = load_reactome(genes, set(data["universe_gene_ids"]))
    coverage_rows = [{"region_group": t["region_group"], "protein_coding_gene_bodies": len(t["genes"]),
                      "gtex_v8_mapped_gene_ids": sum(g in scores for g in t["genes"]),
                      "eligible_null_anchors": len(pools[t["region_group"]])} for t in targets]
    write_tsv(OUT / "positional_annotation_coverage_20.tsv", coverage_rows, list(coverage_rows[0]))
    coverage = {"candidate_regions": len(targets),
                "candidate_union_protein_coding_gene_ids": len(set(g for t in targets for g in t["genes"])),
                "candidate_regions_with_at_least_one_gene": sum(bool(t["genes"]) for t in targets),
                "candidate_regions_with_gtex_gene": sum(x["gtex_v8_mapped_gene_ids"] > 0 for x in coverage_rows),
                "gtex_tissues": len(tissues), "gtex_genes_matching_gencode": gtex_matched,
                "reactome_pathways_in_source": pathway_counts["gmt_pathways"],
                "reactome_pathways_size_eligible": len(pathways),
                "reactome_pathways_outside_size_bounds": pathway_counts["outside_size_bounds"]}
    (OUT / "annotation_feasibility.json").write_text(json.dumps(coverage, indent=2) + "\n")
    failed = [t["region_group"] for t in targets if len(pools[t["region_group"]]) < 10]
    if failed:
        (OUT / "ENRICHMENT_STATUS.md").write_text(
            "# Prespecified matching gate failed\n\nNo enrichment p-values were calculated. "
            f"{len(failed)} of 20 candidate regions had fewer than 10 eligible null anchors: "
            + ", ".join(failed) + ". See `matching_feasibility.tsv`.\n")
        print("MATCHING_GATE_FAILED", len(failed), failed)
        return
    pathway_sets = [set(x["genes"]) for x in pathways]
    target_tissue = [region_tissue(t["genes"], scores, len(tissues)) for t in targets]
    obs_tissue = [sum(x[0][i] for x in target_tissue) for i in range(len(tissues))]
    obs_pathway = [sum(bool(set(t["genes"]) & ps) for t in targets) for ps in pathway_sets]
    tissue_extreme = [0] * len(tissues)
    pathway_extreme = [0] * len(pathways)
    rng = random.Random(20260928)
    accepted = 0
    attempts = 0
    while accepted < 10000 and attempts < 200000:
        attempts += 1
        nulls = select_nulls(targets, pools, rng)
        if nulls is None:
            continue
        accepted += 1
        nt = [region_tissue(n["genes"], scores, len(tissues))[0] for n in nulls]
        for i in range(len(tissues)):
            tissue_extreme[i] += sum(x[i] for x in nt) >= obs_tissue[i] - 1e-12
        null_gene_sets = [set(n["genes"]) for n in nulls]
        for i, ps in enumerate(pathway_sets):
            pathway_extreme[i] += sum(bool(gs & ps) for gs in null_gene_sets) >= obs_pathway[i]
    status = {"accepted_null_replicates": accepted, "attempts": attempts,
              "seed": 20260928, "gtex_genes_matching_gencode": gtex_matched,
              "reactome": pathway_counts, "reactome_pathways_tested": len(pathways),
              "candidate_regions_with_gtex_genes": sum(x[1] > 0 for x in target_tissue),
              "candidate_genes_with_gtex": len(set(g for t in targets for g in t["genes"]) & set(scores))}
    (OUT / "enrichment_run_summary.json").write_text(json.dumps(status, indent=2) + "\n")
    if accepted < 10000:
        (OUT / "ENRICHMENT_STATUS.md").write_text(
            "# Prespecified joint-null gate failed\n\nNo enrichment p-values were calculated. "
            f"Only {accepted} accepted non-overlapping null replicates in {attempts} attempts.\n")
        print(json.dumps(status, indent=2))
        return
    tp = [(x + 1) / (accepted + 1) for x in tissue_extreme]
    tq = benjamini_hochberg(tp)
    tissue_rows = [{"tissue": tissues[i], "candidate_sum_region_mean_log2_tpm_specificity": f"{obs_tissue[i]:.9g}",
                    "candidate_regions_with_gtex_genes": status["candidate_regions_with_gtex_genes"],
                    "null_extreme_count": tissue_extreme[i], "p_empirical_high": f"{tp[i]:.9g}",
                    "q_bh_54": f"{tq[i]:.9g}"} for i in range(len(tissues))]
    write_tsv(OUT / "tissue_bulk_gtex_v8_enrichment.tsv", tissue_rows, list(tissue_rows[0]))
    pp = [(x + 1) / (accepted + 1) for x in pathway_extreme]
    pq = benjamini_hochberg(pp)
    pathway_rows = [{"reactome_id": pw["reactome_id"], "pathway": pw["name"],
                     "mapped_universe_genes": len(pw["genes"]), "candidate_regions_with_pathway_gene": obs_pathway[i],
                     "null_extreme_count": pathway_extreme[i], "p_empirical_high": f"{pp[i]:.9g}",
                     "q_bh_all_tested": f"{pq[i]:.9g}"} for i, pw in enumerate(pathways)]
    write_tsv(OUT / "reactome_positional_enrichment.tsv", pathway_rows,
              list(pathway_rows[0]) if pathway_rows else
              ["reactome_id", "pathway", "mapped_universe_genes", "candidate_regions_with_pathway_gene", "null_extreme_count", "p_empirical_high", "q_bh_all_tested"])
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("prepare", "enrich"))
    args = parser.parse_args()
    {"prepare": prepare, "enrich": enrich}[args.phase]()
