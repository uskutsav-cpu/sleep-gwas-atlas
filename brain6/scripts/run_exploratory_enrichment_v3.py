#!/usr/bin/env python3
"""Feasibility-defined 18-region sensitivity using the frozen v2 background."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import defaultdict
from pathlib import Path

import build_exploratory_enrichment_v2 as v2

OUT = v2.ROOT / "brain6/results/brain6_exploratory_enrichment_v3"
V2 = v2.OUT
EXCLUDED = {"chr3:52219816-53219816", "chr17:43049526-45343136"}
SEED = 20260928
REPLICATES = 10000
MAX_ATTEMPTS = 200000


def read_tsv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sources() -> tuple[dict, list[dict], list[str], list[dict], dict]:
    data = json.loads((V2 / "prepared.json").read_text())
    audit = read_tsv(V2 / "matching_feasibility.tsv")
    tissues, scores, _ = v2.load_gtex(data["genes_by_id"])
    pathways, _ = v2.load_reactome(data["genes_by_id"], set(data["universe_gene_ids"]))
    assert len(audit) == 20 and len(tissues) == 54 and len(pathways) == 1680
    assert {r["region_group"] for r in audit if int(r["eligible_null_anchors"]) < 10} == EXCLUDED
    return data, audit, tissues, pathways, scores


def freeze() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data, audit, tissues, pathways, _ = sources()
    included = [r for r in audit if r["region_group"] not in EXCLUDED]
    assert len(included) == 18 and all(int(r["eligible_null_anchors"]) >= 10 for r in included)
    set_rows = [{"region_group": r["region_group"], "v2_eligible_null_anchors": r["eligible_null_anchors"],
                 "subset_reason": "PRESPECIFIED_V2_MATCHING_GATE_PASSED"} for r in included]
    v2.write_tsv(OUT / "FROZEN_SET_18.tsv", set_rows, list(set_rows[0]))
    tissue_rows = [{"tissue": name} for name in tissues]
    v2.write_tsv(OUT / "FROZEN_TISSUES_54.tsv", tissue_rows, ["tissue"])
    pathway_rows = [{"reactome_id": p["reactome_id"], "pathway": p["name"],
                     "mapped_tested_universe_genes": len(p["genes"])} for p in pathways]
    v2.write_tsv(OUT / "FROZEN_PATHWAYS_1680.tsv", pathway_rows, list(pathway_rows[0]))
    status_rows = [{"region_group": r["region_group"], "pair_candidates": next(
                       x["n_pair_specific_candidates"] for x in v2.load_regions() if x["region_group"] == r["region_group"]),
                    "v2_eligible_null_anchors": r["eligible_null_anchors"],
                    "v3_subset": "INCLUDED_18" if r["region_group"] not in EXCLUDED else "EXCLUDED_MATCHING_GATE",
                    "v3_reason": "v2 >=10 null anchors" if r["region_group"] not in EXCLUDED else "v2 <10 null anchors"}
                   for r in audit]
    v2.write_tsv(OUT / "matching_status_20.tsv", status_rows, list(status_rows[0]))
    binding = {"v2_protocol_commit": "cfb19bc1", "v2_full_set_audit_commit": "bbdead78",
               "v2_prepared_sha256": file_sha(V2 / "prepared.json"),
               "v2_matching_sha256": file_sha(V2 / "matching_feasibility.tsv"),
               "v2_source_receipts_sha256": file_sha(V2 / "source_receipts.tsv"),
               "18_region_ids": [x["region_group"] for x in set_rows],
               "excluded_region_ids": sorted(EXCLUDED), "seed": SEED,
               "target_accepted_replicates": REPLICATES, "maximum_attempts": MAX_ATTEMPTS,
               "source_manifest_sha256": {p.name: file_sha(p) for p in
                   [OUT / "FROZEN_SET_18.tsv", OUT / "FROZEN_TISSUES_54.tsv", OUT / "FROZEN_PATHWAYS_1680.tsv"]}}
    (OUT / "frozen_binding.json").write_text(json.dumps(binding, indent=2) + "\n")
    print(json.dumps({"region_count": len(set_rows), "tissue_tests": len(tissues),
                      "pathway_tests": len(pathways), "excluded": sorted(EXCLUDED)}, indent=2))


def bh(pvals: list[float]) -> list[float]:
    n = len(pvals)
    out = [1.0] * n
    order = sorted(range(n), key=lambda i: pvals[i])
    least = 1.0
    for j in range(n - 1, -1, -1):
        i = order[j]
        least = min(least, pvals[i] * n / (j + 1))
        out[i] = least
    return out


def pairwise_nonoverlap(rows: list[dict]) -> bool:
    for i, a in enumerate(rows):
        if any(v2.intervals_overlap(a, b) for b in rows[i + 1:]):
            return False
    return True


def enrich() -> None:
    data, audit, tissues, pathways, scores = sources()
    binding = json.loads((OUT / "frozen_binding.json").read_text())
    receipts = {r["input"]: r["sha256"] for r in read_tsv(V2 / "source_receipts.tsv")}
    for source in (v2.GTF, v2.GTEX, v2.REACTOME):
        assert file_sha(source) == receipts[str(source)]
    assert file_sha(V2 / "prepared.json") == binding["v2_prepared_sha256"]
    assert file_sha(V2 / "matching_feasibility.tsv") == binding["v2_matching_sha256"]
    assert file_sha(V2 / "source_receipts.tsv") == binding["v2_source_receipts_sha256"]
    for name, digest in binding["source_manifest_sha256"].items():
        assert file_sha(OUT / name) == digest
    assert [x["region_group"] for x in read_tsv(OUT / "FROZEN_SET_18.tsv")] == binding["18_region_ids"]
    assert [x["tissue"] for x in read_tsv(OUT / "FROZEN_TISSUES_54.tsv")] == tissues
    assert [x["reactome_id"] for x in read_tsv(OUT / "FROZEN_PATHWAYS_1680.tsv")] == [p["reactome_id"] for p in pathways]
    targets = [t for t in data["targets"] if t["region_group"] in binding["18_region_ids"]]
    assert len(targets) == 18
    pools: dict[str, list[dict]] = defaultdict(list)
    for p in data["pools"]:
        if p["region_group"] in binding["18_region_ids"]:
            pools[p["region_group"]].append(p)
    assert all(len(pools[t["region_group"]]) >= 10 for t in targets)
    gene_to_pathways: dict[str, list[int]] = defaultdict(list)
    for j, p in enumerate(pathways):
        for gene in p["genes"]:
            gene_to_pathways[gene].append(j)

    def pathway_hits(genes: list[str]) -> set[int]:
        return {j for gene in genes for j in gene_to_pathways.get(gene, ())}

    def tissue_score(genes: list[str]) -> list[float]:
        return v2.region_tissue(genes, scores, 54)[0]

    observed_tissue = [0.0] * 54
    observed_pathway = [0] * len(pathways)
    for t in targets:
        for i, value in enumerate(tissue_score(t["genes"])):
            observed_tissue[i] += value
        for j in pathway_hits(t["genes"]):
            observed_pathway[j] += 1
    prepared_nulls = {}
    for group, rows in pools.items():
        for n in rows:
            key = (group, n["anchor_loc"])
            prepared_nulls[key] = (tissue_score(n["genes"]), pathway_hits(n["genes"]))
    extreme_tissue = [0] * 54
    extreme_pathway = [0] * len(pathways)
    rng = random.Random(SEED)
    accepted = attempts = 0
    while accepted < REPLICATES and attempts < MAX_ATTEMPTS:
        attempts += 1
        nulls = [rng.choice(pools[t["region_group"]]) for t in targets]
        if not pairwise_nonoverlap(nulls):
            continue
        accepted += 1
        null_tissue = [0.0] * 54
        null_pathway = [0] * len(pathways)
        for n in nulls:
            ts, ps = prepared_nulls[(n["region_group"], n["anchor_loc"])]
            for i, value in enumerate(ts):
                null_tissue[i] += value
            for j in ps:
                null_pathway[j] += 1
        for i, value in enumerate(null_tissue):
            extreme_tissue[i] += value >= observed_tissue[i]
        for j, value in enumerate(null_pathway):
            extreme_pathway[j] += value >= observed_pathway[j]
    run = {"status": "PASS" if accepted == REPLICATES else "JOINT_NULL_GATE_FAILED",
           "subset_regions": 18, "excluded_full_set_regions": sorted(EXCLUDED),
           "accepted_null_replicates": accepted, "attempts": attempts,
           "seed": SEED, "tissue_tests": 54, "pathway_tests": len(pathways),
           "interpretation": "exploratory positional 18-region sensitivity; no full-20 enrichment result"}
    (OUT / "run_summary.json").write_text(json.dumps(run, indent=2) + "\n")
    if accepted < REPLICATES:
        print(json.dumps(run, indent=2))
        return
    tp = [(1 + x) / (1 + accepted) for x in extreme_tissue]
    tq = bh(tp)
    tissue_rows = [{"tissue": tissues[i], "region_subset": "MATCHED_18_OF_FROZEN_20",
                    "candidate_sum_region_mean_log2_tpm_specificity": f"{observed_tissue[i]:.12g}",
                    "null_extreme_count": extreme_tissue[i], "p_empirical_high": f"{tp[i]:.12g}",
                    "q_bh_54": f"{tq[i]:.12g}"} for i in range(54)]
    v2.write_tsv(OUT / "tissue_bulk_gtex_v8_enrichment_18.tsv", tissue_rows, list(tissue_rows[0]))
    pp = [(1 + x) / (1 + accepted) for x in extreme_pathway]
    pq = bh(pp)
    pathway_rows = [{"reactome_id": p["reactome_id"], "pathway": p["name"],
                     "region_subset": "MATCHED_18_OF_FROZEN_20",
                     "mapped_tested_universe_genes": len(p["genes"]),
                     "candidate_regions_with_pathway_gene": observed_pathway[j],
                     "null_extreme_count": extreme_pathway[j],
                     "p_empirical_high": f"{pp[j]:.12g}", "q_bh_1680": f"{pq[j]:.12g}"}
                    for j, p in enumerate(pathways)]
    v2.write_tsv(OUT / "reactome_positional_enrichment_18.tsv", pathway_rows, list(pathway_rows[0]))
    print(json.dumps(run, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("freeze", "enrich"))
    args = parser.parse_args()
    {"freeze": freeze, "enrich": enrich}[args.phase]()
