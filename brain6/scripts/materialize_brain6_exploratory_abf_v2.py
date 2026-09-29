#!/usr/bin/env python3
"""Materialize prelocked, allele-oriented binary GWAS ABF locus inputs."""

from __future__ import annotations

import csv
import json
import math
import sqlite3
from pathlib import Path

from audit_brain6_exploratory_finemap_coloc_v1 import orient, sha256, write_tsv
from materialize_brain6_exploratory_finemap_coloc_v1 import (
    CANDIDATE, NORMALIZED, read_candidates, load_eligible_reference,
)


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "brain6/results/brain6_exploratory_finemap_coloc_v2"
WORK = Path("/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2")
MIN_SNP = 500
MIN_REF_COV = 0.50
MIN_GWAS_COV = 0.80


def write_input(path: Path, rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, rows, ["SNP", "CHR", "POS", "A1", "A2", "BETA", "SE", "N_AUDIT"])
    return sha256(path)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (WORK / "inputs").mkdir(parents=True, exist_ok=True)
    candidates = read_candidates()
    reference = load_eligible_reference(candidates)
    connections = {}
    for trait in {t for c in candidates for t in c["pair_id"].split("__")}:
        path = NORMALIZED / f"normalize_{trait}/variants.sqlite"
        connections[trait] = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    trait_manifest = []
    eligible_data = {}
    for loc in candidates:
        cid = loc["candidate_locus_id"]
        chrom, start, stop = loc["chrom"], loc["start"], loc["stop"]
        ref = {s: r for s, r in reference[chrom].items() if start <= int(r["POS"]) <= stop}
        for trait in loc["pair_id"].split("__"):
            raw = list(connections[trait].execute(
                "SELECT snp,chr,bp,a1,a2,beta,se,p,n FROM variants WHERE chr=? AND bp>=? AND bp<=?",
                (chrom, start, stop),
            ))
            valid_source = 0
            matched = []
            for snp, chr0, pos, a1, a2, beta, se, p, n in raw:
                if not all(math.isfinite(x) for x in (beta, se, p, n)) or se <= 0 or n <= 0 or not (0 < p <= 1):
                    continue
                valid_source += 1
                marker = ref.get(snp.lower())
                if marker is None or int(marker["POS"]) != pos:
                    continue
                orientation = orient(a1, a2, marker["A1"], marker["A2"])
                if orientation not in {"EXACT", "SWAP", "STRAND_COMPLEMENT", "STRAND_COMPLEMENT_SWAP"}:
                    continue
                effect = -beta if orientation in {"SWAP", "STRAND_COMPLEMENT_SWAP"} else beta
                matched.append({
                    "SNP": snp.lower(), "CHR": chr0, "POS": pos,
                    "A1": marker["A1"], "A2": marker["A2"],
                    "BETA": effect, "SE": se, "N_AUDIT": n,
                })
            if len({x["SNP"] for x in matched}) != len(matched):
                raise RuntimeError(f"Duplicate oriented SNP: {cid}:{trait}")
            matched.sort(key=lambda x: (x["POS"], x["SNP"]))
            status, reason = "READY", ""
            if trait == "longsleep":
                status, reason = "NOT_RUN_MODEL_SEMANTICS", "LONGSLEEP_BOLT_LMM_EFFECT_SCALE_UNVERIFIED_FOR_BINARY_ABF"
            elif len(matched) < MIN_SNP:
                status, reason = "NOT_RUN_INPUT_QC", "LT_500_ORIENTED_SNPS"
            elif len(matched) / len(ref) < MIN_REF_COV:
                status, reason = "NOT_RUN_INPUT_QC", "LT_50PCT_ELIGIBLE_REFERENCE_COVERAGE"
            elif len(matched) / valid_source < MIN_GWAS_COV:
                status, reason = "NOT_RUN_INPUT_QC", "LT_80PCT_VALID_GWAS_ROW_COVERAGE"
            path = WORK / "inputs" / f"{cid}__{trait}.tsv"
            digest = write_input(path, matched) if status == "READY" else ""
            eligible_data[(cid, trait)] = {x["SNP"]: x for x in matched} if status == "READY" else {}
            trait_manifest.append({
                "analysis_id": f"{cid}__{trait}", "candidate_locus_id": cid,
                "pair_id": loc["pair_id"], "trait_id": trait,
                "chromosome": chrom, "start": start, "stop": stop,
                "status": status, "reason": reason,
                "n_gwas_rows": len(raw), "n_valid_gwas_rows": valid_source,
                "n_reference_eligible": len(ref), "n_oriented": len(matched),
                "reference_coverage_fraction": round(len(matched)/len(ref), 6) if ref else 0,
                "gwas_coverage_fraction": round(len(matched)/valid_source, 6) if valid_source else 0,
                "input_tsv": str(path) if status == "READY" else "",
                "input_sha256": digest,
            })
    by_id = {x["analysis_id"]: x for x in trait_manifest}
    pair_manifest = []
    for loc in candidates:
        cid = loc["candidate_locus_id"]
        t1, t2 = loc["pair_id"].split("__")
        m1, m2 = (by_id[f"{cid}__{t}"] for t in (t1, t2))
        status, reason = "READY", ""
        shared = set()
        paths = ["", ""]
        digests = ["", ""]
        if m1["status"] != "READY" or m2["status"] != "READY":
            status = "NOT_RUN_MODEL_SEMANTICS" if "NOT_RUN_MODEL_SEMANTICS" in {m1["status"], m2["status"]} else "NOT_RUN_INPUT_QC"
            reason = f"{t1}:{m1['status']}:{m1['reason']};{t2}:{m2['status']}:{m2['reason']}"
        else:
            d1, d2 = eligible_data[(cid,t1)], eligible_data[(cid,t2)]
            shared = set(d1) & set(d2)
            if len(shared) < MIN_SNP:
                status, reason = "NOT_RUN_INPUT_QC", "LT_500_SHARED_SNPS"
            else:
                for i, (trait, data) in enumerate(((t1,d1),(t2,d2))):
                    selected = [data[s] for s in shared]
                    selected.sort(key=lambda x: (x["POS"],x["SNP"]))
                    path = WORK / "inputs" / f"{cid}__shared__{trait}.tsv"
                    digests[i] = write_input(path, selected)
                    paths[i] = str(path)
        pair_manifest.append({
            "candidate_locus_id": cid, "pair_id": loc["pair_id"],
            "trait1": t1, "trait2": t2, "status": status, "reason": reason,
            "n_shared": len(shared), "trait1_input_tsv": paths[0],
            "trait1_sha256": digests[0], "trait2_input_tsv": paths[1],
            "trait2_sha256": digests[1],
        })
    write_tsv(OUT / "fine_mapping_input_manifest.tsv", trait_manifest, list(trait_manifest[0]))
    write_tsv(OUT / "trait_coloc_input_manifest.tsv", pair_manifest, list(pair_manifest[0]))
    (OUT / "materialization_provenance.json").write_text(json.dumps({
        "status": "PRE_OUTCOME_INPUT_MATERIALIZATION",
        "candidate_sha256": sha256(CANDIDATE),
        "method_lock_sha256": sha256(OUT / "PRE_OUTCOME_ABF_METHOD_LOCK.md"),
        "script_sha256": sha256(Path(__file__)),
        "fine_mapping_manifest_sha256": sha256(OUT / "fine_mapping_input_manifest.tsv"),
        "trait_coloc_manifest_sha256": sha256(OUT / "trait_coloc_input_manifest.tsv"),
    }, indent=2, sort_keys=True) + "\n")
    for connection in connections.values():
        connection.close()
    print(json.dumps({"fine_mapping_ready": sum(x["status"] == "READY" for x in trait_manifest),
                      "trait_coloc_ready": sum(x["status"] == "READY" for x in pair_manifest)}))


if __name__ == "__main__":
    main()
