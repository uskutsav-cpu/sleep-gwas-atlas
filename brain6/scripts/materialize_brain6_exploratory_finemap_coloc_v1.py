#!/usr/bin/env python3
"""Apply the frozen source/LD admission rule; emit allele-oriented RSS inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sqlite3
import statistics
from collections import defaultdict
from pathlib import Path

from audit_brain6_exploratory_finemap_coloc_v1 import orient, sha256, write_tsv


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "brain6/results/brain6_exploratory_finemap_coloc_v1"
SSD = Path("/Volumes/Extreme SSD/brain6-work")
WORK = SSD / "exploratory-finemap-coloc-v1"
REF = SSD / "lava-ukb-v1.1"
NORMALIZED = SSD / "preparation-v1"
CANDIDATE = REPO / "brain6/results/brain6_bounded_manuscript_v1/candidate_evidence_25.tsv"
MIN_SNP = 500
MIN_REF_COV = 0.50
MIN_GWAS_COV = 0.80
MIN_N_CONSISTENCY = 0.80
N_RATIO = 0.20
MIN_MAF = 0.01
MIN_NOBS = 90000


def read_candidates() -> list[dict]:
    source = list(csv.DictReader(CANDIDATE.open(), delimiter="\t"))
    if len(source) != 25 or len({r["candidate_locus_id"] for r in source}) != 25:
        raise RuntimeError("Protected candidate table must contain 25 unique loci")
    result = []
    for row in source:
        match = re.search(r"_chr(\d+)_(\d+)_(\d+)$", row["candidate_locus_id"])
        if match is None:
            raise RuntimeError(f"Invalid protected interval: {row['candidate_locus_id']}")
        chrom, start, stop = map(int, match.groups())
        result.append(dict(row, chrom=chrom, start=start, stop=stop))
    return result


def load_eligible_reference(candidates: list[dict]) -> dict[int, dict[str, dict]]:
    wanted = defaultdict(list)
    for loc in candidates:
        wanted[loc["chrom"]].append((loc["start"], loc["stop"]))
    result = defaultdict(dict)
    for chrom, intervals in wanted.items():
        path = REF / f"lava-ukb-v1.1_chr{chrom}.info"
        with path.open() as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            for row in rows:
                pos = int(row["POS"])
                if not any(a <= pos <= b for a, b in intervals):
                    continue
                freq = float(row["FREQ"])
                if min(freq, 1 - freq) < MIN_MAF or int(row["NOBS"]) < MIN_NOBS:
                    continue
                snp = row["SNP"].lower()
                if snp in result[chrom]:
                    raise RuntimeError(f"Duplicate reference rsID {snp}")
                result[chrom][snp] = row
    return result


def write_sumstats(path: Path, rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, rows, ["SNP", "CHR", "POS", "A1", "A2", "BETA", "SE", "Z", "N_EFF"])
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
    manifests = []
    data = {}
    for loc in candidates:
        cid, chrom, start, stop = (loc[k] for k in ("candidate_locus_id", "chrom", "start", "stop"))
        eligible_ref = {s: r for s, r in reference[chrom].items() if start <= int(r["POS"]) <= stop}
        for trait in loc["pair_id"].split("__"):
            source_rows = list(connections[trait].execute(
                "SELECT snp,chr,bp,a1,a2,beta,se,p,n FROM variants "
                "WHERE chr=? AND bp>=? AND bp<=?", (chrom, start, stop)
            ))
            valid = []
            for snp, chr0, pos, a1, a2, beta, se, p, n in source_rows:
                if not all(math.isfinite(v) for v in (beta, se, p, n)) or se <= 0 or n <= 0 or not 0 < p <= 1:
                    continue
                ref = eligible_ref.get(snp.lower())
                if ref is None or int(ref["POS"]) != pos:
                    continue
                orientation = orient(a1, a2, ref["A1"], ref["A2"])
                if orientation not in {"EXACT", "SWAP", "STRAND_COMPLEMENT", "STRAND_COMPLEMENT_SWAP"}:
                    continue
                flip = orientation in {"SWAP", "STRAND_COMPLEMENT_SWAP"}
                valid.append({
                    "SNP": snp.lower(), "CHR": chr0, "POS": pos,
                    "A1": ref["A1"], "A2": ref["A2"],
                    "BETA": -beta if flip else beta, "SE": se,
                    "Z": (-beta if flip else beta) / se, "N_EFF": n,
                })
            if len({r["SNP"] for r in valid}) != len(valid):
                raise RuntimeError(f"Duplicate oriented SNP in {cid}:{trait}")
            n_median = statistics.median(r["N_EFF"] for r in valid) if valid else 0
            consistent = [r for r in valid if (1 - N_RATIO) * n_median <= r["N_EFF"] <= (1 + N_RATIO) * n_median]
            consistent.sort(key=lambda r: (r["POS"], r["SNP"]))
            status = "READY"
            reason = ""
            if trait == "longsleep":
                status, reason = "NOT_RUN_MODEL_SEMANTICS", "BOLT_LMM_BINARY_TAIL_VARIANT_N_AND_MODEL_MAPPING_UNVERIFIED"
            elif len(consistent) < MIN_SNP:
                status, reason = "NOT_RUN_INPUT_QC", "LT_500_ORIENTED_N_CONSISTENT_SNPS"
            elif len(consistent) / len(eligible_ref) < MIN_REF_COV:
                status, reason = "NOT_RUN_INPUT_QC", "LT_50PCT_ELIGIBLE_LD_REFERENCE_COVERAGE"
            elif len(consistent) / len(source_rows) < MIN_GWAS_COV:
                status, reason = "NOT_RUN_INPUT_QC", "LT_80PCT_DENSE_GWAS_COVERAGE"
            elif len(consistent) / len(valid) < MIN_N_CONSISTENCY:
                status, reason = "NOT_RUN_INPUT_QC", "LT_80PCT_EFFECTIVE_N_CONSISTENCY"
            path = WORK / "inputs" / f"{cid}__{trait}.tsv"
            digest = write_sumstats(path, consistent) if status == "READY" else ""
            data[(cid, trait)] = consistent if status == "READY" else []
            manifests.append({
                "analysis_id": f"{cid}__{trait}", "candidate_locus_id": cid,
                "pair_id": loc["pair_id"], "trait_id": trait,
                "chromosome": chrom, "start": start, "stop": stop,
                "status": status, "reason": reason,
                "n_gwas_rows": len(source_rows), "n_ref_eligible": len(eligible_ref),
                "n_oriented": len(valid), "n_after_n_gate": len(consistent),
                "ref_coverage_fraction": round(len(consistent) / len(eligible_ref), 6) if eligible_ref else 0,
                "gwas_coverage_fraction": round(len(consistent) / len(source_rows), 6) if source_rows else 0,
                "n_consistency_fraction": round(len(consistent) / len(valid), 6) if valid else 0,
                "n_eff_median": n_median if valid else "",
                "input_tsv": str(path) if status == "READY" else "",
                "input_sha256": digest,
                "ld_rds": str(WORK / "ld" / f"{cid}__{trait}.rds") if status == "READY" else "",
            })
    pair_manifest = []
    by_analysis = {r["analysis_id"]: r for r in manifests}
    for loc in candidates:
        cid = loc["candidate_locus_id"]
        t1, t2 = loc["pair_id"].split("__")
        a = by_analysis[f"{cid}__{t1}"]
        b = by_analysis[f"{cid}__{t2}"]
        status, reason = "READY", ""
        pair_paths = ["", ""]
        hashes = ["", ""]
        shared = set()
        if a["status"] != "READY" or b["status"] != "READY":
            status = "NOT_RUN_MODEL_SEMANTICS" if "NOT_RUN_MODEL_SEMANTICS" in {a["status"], b["status"]} else "NOT_RUN_INPUT_QC"
            reason = f"{t1}:{a['status']}:{a['reason']};{t2}:{b['status']}:{b['reason']}"
        else:
            ia = {r["SNP"]: r for r in data[(cid, t1)]}
            ib = {r["SNP"]: r for r in data[(cid, t2)]}
            shared = set(ia) & set(ib)
            if len(shared) < MIN_SNP:
                status, reason = "NOT_RUN_INPUT_QC", "LT_500_SHARED_ORIENTED_SNPS"
            else:
                for j, (trait, dictionary) in enumerate(((t1, ia), (t2, ib))):
                    path = WORK / "inputs" / f"{cid}__shared__{trait}.tsv"
                    subset = [dictionary[s] for s in shared]
                    subset.sort(key=lambda r: (r["POS"], r["SNP"]))
                    hashes[j] = write_sumstats(path, subset)
                    pair_paths[j] = str(path)
        pair_manifest.append({
            "candidate_locus_id": cid, "pair_id": loc["pair_id"],
            "trait1": t1, "trait2": t2, "status": status, "reason": reason,
            "n_shared": len(shared), "trait1_input_tsv": pair_paths[0],
            "trait1_sha256": hashes[0], "trait2_input_tsv": pair_paths[1],
            "trait2_sha256": hashes[1],
            "ld_rds": str(WORK / "ld" / f"{cid}__shared.rds") if status == "READY" else "",
        })
    write_tsv(OUT / "fine_mapping_input_manifest.tsv", manifests, list(manifests[0]))
    write_tsv(OUT / "trait_coloc_input_manifest.tsv", pair_manifest, list(pair_manifest[0]))
    provenance = {
        "status": "PRE_OUTCOME_MATERIALIZATION",
        "candidate_sha256": sha256(CANDIDATE),
        "method_lock_sha256": sha256(OUT / "PRE_OUTCOME_METHOD_LOCK.md"),
        "script_sha256": sha256(Path(__file__)),
        "input_audit_provenance_sha256": sha256(OUT / "input_audit_provenance.json"),
        "output_sha256": {name: sha256(OUT / name) for name in (
            "fine_mapping_input_manifest.tsv", "trait_coloc_input_manifest.tsv")},
        "no_association_outcomes_inspected": True,
    }
    (OUT / "materialization_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    for connection in connections.values():
        connection.close()
    print(json.dumps({
        "fine_mapping_ready": sum(x["status"] == "READY" for x in manifests),
        "trait_coloc_ready": sum(x["status"] == "READY" for x in pair_manifest),
    }))


if __name__ == "__main__":
    main()
