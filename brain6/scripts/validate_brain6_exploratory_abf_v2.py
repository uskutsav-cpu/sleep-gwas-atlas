#!/usr/bin/env python3
"""Independent family-wide integrity audit of the frozen single-signal ABF branch."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from audit_brain6_exploratory_finemap_coloc_v1 import sha256


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "brain6/results/brain6_exploratory_finemap_coloc_v2"
V1 = REPO / "brain6/results/brain6_exploratory_finemap_coloc_v1"
SSD = Path("/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2")
CANDIDATES = REPO / "brain6/results/brain6_bounded_manuscript_v1/candidate_evidence_25.tsv"
LIB = Path("/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib")


def read(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def approx(a: float, b: float, tolerance: float = 1e-6) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tolerance


def logsumexp(values: list[float]) -> float:
    high = max(values)
    return high + math.log(math.fsum(math.exp(x - high) for x in values))


def log_bf(beta: float, se: float) -> float:
    variance = se * se
    ratio = 0.2**2 / (0.2**2 + variance)
    return 0.5 * (math.log1p(-ratio) + ratio * (beta / se) ** 2)


def expected_coloc(a: list[dict[str,str]], b: list[dict[str,str]], p12: float) -> list[float]:
    l1 = [log_bf(float(x["BETA"]),float(x["SE"])) for x in a]
    l2 = [log_bf(float(x["BETA"]),float(x["SE"])) for x in b]
    lsum1, lsum2 = logsumexp(l1), logsumexp(l2)
    lsum12 = logsumexp([x+y for x,y in zip(l1,l2)])
    relation = math.exp(lsum12-lsum1-lsum2)
    if not 0 <= relation < 1:
        raise RuntimeError("Independent coloc H3 summation unstable")
    log_h3 = math.log(1e-4) * 2 + lsum1 + lsum2 + math.log1p(-relation)
    log_h = [0.0, math.log(1e-4)+lsum1, math.log(1e-4)+lsum2,
             log_h3, math.log(p12)+lsum12]
    denominator = logsumexp(log_h)
    return [math.exp(x-denominator) for x in log_h]


def main() -> None:
    candidates = read(CANDIDATES)
    expected_pairs = {r["candidate_locus_id"]: r["pair_id"] for r in candidates}
    if len(expected_pairs) != 25 or len(candidates) != 25:
        raise RuntimeError("Protected candidate family is not 25 unique rows")
    manifest = read(OUT / "fine_mapping_input_manifest.tsv")
    pair_manifest = read(OUT / "trait_coloc_input_manifest.tsv")
    fine = read(OUT / "fine_mapping_50.tsv")
    coloc = read(OUT / "trait_coloc_25.tsv")
    cs = read(OUT / "fine_mapping_credible_sets.tsv")
    pip = read(SSD / "variant_pip.tsv")
    sensitivity = read(OUT / "trait_coloc_prior_sensitivity.tsv")
    if not (len(manifest) == len(fine) == 50 and len(pair_manifest) == len(coloc) == 25):
        raise RuntimeError("Incomplete 50/25 family")
    expected_traits = {(cid, trait) for cid, pair in expected_pairs.items() for trait in pair.split("__")}
    manifest_map = {(r["candidate_locus_id"], r["trait_id"]): r for r in manifest}
    fine_map = {(r["candidate_locus_id"], r["trait_id"]): r for r in fine}
    if set(manifest_map) != expected_traits or set(fine_map) != expected_traits:
        raise RuntimeError("Fine-map candidate-trait identity mismatch")
    if {r["candidate_locus_id"] for r in pair_manifest} != set(expected_pairs) or \
       {r["candidate_locus_id"] for r in coloc} != set(expected_pairs):
        raise RuntimeError("Coloc candidate identity mismatch")

    groups = defaultdict(list)
    for row in pip:
        key = (row["candidate_locus_id"], row["trait_id"])
        groups[key].append(row)
    cs_groups = defaultdict(list)
    for row in cs:
        cs_groups[(row["candidate_locus_id"], row["trait_id"])].append(row)
    n_admitted = 0
    for key in expected_traits:
        m, f = manifest_map[key], fine_map[key]
        if m["pair_id"] != expected_pairs[key[0]] or f["pair_id"] != expected_pairs[key[0]]:
            raise RuntimeError(f"Pair identity mismatch: {key}")
        if m["status"] != "READY":
            if f["status"] != m["status"] or groups[key] or cs_groups[key]:
                raise RuntimeError(f"NOT_RUN mismatch: {key}")
            continue
        n_admitted += 1
        if sha256(Path(m["input_tsv"])) != m["input_sha256"]:
            raise RuntimeError(f"Input SHA mismatch: {key}")
        source = read(Path(m["input_tsv"]))
        if f["status"] not in {"EXPLORATORY_ABF_CREDIBLE_SET", "NO_ASSOCIATION_SUPPORTED_CREDIBLE_SET"}:
            raise RuntimeError(f"Admitted input missing ABF result: {key}")
        if len(source) != int(f["n_snps"]) or len(source) != len(groups[key]):
            raise RuntimeError(f"ABF SNP count mismatch: {key}")
        if {r["SNP"] for r in source} != {r["snp"] for r in groups[key]}:
            raise RuntimeError(f"ABF SNP identity mismatch: {key}")
        assoc = float(f["PP_association"])
        null = float(f["PP_null"])
        log_bfs = {r["SNP"]: log_bf(float(r["BETA"]), float(r["SE"])) for r in source}
        independent_null = math.exp(math.log(1-len(source)*1e-4) -
                                    logsumexp([math.log(1-len(source)*1e-4),
                                               math.log(1e-4)+logsumexp(list(log_bfs.values()))]))
        if not approx(independent_null,null):
            raise RuntimeError(f"Independent ABF null posterior mismatch: {key}")
        if not approx(assoc + null, 1) or not approx(sum(float(r["SNP.PP"]) for r in groups[key]), assoc):
            raise RuntimeError(f"ABF posterior mass mismatch: {key}")
        if not approx(sum(float(r["conditional_pip"]) for r in groups[key]), 1):
            raise RuntimeError(f"Conditional PIP mass mismatch: {key}")
        ordered = sorted(groups[key], key=lambda r: (-float(r["conditional_pip"]), r["snp"]))
        top_value = next((float(r["conditional_pip"]) for r in groups[key] if r["snp"]==f["top_snp"]), None)
        if top_value is None or not approx(top_value,float(ordered[0]["conditional_pip"])) or \
           not approx(top_value,float(f["top_snp_conditional_pip"])):
            raise RuntimeError(f"ABF top-SNP mismatch: {key}")
        if f["status"] == "NO_ASSOCIATION_SUPPORTED_CREDIBLE_SET":
            if assoc >= 0.8 or cs_groups[key] or int(f["credible_set_size"]) != 0:
                raise RuntimeError(f"No-credible-set rule violated: {key}")
        else:
            subset = cs_groups[key]
            if assoc < 0.8 or len(subset) != int(f["credible_set_size"]):
                raise RuntimeError(f"Credible-set rule violated: {key}")
            if [r["snp"] for r in subset] != [r["snp"] for r in ordered[:len(subset)]]:
                raise RuntimeError(f"Credible-set prefix mismatch: {key}")
            cover = sum(float(r["conditional_pip"]) for r in subset)
            previous = cover - float(subset[-1]["conditional_pip"])
            if cover < 0.95 - 1e-12 or previous >= 0.95 or not approx(cover,float(f["credible_set_cumulative_pip"])):
                raise RuntimeError(f"Credible-set coverage mismatch: {key}")

    pair_map = {r["candidate_locus_id"]: r for r in pair_manifest}
    executed = 0
    for row in coloc:
        cid = row["candidate_locus_id"]
        m = pair_map[cid]
        if row["pair_id"] != expected_pairs[cid]:
            raise RuntimeError(f"Coloc pair mismatch: {cid}")
        if m["status"] != "READY":
            if row["status"] != m["status"]:
                raise RuntimeError(f"Coloc NOT_RUN mismatch: {cid}")
            continue
        executed += 1
        if row["status"] != "EXPLORATORY_TRAIT_COLOC_ABF" or int(row["n_shared"]) != int(m["n_shared"]):
            raise RuntimeError(f"Missing coloc result: {cid}")
        if sha256(Path(m["trait1_input_tsv"])) != m["trait1_sha256"] or \
           sha256(Path(m["trait2_input_tsv"])) != m["trait2_sha256"]:
            raise RuntimeError(f"Coloc input SHA mismatch: {cid}")
        a,b = (read(Path(m[k])) for k in ("trait1_input_tsv","trait2_input_tsv"))
        if [x["SNP"] for x in a] != [x["SNP"] for x in b] or \
           [x["A1"] for x in a] != [x["A1"] for x in b]:
            raise RuntimeError(f"Coloc allele mismatch: {cid}")
        values = [float(row[f"PP.H{i}"]) for i in range(5)]
        expected = expected_coloc(a,b,1e-5)
        if any(not approx(x,y) for x,y in zip(values,expected)):
            raise RuntimeError(f"Independent coloc ABF recomputation mismatch: {cid}")
        if any(x<0 or x>1 for x in values) or not approx(sum(values),1):
            raise RuntimeError(f"Coloc posterior vector invalid: {cid}")
        conditional = values[4]/(values[3]+values[4]) if values[3]+values[4]>0 else math.nan
        if not approx(conditional,float(row["PP.H4_conditional_H3H4"])):
            raise RuntimeError(f"Coloc conditional H4 invalid: {cid}")
        support = values[4]>=0.8 and conditional>=0.8 and \
                  float(row["trait1_PP_association_shared"])>=0.8 and \
                  float(row["trait2_PP_association_shared"])>=0.8
        if row["descriptive_support"] != ("ABF_MODEL_SUPPORT" if support else "ABF_MODEL_SUPPORT_NOT_MET"):
            raise RuntimeError(f"Coloc support label mismatch: {cid}")
    if n_admitted != 28 or executed != 6 or len(sensitivity) != 12:
        raise RuntimeError("Unexpected executed family size")
    if Counter((r["candidate_locus_id"],r["p12"]) for r in sensitivity) != \
       Counter({(cid,p):1 for cid,m in pair_map.items() if m["status"]=="READY" for p in ("1e-06","5e-05")}):
        raise RuntimeError("Prior-sensitivity family mismatch")
    for row in sensitivity:
        if not approx(sum(float(row[f"PP.H{i}"]) for i in range(5)),1):
            raise RuntimeError("Prior sensitivity posterior mass invalid")
    ld_qc = read(V1 / "ld_numeric_qc.tsv")
    if len(ld_qc) != 32 or any(r["valid_for_v1_susie"] != "NO" for r in ld_qc):
        raise RuntimeError("v1 SuSiE LD gate unexpectedly changed")
    description = (LIB / "coloc/DESCRIPTION").read_text()
    if "Version: 5.2.3" not in description:
        raise RuntimeError("Pinned coloc package version changed")
    result_files = ("fine_mapping_50.tsv", "fine_mapping_credible_sets.tsv", "trait_coloc_25.tsv",
                    "trait_coloc_prior_sensitivity.tsv")
    report = {
        "status": "PASS_INTEGRITY_ONLY",
        "method": "EXPLORATORY_SINGLE_SIGNAL_ABF",
        "single_causal_variant_assumption_verified": False,
        "candidate_count": 25,
        "fine_mapping_unit_count": 50,
        "fine_mapping_admitted_count": n_admitted,
        "fine_mapping_credible_set_count": sum(x["status"]=="EXPLORATORY_ABF_CREDIBLE_SET" for x in fine),
        "trait_coloc_executed_count": executed,
        "trait_coloc_abf_model_support_count": sum(x["descriptive_support"]=="ABF_MODEL_SUPPORT" for x in coloc),
        "v1_susie_ld_matrices_passing": 0,
        "source_candidate_sha256": sha256(CANDIDATES),
        "v1_input_audit_provenance_sha256": sha256(V1 / "input_audit_provenance.json"),
        "method_lock_sha256": sha256(OUT / "PRE_OUTCOME_ABF_METHOD_LOCK.md"),
        "input_materialization_provenance_sha256": sha256(OUT / "materialization_provenance.json"),
        "materializer_script_sha256": sha256(REPO / "brain6/scripts/materialize_brain6_exploratory_abf_v2.py"),
        "runner_script_sha256": sha256(REPO / "brain6/scripts/run_brain6_exploratory_abf_v2.R"),
        "validation_script_sha256": sha256(Path(__file__)),
        "result_sha256": {name:sha256(OUT / name) for name in result_files},
        "results_report_sha256": sha256(OUT / "RESULTS_AND_LIMITATIONS.md"),
        "full_variant_pip_path": str(SSD / "variant_pip.tsv"),
        "full_variant_pip_sha256": sha256(SSD / "variant_pip.tsv"),
        "coloc_version": "5.2.3",
    }
    (OUT / "integrity_provenance.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:report[k] for k in ("candidate_count","fine_mapping_admitted_count",
                     "fine_mapping_credible_set_count","trait_coloc_executed_count",
                     "trait_coloc_abf_model_support_count")},sort_keys=True))


if __name__ == "__main__":
    main()
