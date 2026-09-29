#!/usr/bin/env python3
"""Apply frozen source/coverage/signal gates to regional Brain6 GWAS–QTL pairs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/brain6_exploratory_functional_v1"
SSD = Path("/Volumes/Extreme SSD/brain6-work/brain6-exploratory-functional-v1")
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
DATASETS = {"QTD000171": "cortex_eqtl", "QTD000175": "cortex_sqtl",
            "QTD000176": "frontal_eqtl", "QTD000180": "frontal_sqtl"}
GATE_COLUMNS = ("candidate_locus_id", "region_grch37", "gwas_trait", "qtl_dataset_id", "qtl_context",
                "molecular_trait_id", "gene_id", "n_qtl_rows_region", "n_unique_shared_variants", "min_qtl_p_shared",
                "min_gwas_p_shared", "n_gwas_min", "n_gwas_max", "n_qtl_min", "n_qtl_max", "status", "input_path")
INPUT_COLUMNS = ("variant", "gwas_beta", "gwas_se", "gwas_p", "gwas_n_total", "gwas_n_case", "gwas_n_control",
                 "qtl_beta", "qtl_se", "qtl_p", "qtl_maf", "qtl_n")


def read_tsv(path: Path):
    with gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def key(chrom: str, pos: str, a: str, b: str) -> tuple[str, str, tuple[str, str]]:
    return (chrom, pos, tuple(sorted((a.upper(), b.upper()))))


def finite(row: dict, names: tuple[str, ...]) -> bool:
    try:
        return all(math.isfinite(float(row[name])) for name in names)
    except (ValueError, KeyError):
        return False


def main() -> None:
    if (OUT / "coloc_gate.tsv").exists():
        raise FileExistsError("Frozen gate output already exists")
    members = [r for r in read_tsv(MEMBERS) if r["pair_id"] == "insomnia__adhd"]
    if len(members) != 6:
        raise ValueError("Frozen eligible candidate count changed")
    (OUT / "coloc_inputs").mkdir(parents=True, exist_ok=True)
    result = []
    input_receipts = {}
    for candidate in members:
        region = candidate["region_group"]
        stem = region.replace(":", "_").replace("-", "_")
        for trait in ("insomnia", "adhd"):
            gwpath = SSD / "gwas_regional" / f"{trait}_{stem}.tsv.gz"
            gwr = json.loads((SSD / "gwas_regional" / f"{trait}_{stem}.receipt.json").read_text())
            if sha(gwpath) != gwr["sha256"]:
                raise ValueError(f"GWAS slice hash mismatch: {gwpath}")
            unique = defaultdict(list)
            for g in read_tsv(gwpath):
                if g["mapping_status"] == "UNIQUE_PLUS_SNP" and finite(g, ("beta_log_or", "se_log_or", "p", "n_total")):
                    if float(g["se_log_or"]) > 0 and float(g["n_total"]) > 0 and 0 < float(g["p"]) <= 1:
                        unique[key(g["chr38"], g["bp38"], g["a1"], g["a2"])].append(g)
            gwas = {k: v[0] for k, v in unique.items() if len(v) == 1}
            for dataset, context in DATASETS.items():
                qtpath = SSD / "qtl_regional" / f"{dataset}_{stem}.tsv.gz"
                qtr = json.loads((SSD / "qtl_regional" / f"{dataset}_{stem}.receipt.json").read_text())
                if sha(qtpath) != qtr["sha256"]:
                    raise ValueError(f"QTL slice hash mismatch: {qtpath}")
                groups = defaultdict(list)
                for q in read_tsv(qtpath):
                    if not finite(q, ("beta", "se", "pvalue", "maf", "an")):
                        continue
                    if not (float(q["se"]) > 0 and 0 < float(q["pvalue"]) <= 1 and
                            0 < float(q["maf"]) < 0.5 and float(q["an"]) > 0):
                        continue
                    groups[(q["molecular_trait_id"], q["gene_id"])].append(q)
                if not groups:
                    result.append({"candidate_locus_id": candidate["candidate_locus_id"], "region_grch37": region,
                                   "gwas_trait": trait, "qtl_dataset_id": dataset, "qtl_context": context,
                                   "molecular_trait_id": "NA", "gene_id": "NA", "n_qtl_rows_region": 0,
                                   "n_unique_shared_variants": 0, "min_qtl_p_shared": "NA", "min_gwas_p_shared": "NA",
                                   "n_gwas_min": "NA", "n_gwas_max": "NA", "n_qtl_min": "NA", "n_qtl_max": "NA",
                                   "status": "NO_QTL_ROWS_IN_REGION", "input_path": "NA"})
                for (molecular, gene), rows in sorted(groups.items()):
                    matched = defaultdict(list)
                    for q in rows:
                        k = key(q["chromosome"], q["position"], q["ref"], q["alt"])
                        if k in gwas:
                            matched[k].append(q)
                    shared = []
                    for k, qr in matched.items():
                        if len(qr) > 1:
                            measurements = {(x["beta"], x["se"], x["pvalue"], x["maf"], x["an"]) for x in qr}
                            if len(measurements) != 1:
                                continue
                        q, g = qr[0], gwas[k]
                        shared.append({"variant": q["variant"], "gwas_beta": g["beta_log_or"],
                                       "gwas_se": g["se_log_or"], "gwas_p": g["p"],
                                       "gwas_n_total": g["n_total"], "gwas_n_case": g["n_case"],
                                       "gwas_n_control": g["n_control"], "qtl_beta": q["beta"],
                                       "qtl_se": q["se"], "qtl_p": q["pvalue"], "qtl_maf": q["maf"],
                                       "qtl_n": float(q["an"]) / 2})
                    qmin = min((float(x["qtl_p"]) for x in shared), default=1.0)
                    gmin = min((float(x["gwas_p"]) for x in shared), default=1.0)
                    status = ("LOW_SHARED_VARIANT_COUNT" if len(shared) < 500 else
                              "NO_STRONG_QTL_SIGNAL" if qmin > 5e-8 else
                              "NO_STRONG_GWAS_SIGNAL" if gmin > 5e-8 else "ELIGIBLE_SINGLE_SIGNAL_ABF")
                    input_path = "NA"
                    if status == "ELIGIBLE_SINGLE_SIGNAL_ABF":
                        filename = f"{trait}_{dataset}_{stem}_{molecular.replace(':', '_')}.tsv.gz"
                        path = OUT / "coloc_inputs" / filename
                        with gzip.open(path, "xt", newline="") as stream:
                            writer = csv.DictWriter(stream, fieldnames=INPUT_COLUMNS, delimiter="\t", lineterminator="\n")
                            writer.writeheader()
                            writer.writerows(sorted(shared, key=lambda x: x["variant"]))
                        input_path = str(path.relative_to(ROOT))
                        input_receipts[input_path] = sha(path)
                    result.append({"candidate_locus_id": candidate["candidate_locus_id"], "region_grch37": region,
                                   "gwas_trait": trait, "qtl_dataset_id": dataset, "qtl_context": context,
                                   "molecular_trait_id": molecular, "gene_id": gene, "n_qtl_rows_region": len(rows),
                                   "n_unique_shared_variants": len(shared), "min_qtl_p_shared": qmin if shared else "NA",
                                   "min_gwas_p_shared": gmin if shared else "NA",
                                   "n_gwas_min": min((float(x["gwas_n_total"]) for x in shared), default="NA"),
                                   "n_gwas_max": max((float(x["gwas_n_total"]) for x in shared), default="NA"),
                                   "n_qtl_min": min((float(x["qtl_n"]) for x in shared), default="NA"),
                                   "n_qtl_max": max((float(x["qtl_n"]) for x in shared), default="NA"),
                                   "status": status, "input_path": input_path})
    with (OUT / "coloc_gate.tsv").open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=GATE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    with (OUT / "coloc_input_receipts.json").open("x") as f:
        json.dump({"source_gate_sha256": sha(OUT / "coloc_gate.tsv"), "input_sha256": input_receipts}, f,
                  indent=2, sort_keys=True)
        f.write("\n")
    print(json.dumps({"gate_rows": len(result), "eligible": sum(r["status"] == "ELIGIBLE_SINGLE_SIGNAL_ABF" for r in result),
                      "status_counts": {x: sum(r["status"] == x for r in result) for x in sorted({r["status"] for r in result})}},
                     indent=2))


if __name__ == "__main__":
    main()
