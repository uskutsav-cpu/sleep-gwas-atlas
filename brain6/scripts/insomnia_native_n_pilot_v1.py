#!/usr/bin/env python3
"""Freeze and materialize exact-source Jansen insomnia native-N diagnostic inputs."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/confirmatory_source_rescue_20260927/other_traits"
CONFIG = ROOT / "brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_v1.json"
RULE = CONFIG.with_name("insomnia_native_n_pilot_decision_rule_v1.json")
PROTOCOL = BASE / "INSOMNIA_NATIVE_N_PILOT_PROTOCOL_V1.md"
README = BASE / "Insomnia_sumstats_Jansenetal.readme_220525.txt"
SELECTION = BASE / "insomnia_native_n_88_loci.tsv"
SELECTION_RECEIPT = BASE / "insomnia_native_n_88_loci.selection.json"
SOURCE = Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/Insomnia_sumstats_Jansenetal.txt.gz")
CANONICAL = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1")
PROVENANCE = CANONICAL / "provenance.json"
LOCI = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
CANONICAL_AGGREGATE = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
REFERENCE_PROVENANCE = Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json")
REFERENCE_PREFIX = ROOT / "ref/lava/ukb_v1.1/lava-ukb-v1.1"
R_RUNNER = ROOT / "brain6/scripts/run_insomnia_native_n_pilot_v1.R"
PYTHON_RUNNER = ROOT / "brain6/scripts/run_insomnia_native_n_pilot_v1.py"
RSCRIPT = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")
OUTPUT = Path("/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v1")
HEADER = ["SNP", "UNIQUE_ID", "CHR", "BP", "A1", "A2", "MAF", "OR", "SE", "P", "N", "INFO"]
SOURCE_SHA = "32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b"
README_SHA = "a5e92e264da88054173fc441e05f4587a6c6134078cc12b69e7294986fcc84a0"
SELECTION_SHA = "95b87274638396ddc9219177dd7209f92ba6c559e43845aed19c896982d6b1ae"
SELECTION_RECEIPT_SHA = "e29de75678dd8331ee0e7ff26fee4641a93aeb41d4fcb0bc33cca2d3b752093b"
PROVENANCE_SHA = "c50dae720237b548a282d931118f61348821a139f198c64d6d6876a015ce3d02"
LOCI_SHA = "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882"
CANONICAL_SHA = "ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352"
REFERENCE_SHA = "35ab371935b9c260dab21be248f4907b693565787ef9f62334e7d8030204270e"
RSCRIPT_SHA = "0948a03d938dd0d90b2436c87eca8e1dd59e6d0449bbf5d42bad8e9909102108"
CASES, CONTROLS = 109402, 277131
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check(path: Path, expected: str) -> None:
    if sha(path) != expected:
        raise ValueError(f"Frozen hash mismatch: {path}")


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def save_new(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write("\n")


def common_check(*, source: bool) -> None:
    fixed = {README: README_SHA, SELECTION: SELECTION_SHA,
             SELECTION_RECEIPT: SELECTION_RECEIPT_SHA, PROVENANCE: PROVENANCE_SHA,
             LOCI: LOCI_SHA, CANONICAL_AGGREGATE: CANONICAL_SHA,
             REFERENCE_PROVENANCE: REFERENCE_SHA, RSCRIPT: RSCRIPT_SHA}
    if source:
        fixed[SOURCE] = SOURCE_SHA
    for path, expected in fixed.items():
        check(path, expected)
    sel = json.loads(SELECTION_RECEIPT.read_text())
    if sel["n_loci"] != 88 or sel["table_sha256"] != SELECTION_SHA or sel["status_counts"] != {"NOT_RUN": 44, "TESTED": 44}:
        raise ValueError("Frozen locus selection receipt changed")


def selected_rows() -> list[dict[str, str]]:
    rows = table(SELECTION)
    if len(rows) != 88 or len({r["LOC"] for r in rows}) != 88:
        raise ValueError("Expected 88 unique selected loci")
    if Counter((r["CHR"], r["canonical_status"]) for r in rows) != {
        (str(chrom), status): 2 for chrom in range(1, 23) for status in ("TESTED", "NOT_RUN")
    }:
        raise ValueError("Frozen chromosome/status strata changed")
    return rows


def freeze() -> None:
    if CONFIG.exists() or OUTPUT.exists():
        raise FileExistsError("Frozen config or exclusive output root already exists")
    common_check(source=True)
    rows = selected_rows()
    if not REFERENCE_PREFIX.with_name(REFERENCE_PREFIX.name + "_chr1.info").exists():
        raise FileNotFoundError("Pinned LAVA reference prefix unavailable")
    data = {
        "schema_version": 1, "analysis_id": "brain6_insomnia_native_n_trait_only_pilot_v1",
        "scope": "ADMISSIBLE_SENSITIVITY_ONLY; trait-only technical diagnostic; no canonical or family promotion",
        "source_path": str(SOURCE), "source_sha256": SOURCE_SHA,
        "source_readme_sha256": README_SHA, "source_expected_rows": 10862567,
        "source_n_min": 366461, "source_n_max": 386533,
        "source_case_count": CASES, "source_control_count": CONTROLS,
        "source_n_semantics": "official CNCR README: literal per-SNP sample size",
        "binary_fraction_approximation": "fixed source-cohort cases/(cases+controls); variant-specific counts unavailable",
        "selection_table_sha256": SELECTION_SHA,
        "selection_receipt_sha256": SELECTION_RECEIPT_SHA,
        "locus_ids": [row["LOC"] for row in rows], "expected_loci": 88,
        "canonical_provenance_sha256": PROVENANCE_SHA,
        "canonical_aggregate_sha256": CANONICAL_SHA,
        "loci_file_sha256": LOCI_SHA,
        "reference_provenance_sha256": REFERENCE_SHA,
        "protocol_sha256": sha(PROTOCOL), "decision_rule_sha256": sha(RULE),
        "materializer_sha256": sha(Path(__file__)),
        "python_runner_sha256": sha(PYTHON_RUNNER), "r_runner_sha256": sha(R_RUNNER),
        "rscript_sha256": RSCRIPT_SHA,
        "worker_count": 4, "seed": 20260928, "strict_gate_p": 0.05 / 17465,
        "locus_processing": {"min_K": 2, "prune_thresh": 99, "max_prop_K": 0.75,
                             "drop_failed": True, "max_block_size": 3000, "cap_estimates": True},
        "univariate": {"cap_estimates": True},
        "output_root": str(OUTPUT),
    }
    save_new(CONFIG, data)
    print(json.dumps({"status": "FROZEN_BEFORE_OUTCOMES", "config": str(CONFIG),
                      "sha256": sha(CONFIG), "loci": 88, "workers": 4}, sort_keys=True))


def sign(a1: str, a2: str, raw1: str, raw2: str) -> int | None:
    a1, a2, raw1, raw2 = [x.upper() for x in (a1, a2, raw1, raw2)]
    if any(len(x) != 1 or x not in "ACGT" for x in (a1, a2, raw1, raw2)):
        return None
    if (a1, a2) == (raw1, raw2) or (a1, a2) == (raw1.translate(COMPLEMENT), raw2.translate(COMPLEMENT)):
        return 1
    if (a1, a2) == (raw2, raw1) or (a1, a2) == (raw2.translate(COMPLEMENT), raw1.translate(COMPLEMENT)):
        return -1
    return None


def write_shard(path: Path, rows: list[dict[str, str]], native_n: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=False)
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("SNP", "A1", "A2", "Z", "N"))
                for row in rows:
                    writer.writerow((row["SNP"], row["A1"], row["A2"], row["Z"], native_n[row["SNP"]]))


def materialize() -> None:
    config = json.loads(CONFIG.read_text())
    if config["materializer_sha256"] != sha(Path(__file__)) or config["python_runner_sha256"] != sha(PYTHON_RUNNER) or config["r_runner_sha256"] != sha(R_RUNNER):
        raise ValueError("Frozen pilot code changed")
    if config["protocol_sha256"] != sha(PROTOCOL) or config["decision_rule_sha256"] != sha(RULE):
        raise ValueError("Frozen protocol or advancement rule changed")
    if config["source_sha256"] != SOURCE_SHA or config["output_root"] != str(OUTPUT) or config["worker_count"] != 4:
        raise ValueError("Frozen source/output/workers changed")
    common_check(source=True)
    selected = selected_rows()
    if [r["LOC"] for r in selected] != config["locus_ids"]:
        raise ValueError("Selected loci differ from frozen config")
    if OUTPUT.exists():
        raise FileExistsError("Exclusive output root exists; never overwrite or restart")
    provenance = json.loads(PROVENANCE.read_text())
    original_records = {str(r["locus_id"]): r for r in provenance["records"]
                        if r.get("kind") == "sumstats" and r.get("trait") == "insomnia"}
    original: dict[str, list[dict[str, str]]] = {}
    targets: set[str] = set()
    for loc in selected:
        key = loc["LOC"]
        record = original_records.get(key)
        shard = CANONICAL / f"locus_{key}/insomnia.sumstats.tsv.gz"
        if not record or record["path"] != str(shard) or sha(shard) != record["sha256"]:
            raise ValueError(f"Original shard receipt mismatch: {key}")
        with gzip.open(shard, "rt", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames != ["SNP", "A1", "A2", "Z", "N"]:
                raise ValueError(f"Original shard schema mismatch: {key}")
            rows = list(reader)
        if len(rows) != record["rows"] or len({r["SNP"] for r in rows}) != len(rows):
            raise ValueError(f"Original shard row identity mismatch: {key}")
        if any(not math.isclose(float(r["N"]), 313750.03595553286, rel_tol=1e-12) for r in rows):
            raise ValueError(f"Original fixed effective N mismatch: {key}")
        original[key] = rows
        targets.update(r["SNP"] for r in rows)
    found: dict[str, list[str]] = {}
    duplicates: set[str] = set()
    source_rows = 0
    with gzip.open(SOURCE, "rt", encoding="utf-8", newline="") as stream:
        if stream.readline().rstrip("\r\n").split("\t") != HEADER:
            raise ValueError("Exact Jansen source header changed")
        for line in stream:
            source_rows += 1
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != len(HEADER):
                raise ValueError(f"Malformed source row {source_rows}")
            snp = parts[0]
            if snp in targets:
                if snp in found:
                    duplicates.add(snp)
                else:
                    found[snp] = parts
    if source_rows != config["source_expected_rows"] or duplicates or set(found) != targets:
        raise ValueError(f"Source row count, duplicate or missing target mismatch: rows={source_rows}, duplicated={len(duplicates)}, missing={len(targets-set(found))}")
    native_n: dict[str, int] = {}
    for loc in selected:
        key = loc["LOC"]
        for row in original[key]:
            src = found[row["SNP"]]
            try:
                chrom, bp, maf, odds, se, p, n, info = (int(src[2]), int(src[3]), float(src[6]),
                    float(src[7]), float(src[8]), float(src[9]), int(src[10]), float(src[11]))
                z = float(row["Z"])
            except (ValueError, OverflowError) as exc:
                raise ValueError(f"Invalid source/canonical numeric at {key}:{row['SNP']}") from exc
            orientation = sign(row["A1"], row["A2"], src[4], src[5])
            if chrom != int(loc["CHR"]) or not int(loc["START"]) <= bp <= int(loc["STOP"]):
                raise ValueError(f"Source position mismatch at {key}:{row['SNP']}")
            if orientation is None:
                raise ValueError(f"Source allele mismatch at {key}:{row['SNP']}")
            if not all(map(math.isfinite, (maf, odds, se, p, info, z))) or odds <= 0 or se <= 0 or not 0 < p <= 1:
                raise ValueError(f"Source statistic invalid at {key}:{row['SNP']}")
            if not (0.01 < maf <= 0.5 and 0.6 < info <= 1):
                raise ValueError(f"Source QC regression at {key}:{row['SNP']}")
            if not config["source_n_min"] <= n <= config["source_n_max"]:
                raise ValueError(f"Source N out of frozen range at {key}:{row['SNP']}")
            if not math.isclose(z, orientation * math.log(odds) / se, rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"Source Z mismatch at {key}:{row['SNP']}")
            old = native_n.setdefault(row["SNP"], n)
            if old != n:
                raise ValueError(f"Source SNP N inconsistency: {row['SNP']}")
    OUTPUT.mkdir(parents=True, exist_ok=False)
    (OUTPUT / "inputs").mkdir()
    (OUTPUT / "results").mkdir()
    loci_out = OUTPUT / "selected.loci"
    with loci_out.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("LOC", "CHR", "START", "STOP"))
        for row in selected:
            writer.writerow((row["LOC"], row["CHR"], row["START"], row["STOP"]))
    records = {}
    for loc in selected:
        key = loc["LOC"]
        directory = OUTPUT / "inputs" / f"locus_{key}"
        shard = directory / "insomnia.sumstats.tsv.gz"
        write_shard(shard, original[key], native_n)
        info = directory / "input_info.tsv"
        info.write_text("phenotype\tcases\tcontrols\tfilename\n" +
                        f"insomnia\t{CASES}\t{CONTROLS}\t{shard.name}\n")
        records[key] = {"chromosome": int(loc["CHR"]), "canonical_status": loc["canonical_status"],
                        "original_rows": len(original[key]), "rows": len(original[key]),
                        "original_shard_sha256": original_records[key]["sha256"],
                        "sumstats_sha256": sha(shard), "input_info_sha256": sha(info)}
    receipt = {"analysis_id": config["analysis_id"], "status": "PASS_EXACT_ROW_SET_NATIVE_N_MATERIALIZATION",
               "scope": config["scope"], "config_sha256": sha(CONFIG), "source_sha256": SOURCE_SHA,
               "source_rows_scanned": source_rows, "target_source_ids": len(targets),
               "source_n_min_observed_selected": min(native_n.values()),
               "source_n_max_observed_selected": max(native_n.values()),
               "selection_sha256": SELECTION_SHA, "selected_loci_sha256": sha(loci_out),
               "canonical_provenance_sha256": PROVENANCE_SHA,
               "reference_provenance_sha256": REFERENCE_SHA,
               "materializer_sha256": config["materializer_sha256"], "records": records,
               "canonical_modified": False, "lava_run_started": False}
    receipt_path = OUTPUT / "materialization.receipt.json"
    save_new(receipt_path, receipt)
    for worker in range(4):
        worker_ids = config["locus_ids"][worker::4]
        worker_cfg = {"analysis_id": config["analysis_id"], "scope": config["scope"],
                      "worker_id": worker+1, "locus_ids": worker_ids, "expected_loci": len(worker_ids),
                      "input_root": str(OUTPUT / "inputs"), "loci_file": str(loci_out),
                      "reference_prefix": str(REFERENCE_PREFIX), "trait_id": "insomnia",
                      "random_seed": config["seed"], "strict_gate_p": config["strict_gate_p"],
                      "output_tsv": str(OUTPUT / "results" / f"worker_{worker+1}.tsv"),
                      "summary_json": str(OUTPUT / "results" / f"worker_{worker+1}.summary.json"),
                      "pilot_config_sha256": sha(CONFIG),
                      "materialization_receipt_sha256": sha(receipt_path),
                      "execution_policy": {"locus_processing": config["locus_processing"],
                                           "univariate": config["univariate"],
                                           "worker_count": 4, "blas_threads": 1,
                                           "openmp_threads": 1,
                                           "numerical_stability_patch": "LAVA015_BLOCK_REDUCED_SYMMETRY_V1"}}
        save_new(OUTPUT / f"worker_{worker+1}.config.json", worker_cfg)
    print(json.dumps({"status": receipt["status"], "loci": 88, "rows": sum(r["rows"] for r in records.values()),
                      "workers": 4, "receipt_sha256": sha(receipt_path)}, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "materialize"))
    args = parser.parse_args()
    if args.action == "freeze":
        freeze()
    else:
        materialize()


if __name__ == "__main__":
    main()
