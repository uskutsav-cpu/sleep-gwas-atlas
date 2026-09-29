#!/usr/bin/env python3
"""Freeze and materialize an isolated MDD2025 trait-only LAVA pilot.

The freeze step selects loci without reading association statistics. The
materialization step refuses to run unless the frozen config and source match.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

from screen_mdd2025_lava_source_coverage_v1 import digest, load_reference

ROOT = Path(__file__).resolve().parents[2]
SOURCE_CONFIG = ROOT / "brain6/config/lava_multitrait_feasibility_v1/mdd2025_source_screen.json"
PILOT_CONFIG = ROOT / "brain6/config/lava_multitrait_feasibility_v1/mdd2025_trait_pilot_v1.json"
RUNNER = ROOT / "brain6/scripts/run_mdd2025_trait_pilot_v1.R"
OUTPUT = Path("/Volumes/Extreme SSD/brain6-work/lava-multitrait-feasibility-v1/mdd2025_trait_pilot_v1")
LOCI_FILE = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
REFERENCE = ROOT / "ref/lava/ukb_v1.1"
REFERENCE_PREFIX = REFERENCE / "lava-ukb-v1.1"
PROVENANCE = REFERENCE / "reference.provenance.json"
BASES = set("ACGT")
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
HEADER = ("#CHROM", "POS", "ID", "EA", "NEA", "BETA", "SE", "PVAL", "FCAS",
          "FCON", "IMPINFO", "NEFF", "NCAS", "NCON", "HETI", "HETDF", "HETPVAL")


def write_json_exclusive(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def loci_by_chromosome() -> dict[int, list[tuple[int, int, str]]]:
    result: dict[int, list[tuple[int, int, str]]] = {chrom: [] for chrom in range(1, 23)}
    with LOCI_FILE.open(encoding="utf-8") as stream:
        if stream.readline().split() != ["LOC", "CHR", "START", "STOP"]:
            raise ValueError("Frozen LAVA locus schema changed")
        for line in stream:
            if line.strip():
                loc, chrom, start, stop = line.split()
                result[int(chrom)].append((int(start), int(stop), loc))
    if sum(map(len, result.values())) != 2495:
        raise ValueError("Frozen LAVA locus count changed")
    return result


def freeze() -> None:
    if PILOT_CONFIG.exists():
        raise FileExistsError(PILOT_CONFIG)
    source_cfg = json.loads(SOURCE_CONFIG.read_text())
    if source_cfg["source_sha256"] != digest(Path(source_cfg["source_path"])):
        raise ValueError("MDD source identity changed")
    by_chrom = loci_by_chromosome()
    selected = []
    for chrom, loci in by_chrom.items():
        if len(loci) < 5:
            raise ValueError(f"Fewer than five frozen loci on chr{chrom}")
        ranked = sorted(loci, key=lambda loc: hashlib.sha256(f"brain6-mdd2025-pilot-v1:{loc[2]}".encode()).digest())
        selected.extend(loc[2] for loc in ranked[:5])
    selected.append("950")  # A previously documented seven-source MHC gap control.
    selected = sorted(set(selected), key=int)
    config = {
        "schema_version": 1,
        "analysis_id": "brain6_mdd2025_trait_only_pilot_v1",
        "scope": "Outcome-blinded, trait-only diagnostic; no canonical/family promotion",
        "selection": "Five lowest SHA-256-ranked frozen locus IDs per chromosome using literal brain6-mdd2025-pilot-v1 prefix, plus previously documented locus 950 gap control",
        "locus_ids": selected,
        "expected_loci": len(selected),
        "source_config_sha256": digest(SOURCE_CONFIG),
        "source_sha256": source_cfg["source_sha256"],
        "cohort_sidecar_md5": source_cfg["cohort_sidecar_md5"],
        "source_case_count": 412305,
        "source_control_count": 1588397,
        "source_N_semantics": "per-variant NEFF; source header explicitly defines NEFF as effective sample size",
        "loci_file_sha256": digest(LOCI_FILE),
        "reference_provenance_sha256": digest(PROVENANCE),
        "materializer_sha256": digest(Path(__file__)),
        "r_runner_sha256": digest(RUNNER),
        "minimum_info": 0.9,
        "minimum_maf": 0.01,
        "minimum_neff": 100,
        "exclude_palindromic": True,
        "require_exact_reference_snp_id_position_alleles": True,
        "quarantine_duplicate_ids_and_ambiguous_positions": True,
        "worker_count": 4,
        "seed": 20260922,
        "strict_gate_p": 0.05 / 17465,
        "locus_processing": {"min_K": 2, "prune_thresh": 99, "max_prop_K": 0.75,
                             "drop_failed": True, "max_block_size": 3000, "cap_estimates": True},
        "univariate": {"cap_estimates": True},
        "output_root": str(OUTPUT),
    }
    write_json_exclusive(PILOT_CONFIG, config)
    print(json.dumps({"frozen_config": str(PILOT_CONFIG), "sha256": digest(PILOT_CONFIG),
                      "loci": len(selected), "workers": 4}, sort_keys=True))


def materialize() -> None:
    config = json.loads(PILOT_CONFIG.read_text())
    source_cfg = json.loads(SOURCE_CONFIG.read_text())
    source = Path(source_cfg["source_path"])
    sidecar = source.parent / "pgc-mdd2025_no23andMe_eur_v3.49.24.11.txt"
    if config["materializer_sha256"] != digest(Path(__file__)) or config["r_runner_sha256"] != digest(RUNNER):
        raise ValueError("Frozen pilot code changed")
    if config["source_config_sha256"] != digest(SOURCE_CONFIG) or config["source_sha256"] != digest(source):
        raise ValueError("Frozen source identity changed")
    if config["cohort_sidecar_md5"] != digest(sidecar, "md5"):
        raise ValueError("Frozen cohort sidecar changed")
    if config["loci_file_sha256"] != digest(LOCI_FILE) or config["reference_provenance_sha256"] != digest(PROVENANCE):
        raise ValueError("Frozen LAVA loci or reference provenance changed")
    by_chrom = loci_by_chromosome()
    chosen = set(config["locus_ids"])
    if len(chosen) != config["expected_loci"] or len(chosen) < 110:
        raise ValueError("Frozen pilot selection invalid")
    selected = {chrom: sorted((loc for loc in loci if loc[2] in chosen)) for chrom, loci in by_chrom.items()}
    selected_starts = {chrom: [loc[0] for loc in loci] for chrom, loci in selected.items()}
    if sum(map(len, selected.values())) != len(chosen):
        raise ValueError("Frozen pilot selection contains unknown locus")
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite pilot root: {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    (OUTPUT / "inputs").mkdir()
    (OUTPUT / "results").mkdir()
    records = {}
    counters = Counter()
    current_chr = 0
    rows: dict[str, dict[str, tuple[int, str, str, float, float]]] = {}
    positions: dict[str, Counter[int]] = {}
    reference = {}

    def flush(chrom: int) -> None:
        for _, _, loc in selected[chrom]:
            loc_rows = rows.get(loc, {})
            pos_counts = positions.get(loc, Counter())
            accepted = [((pos, snp), snp, a1, a2, z, n) for snp, (pos, a1, a2, z, n) in loc_rows.items()
                        if snp and pos_counts[pos] == 1]
            accepted.sort()
            directory = OUTPUT / "inputs" / f"locus_{loc}"
            directory.mkdir()
            sumstats = directory / "mdd2025_no23andme.sumstats.tsv.gz"
            with sumstats.open("xb") as raw:
                with gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw, mtime=0) as packed:
                    with io.TextIOWrapper(packed, encoding="utf-8", newline="") as stream:
                        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                        writer.writerow(("SNP", "A1", "A2", "Z", "N"))
                        for _, snp, a1, a2, z, n in accepted:
                            writer.writerow((snp, a1, a2, format(z, ".17g"), format(n, ".17g")))
            info = directory / "input_info.tsv"
            with info.open("x", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("phenotype", "cases", "controls", "filename"))
                writer.writerow(("mdd2025_no23andme", config["source_case_count"],
                                 config["source_control_count"], sumstats.name))
            records[loc] = {"rows": len(accepted), "sumstats_sha256": digest(sumstats),
                            "input_info_sha256": digest(info), "chromosome": chrom}
            counters["written_rows"] += len(accepted)

    with gzip.open(source, "rt", encoding="utf-8", newline="") as stream:
        metadata = {}
        for line in stream:
            if line.startswith("##"):
                if line.startswith("##nCase=") or line.startswith("##nControl="):
                    name, value = line[2:].rstrip("\r\n").split("=", 1)
                    metadata[name] = int(value)
                continue
            if tuple(line.rstrip("\r\n").split("\t")) != HEADER:
                raise ValueError("PGC MDD2025 source header changed")
            break
        else:
            raise ValueError("Missing PGC source table")
        if metadata != {"nCase": config["source_case_count"], "nControl": config["source_control_count"]}:
            raise ValueError("Source case/control counts changed")
        col = {name: index for index, name in enumerate(HEADER)}
        for line in stream:
            counters["source_rows"] += 1
            parts = line.rstrip("\r\n").split("\t")
            try:
                chrom = int(parts[col["#CHROM"]])
                pos = int(parts[col["POS"]])
            except (ValueError, IndexError):
                continue
            if not 1 <= chrom <= 22:
                continue
            if chrom != current_chr:
                if chrom <= current_chr:
                    raise ValueError("Source chromosome order changed")
                if current_chr:
                    flush(current_chr)
                current_chr = chrom
                rows = {loc[2]: {} for loc in selected[chrom]}
                positions = {loc[2]: Counter() for loc in selected[chrom]}
                reference = load_reference(REFERENCE / f"lava-ukb-v1.1_chr{chrom}.info", chrom)
            loci = selected[chrom]
            index = bisect.bisect_right(selected_starts[chrom], pos) - 1
            if index < 0 or pos > loci[index][1]:
                continue
            loc = loci[index][2]
            counters["pilot_region_rows"] += 1
            snp = parts[col["ID"]].lower()
            hit = reference.get(snp)
            if not hit or hit[0] != pos:
                counters["reference_id_or_position_mismatch"] += 1
                continue
            a1, a2 = parts[col["EA"]].upper(), parts[col["NEA"]].upper()
            alleles = frozenset((a1, a2))
            if alleles != hit[1] or len(alleles) != 2 or not alleles <= BASES or alleles in PALINDROMIC:
                counters["reference_allele_or_palindromic"] += 1
                continue
            try:
                beta, se, p, fcase, fcontrol, info, neff, ncase, ncontrol = (
                    float(parts[col[k]]) for k in ("BETA", "SE", "PVAL", "FCAS", "FCON", "IMPINFO",
                                                   "NEFF", "NCAS", "NCON"))
                freq = (fcase * ncase + fcontrol * ncontrol) / (ncase + ncontrol)
                z = beta / se
            except (ValueError, ZeroDivisionError):
                counters["invalid_numeric"] += 1
                continue
            if not all(map(math.isfinite, (z, p, freq, info, neff))) or se <= 0 or not 0 <= p <= 1 or not 0 < freq < 1:
                counters["invalid_numeric"] += 1
                continue
            if info < config["minimum_info"] or min(freq, 1 - freq) < config["minimum_maf"] or neff < config["minimum_neff"]:
                counters["below_qc_threshold"] += 1
                continue
            positions[loc][pos] += 1
            if snp in rows[loc]:
                rows[loc][snp] = (pos, a1, a2, z, neff) if rows[loc][snp][0] == -1 else (-1, a1, a2, z, neff)
                counters["duplicate_id_rows"] += 1
            else:
                rows[loc][snp] = (pos, a1, a2, z, neff)
        if current_chr:
            flush(current_chr)
    if len(records) != len(chosen):
        raise ValueError("Pilot materialization omitted loci")
    receipt = {"analysis_id": config["analysis_id"], "status": "PASS_MATERIALIZED_DIAGNOSTIC_INPUTS",
               "scope": config["scope"], "config_sha256": digest(PILOT_CONFIG),
               "source_sha256": config["source_sha256"], "sidecar_sha256": digest(sidecar),
               "reference_provenance_sha256": config["reference_provenance_sha256"],
               "locus_file_sha256": config["loci_file_sha256"],
               "materializer_sha256": config["materializer_sha256"],
               "runner_sha256": config["r_runner_sha256"],
               "source_counters": dict(counters), "records": records}
    write_json_exclusive(OUTPUT / "materialization.receipt.json", receipt)
    materialization_sha = digest(OUTPUT / "materialization.receipt.json")
    ids = config["locus_ids"]
    for worker in range(config["worker_count"]):
        worker_ids = ids[worker::config["worker_count"]]
        worker_cfg = {
            "analysis_id": config["analysis_id"], "scope": config["scope"],
            "worker_id": worker + 1, "locus_ids": worker_ids,
            "expected_loci": len(worker_ids), "input_root": str(OUTPUT / "inputs"),
            "loci_file": str(LOCI_FILE), "reference_prefix": str(REFERENCE_PREFIX),
            "trait_id": "mdd2025_no23andme", "random_seed": config["seed"],
            "strict_gate_p": config["strict_gate_p"],
            "output_tsv": str(OUTPUT / "results" / f"worker_{worker + 1}.tsv"),
            "summary_json": str(OUTPUT / "results" / f"worker_{worker + 1}.summary.json"),
            "pilot_config_sha256": digest(PILOT_CONFIG),
            "materialization_receipt_sha256": materialization_sha,
            "execution_policy": {"locus_processing": config["locus_processing"],
                                 "univariate": config["univariate"],
                                 "worker_count": config["worker_count"],
                                 "blas_threads": 1, "openmp_threads": 1,
                                 "numerical_stability_patch": "LAVA015_BLOCK_REDUCED_SYMMETRY_V1"},
        }
        write_json_exclusive(OUTPUT / f"worker_{worker + 1}.config.json", worker_cfg)
    print(json.dumps({"status": receipt["status"], "pilot_loci": len(records),
                      "written_rows": counters["written_rows"],
                      "workers": config["worker_count"], "output": str(OUTPUT)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "materialize"))
    args = parser.parse_args()
    if args.action == "freeze":
        freeze()
    else:
        materialize()
