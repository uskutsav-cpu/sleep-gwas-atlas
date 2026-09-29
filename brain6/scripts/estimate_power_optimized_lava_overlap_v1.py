#!/usr/bin/env python3
"""Estimate source-specific LDSC intercepts for the duration LAVA sensitivity.

This creates inputs and outputs only beneath a new external sensitivity root.
It never reads or writes canonical LAVA result directories and refuses to
overwrite an existing output root.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CANDIDATE = ROOT / "brain6/results/power_optimized_sensitivity_v1/normalized/sleep_duration_continuous_dashti_2019/sumstats.tsv.gz"
DEFAULT_CANDIDATE_RECEIPT = ROOT / "brain6/results/power_optimized_sensitivity_v1/normalized/sleep_duration_continuous_dashti_2019/receipt.json"
DEFAULT_QC_MANIFEST = ROOT / "brain6/qc/dense_source_qc_master.tsv"
DEFAULT_REFERENCE_ROOT = Path("/Volumes/Extreme SSD/brain6-work/lava-power-optimized-sensitivity-v1/reference/eur_w_ld_chr")
DEFAULT_OUTPUT = Path("/Volumes/Extreme SSD/brain6-work/lava-power-optimized-sensitivity-v1/overlap_v1_retry")
DISORDERS = ("scz", "bipolar", "parkinson")
MUNGE_FIELDS = ("SNP", "A1", "A2", "Z", "P", "N")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def source_rows(manifest: Path) -> dict[str, dict[str, str]]:
    rows = {row["trait_id"]: row for row in table(manifest)}
    if any(trait not in rows for trait in DISORDERS):
        raise ValueError("Source QC manifest does not contain the three locked disorder inputs")
    return rows


def read_header(path: Path) -> tuple[csv.DictReader, object]:
    raw = gzip.open(path, "rt", encoding="utf-8", newline="")
    reader = csv.DictReader(raw, delimiter="\t")
    required = {"SNP", "A1", "A2", "BETA", "SE", "N"}
    if not required.issubset(set(reader.fieldnames or [])):
        raw.close()
        raise ValueError(f"Missing required GWAS fields in {path}: {sorted(required - set(reader.fieldnames or []))}")
    return reader, raw


def make_ldsc_input(source: Path, output: Path) -> dict[str, int | str]:
    if output.exists():
        raise FileExistsError(output)
    reader, raw = read_header(source)
    rows = 0
    with output.open("xb") as raw_out:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw_out, mtime=0) as zipped:
            with __import__("io").TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=MUNGE_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for row in reader:
                    try:
                        beta, se, n = float(row["BETA"]), float(row["SE"]), float(row["N"])
                        z = beta / se
                    except (TypeError, ValueError, ZeroDivisionError):
                        continue
                    snp, a1, a2 = row["SNP"].strip(), row["A1"].strip().upper(), row["A2"].strip().upper()
                    if (not snp or a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"}
                            or a1 == a2 or not all(map(math.isfinite, (z, n))) or n <= 0):
                        continue
                    p = math.erfc(abs(z) / math.sqrt(2.0))
                    if p == 0:
                        p = sys.float_info.min
                    writer.writerow({"SNP": snp, "A1": a1, "A2": a2,
                                     "Z": format(z, ".17g"), "P": format(p, ".17g"),
                                     "N": format(n, ".17g")})
                    rows += 1
    raw.close()
    if rows < 100_000:
        raise ValueError(f"Unexpectedly few usable summary-statistic rows from {source}: {rows}")
    return {"source": str(source), "source_sha256": sha256(source),
            "ldsc_input": str(output), "ldsc_input_sha256": sha256(output), "rows_written": rows}


def run(command: list[str], log: Path) -> None:
    with log.open("xb") as stream:
        result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                check=False, text=True)
    if result.returncode:
        raise RuntimeError(f"Command exited {result.returncode}; see {log}")


def intercept(path: Path, label: str) -> tuple[float, float | None]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if label == "Cross trait intercept":
        # LDSC 3 reports the cross-trait sampling intercept in its Genetic
        # Covariance section as the section's plain "Intercept". Older
        # releases label the same quantity "Cross trait intercept".
        match = re.search(r"(?m)^\s*Cross[- ]trait intercept:\s*([-+0-9.eE]+)(?:\s*\(([-+0-9.eE]+)\))?", text)
        if not match:
            section = re.search(r"(?s)Genetic Covariance\s*[-]+\s*(.*?)(?:\n\s*Genetic Correlation|\Z)", text)
            match = re.search(r"(?m)^\s*Intercept:\s*([-+0-9.eE]+)(?:\s*\(([-+0-9.eE]+)\))?",
                              section.group(1)) if section else None
        if not match:
            raise ValueError(f"Could not parse {label} from the LDSC Genetic Covariance section of {path}")
        return float(match.group(1)), float(match.group(2)) if match.group(2) else None
    match = re.search(rf"(?m)^\s*{re.escape(label)}:\s*([-+0-9.eE]+)(?:\s*\(([-+0-9.eE]+)\))?", text)
    if not match:
        raise ValueError(f"Could not parse {label} from {path}")
    return float(match.group(1)), float(match.group(2)) if match.group(2) else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--candidate-receipt", type=Path, default=DEFAULT_CANDIDATE_RECEIPT)
    parser.add_argument("--qc-manifest", type=Path, default=DEFAULT_QC_MANIFEST)
    parser.add_argument("--reference-root", type=Path, default=DEFAULT_REFERENCE_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--python", type=Path, default=ROOT / ".ldsc-env/bin/python")
    parser.add_argument("--ldsc-dir", type=Path, default=ROOT / "ldsc")
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite sensitivity LDSC output root: {out}")
    required = (args.candidate, args.candidate_receipt, args.qc_manifest,
                args.reference_root / "1.l2.ldscore.gz",
                args.reference_root / "1.l2.M_5_50",
                args.reference_root / "w_hm3.snplist",
                args.ldsc_dir / "munge_sumstats.py", args.ldsc_dir / "ldsc.py")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing inputs: " + ", ".join(missing))
    receipt = json.loads(args.candidate_receipt.read_text(encoding="utf-8"))
    if receipt.get("scientific_status") != "PASS" or receipt.get("status") != "COMPLETE":
        raise ValueError("Candidate normalization receipt does not pass")
    source_digest = sha256(args.candidate)
    if not any((args.candidate_receipt.parent / item.get("path", "")).resolve() == args.candidate.resolve()
               and item.get("sha256") == source_digest for item in receipt.get("outputs", [])):
        raise ValueError("Candidate summary statistics do not match their normalization receipt")
    rows = source_rows(args.qc_manifest)
    for trait in DISORDERS:
        source = Path(rows[trait]["source_path"])
        if not source.is_file() or sha256(source) != rows[trait]["source_sha256"]:
            raise ValueError(f"Source file absent or hash-mismatched for {trait}")
    out.mkdir(parents=True)
    (out / "inputs").mkdir()
    (out / "munge").mkdir()
    (out / "ldsc").mkdir()
    input_sources = {"sleep_duration_continuous_dashti_2019": args.candidate}
    input_sources.update({trait: Path(rows[trait]["source_path"]) for trait in DISORDERS})
    source_records = {}
    for trait, source in input_sources.items():
        source_records[trait] = make_ldsc_input(source, out / "inputs" / f"{trait}.tsv.gz")
    hm3 = args.reference_root / "w_hm3.snplist"
    python = str(args.python.resolve())
    munged: dict[str, Path] = {}
    for trait in input_sources:
        prefix = out / "munge" / trait
        command = [python, str(args.ldsc_dir / "munge_sumstats.py"),
                   "--sumstats", str(out / "inputs" / f"{trait}.tsv.gz"),
                   "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
                   "--signed-sumstats", "Z,0", "--N-col", "N",
                   "--merge-alleles", str(hm3), "--out", str(prefix)]
        run(command, out / "munge" / f"{trait}.log")
        munged[trait] = prefix.with_suffix(".sumstats.gz")
    ld_prefix = str(args.reference_root) + "/"
    for trait, path in munged.items():
        command = [python, str(args.ldsc_dir / "ldsc.py"), "--h2", str(path),
                   "--ref-ld-chr", ld_prefix, "--w-ld-chr", ld_prefix,
                   "--out", str(out / "ldsc" / f"h2_{trait}")]
        run(command, out / "ldsc" / f"h2_{trait}.stdout.log")
    pair_records = []
    sleep = "sleep_duration_continuous_dashti_2019"
    for trait in DISORDERS:
        command = [python, str(args.ldsc_dir / "ldsc.py"), "--rg",
                   f"{munged[sleep]},{munged[trait]}", "--ref-ld-chr", ld_prefix,
                   "--w-ld-chr", ld_prefix, "--out", str(out / "ldsc" / f"rg_{sleep}__{trait}")]
        run(command, out / "ldsc" / f"rg_{sleep}__{trait}.stdout.log")
        rg_log = out / "ldsc" / f"rg_{sleep}__{trait}.log"
        cross, cross_se = intercept(rg_log, "Cross trait intercept")
        h1, _ = intercept(out / "ldsc" / f"h2_{sleep}.log", "Intercept")
        h2, _ = intercept(out / "ldsc" / f"h2_{trait}.log", "Intercept")
        if h1 <= 0 or h2 <= 0 or cross_se is None or abs(cross) >= math.sqrt(h1 * h2):
            raise ValueError(f"Invalid LDSC overlap covariance for {trait}: h1={h1}, h2={h2}, cross={cross}")
        rho = cross / math.sqrt(h1 * h2)
        pair_id = f"sleep_duration_continuous_dashti_2019__{trait}"
        matrix = out / f"{pair_id}.txt"
        with matrix.open("x", encoding="utf-8") as stream:
            stream.write(f" {sleep} {trait}\n{sleep} {h1:.8f} {cross:.8f}\n{trait} {cross:.8f} {h2:.8f}\n")
        pair_records.append({"pair_id": pair_id, "trait1": sleep, "trait2": trait,
                             "sleep_intercept": h1, "disorder_intercept": h2,
                             "cross_trait_intercept": cross, "cross_trait_intercept_se": cross_se,
                             "error_correlation": rho, "matrix_path": str(matrix),
                             "matrix_sha256": sha256(matrix), "status": "PASS_LDSC_INTERCEPT_DERIVED"})
    result = {"schema_version": 1, "analysis_id": "brain6_power_optimized_sensitivity_v1_overlap",
              "status": "PASS_LDSC_INTERCEPTS_AND_MATRICES_CREATED",
              "source_candidates": source_records,
              "candidate_normalization_receipt_sha256": sha256(args.candidate_receipt),
              "reference_archive": {
                  "url": "https://zenodo.org/records/14993076/files/eur_w_ld_chr.tgz?download=1",
                  "doi": "10.5281/zenodo.14993076",
                  "archive_md5": "76c1890c8cf22d99d05c6707cc8441b4",
                  "archive_sha256": sha256(args.reference_root.parent / "eur_w_ld_chr.tgz"),
                  "downloaded_archive": str(args.reference_root.parent / "eur_w_ld_chr.tgz"),
                  "reference_files_root": str(args.reference_root),
                  "hapmap3_snplist_sha256": sha256(hm3)},
              "method": "LDSC cross-trait intercept from HM3-merged, allele-aware summary statistics; univariate LDSC intercepts supply the pair matrix diagonal; no zero-overlap assumption.",
              "pairs": pair_records}
    canonical_json(out / "overlap.provenance.json", result)
    print(json.dumps({"status": result["status"], "pairs": pair_records}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
