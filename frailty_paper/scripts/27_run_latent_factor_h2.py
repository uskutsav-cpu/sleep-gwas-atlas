#!/usr/bin/env python3
"""Resumable source-QC -> harmonization -> HapMap3 munging -> LDSC h2.

Runs only the seven Catalog-mapped latent frailty factors in the frozen
secondary family. All large derived inputs are routed to external storage;
small logs and the result table stay in the repository.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
INDEX = REPO / "frailty_paper/manifests/latent_frailty_accession_index.tsv"
ACQUIRED = REPO / "frailty_paper/manifests/all_acquired_resources.tsv"
CONFIG = REPO / "frailty_paper/config/secondary_frailty_factors_v1.tsv"
EXPECTED = {
    "GCST90624046": "frailty_general",
    "GCST90624047": "frailty_factor_1",
    "GCST90624048": "frailty_factor_2",
    "GCST90624049": "frailty_factor_3",
    "GCST90624050": "frailty_factor_4",
    "GCST90624051": "frailty_factor_5",
    "GCST90624052": "frailty_factor_6",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gzip_readable(path: Path) -> bool:
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for _ in handle:
                pass
        return True
    except (OSError, EOFError, UnicodeError):
        return False


def run(command: list[str], *, cwd: Path = REPO) -> None:
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def run_family(external_root: Path, repo: Path = REPO) -> Path:
    external_root = external_root.resolve()
    workspace = external_root / "analysis-workspace/frailty_v1"
    harmonized = workspace / "harmonized_secondary"
    munged = workspace / "munged_secondary"
    for path in (harmonized, munged):
        if not path.is_dir():
            path.mkdir(parents=True, exist_ok=True)
    logs = repo / "frailty_paper/results/frailty_v1/logs_latent_factors"
    logs.mkdir(parents=True, exist_ok=True)
    progress = repo / "frailty_paper/results/frailty_v1/latent_factor_run_progress.tsv"

    acquired = {row["accession"]: row for row in read_tsv(ACQUIRED)
                if row.get("resource_type") == "GWAS summary statistics"}
    index = {row["accession"]: row for row in read_tsv(INDEX)}
    missing = set(EXPECTED).difference(index) | set(EXPECTED).difference(acquired)
    if missing:
        raise RuntimeError(f"accession registry or acquired manifest is missing {sorted(missing)}")
    if set(EXPECTED).difference(index):
        raise RuntimeError("expected accession metadata is incomplete")

    ref = external_root / "analysis-workspace/reference"
    python = repo / ".ldsc-env/bin/python"
    if not python.is_file():
        raise RuntimeError(f"LDSC Python environment is missing: {python}")
    rows: list[dict[str, str]] = []
    for accession, trait in EXPECTED.items():
        source = repo / acquired[accession]["file"]
        if not source.is_file():
            raise RuntimeError(f"source file is missing: {source}")
        observed_sha = sha256(source)
        expected_sha = acquired[accession]["sha256"].lower()
        if observed_sha != expected_sha:
            raise RuntimeError(f"source checksum mismatch for {accession}")
        if index[accession]["label_status"] != "VERIFIED" or index[accession]["genome_assembly"] != "GRCh37" or index[accession]["sample_ancestry"] != "European":
            raise RuntimeError(f"metadata/build/ancestry gate failed for {accession}")

        hfile = harmonized / f"{trait}.harmonized.tsv.gz"
        qcfile = harmonized / f"{trait}.qc.txt"
        if hfile.is_file() and qcfile.is_file() and gzip_readable(hfile):
            qc = dict(line.rstrip("\n").split("\t", 1) for line in qcfile.open(encoding="utf-8") if "\t" in line)
            if qc.get("infile_sha256") != observed_sha or qc.get("output_build") != "hg19":
                raise RuntimeError(f"existing harmonized output has inconsistent provenance: {trait}")
            print(f"RESUME: verified harmonized output for {trait}", flush=True)
        else:
            run([sys.executable, "scripts/01_harmonize.py", "--trait", trait,
                 "--config", str(CONFIG), "--infile", str(source), "--outdir", str(harmonized),
                 "--source-build", "GRCh37"], cwd=repo)
            if not hfile.is_file() or not qcfile.is_file() or not gzip_readable(hfile):
                raise RuntimeError(f"harmonizer did not produce a complete output for {trait}")

        sumstats = munged / f"{trait}.sumstats.gz"
        munge_log = munged / f"{trait}.log"
        if sumstats.is_file() and munge_log.is_file() and "Conversion finished" in munge_log.read_text(errors="replace") and gzip_readable(sumstats):
            print(f"RESUME: verified munged output for {trait}", flush=True)
        else:
            prefix = munged / trait
            run([str(python), str(repo / "ldsc/munge_sumstats.py"),
                 "--sumstats", str(hfile), "--out", str(prefix),
                 "--merge-alleles", str(ref / "w_hm3.snplist"), "--chunksize", "500000"], cwd=repo)
            if not sumstats.is_file() or not munge_log.is_file() or not gzip_readable(sumstats):
                raise RuntimeError(f"munge did not produce a complete output for {trait}")

        h2log = logs / f"h2_{trait}.log"
        if h2log.is_file() and "Total Observed scale h2:" in h2log.read_text(errors="replace"):
            print(f"RESUME: verified h2 output for {trait}", flush=True)
        else:
            run([str(python), str(repo / "ldsc/ldsc.py"), "--h2", str(sumstats),
                 "--ref-ld-chr", str(ref / "eur_w_ld_chr") + "/",
                 "--w-ld-chr", str(ref / "eur_w_ld_chr") + "/",
                 "--out", str(logs / f"h2_{trait}")], cwd=repo)
        qc_text = qcfile.read_text(encoding="utf-8")
        qc = dict(line.rstrip("\n").split("\t", 1) for line in qc_text.splitlines() if "\t" in line)
        rows.append({
            "accession": accession,
            "trait": trait,
            "source_file": str(source),
            "source_bytes": str(source.stat().st_size),
            "source_sha256": observed_sha,
            "harmonized_file": str(hfile),
            "harmonized_bytes": str(hfile.stat().st_size),
            "harmonized_sha256": sha256(hfile),
            "qc_file": str(qcfile),
            "qc_sha256": sha256(qcfile),
            "qc_rows_in": qc.get("rows_in", ""),
            "qc_rows_out": qc.get("rows_out", ""),
            "munged_file": str(sumstats),
            "munged_bytes": str(sumstats.stat().st_size),
            "munged_sha256": sha256(sumstats),
            "munge_log_sha256": sha256(munge_log),
            "h2_log": str(h2log),
            "h2_log_sha256": sha256(h2log),
            "runner": str(Path(__file__).relative_to(repo)),
            "status": "H2_LOG_COMPLETE",
        })
        progress.parent.mkdir(parents=True, exist_ok=True)
        with progress.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"COMPLETE: {trait}", flush=True)

    out = repo / "frailty_paper/results/frailty_v1/h2_latent_factors.tsv"
    run([sys.executable, "scripts/05_collate.py", "--mode", "h2", "--config", str(CONFIG),
         "--logdir", str(logs), "--out", str(out)], cwd=repo)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path,
                        default=Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1"))
    args = parser.parse_args()
    output = run_family(args.external_root)
    print(f"WROTE {output}", flush=True)


if __name__ == "__main__":
    main()
