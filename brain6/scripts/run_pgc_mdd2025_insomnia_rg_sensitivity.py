#!/usr/bin/env python3
"""Run a reproducible LDSC insomnia–PGC MDD2025 external sensitivity.

PGC source data and all derived summary-statistics files stay on the external
data drive. The output directory is created exclusively and is never reused.
This analysis estimates genome-wide LDSC rg; it does not replace local LAVA.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PGC = Path("/Volumes/Extreme SSD/brain6-work/replication-pgc-mdd2025/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv.gz")
DEFAULT_README = Path("/Volumes/Extreme SSD/brain6-work/replication-pgc-mdd2025/README.pgc-mdd2025.txt")
DEFAULT_OUTPUT = Path("/Volumes/Extreme SSD/brain6-work/replication-pgc-mdd2025/insomnia_rg_sensitivity_v1")
DEFAULT_REFERENCE = Path("/Volumes/Extreme SSD/brain6-work/lava-power-optimized-sensitivity-v1/reference/eur_w_ld_chr")
DEFAULT_REFERENCE_ARCHIVE = DEFAULT_REFERENCE.parent / "eur_w_ld_chr.tgz"
PGC_MD5 = "29d3ce57cfce19ca28eb1643d9e2e428"
README_MD5 = "d743266ab7da04276517d3d4e033052a"
REFERENCE_MD5 = "76c1890c8cf22d99d05c6707cc8441b4"
INSOMNIA_TRAIT_ID = "insomnia"
FIELDS = ("SNP", "A1", "A2", "Z", "P", "N")


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_ldsc_input(source: Path, output: Path, source_kind: str) -> dict[str, object]:
    if output.exists():
        raise FileExistsError(output)
    rows = 0
    chromosomes: dict[str, int] = {}
    with gzip.open(source, "rt", encoding="utf-8", newline="") as raw:
        if source_kind == "pgc_mdd2025":
            for line in raw:
                if not line.startswith("##"):
                    header = line.rstrip("\r\n").lstrip("#").split("\t")
                    break
            else:
                raise ValueError("PGC source has no tabular header")
            reader = csv.DictReader(raw, fieldnames=header, delimiter="\t")
            cols = {"SNP": "ID", "A1": "EA", "A2": "NEA", "BETA": "BETA", "SE": "SE", "N": "NEFF"}
            chr_field = "CHROM"
        elif source_kind == "insomnia":
            reader = csv.DictReader(raw, delimiter="\t")
            cols = {name: name for name in ("SNP", "A1", "A2", "BETA", "SE", "N")}
            chr_field = "CHR"
        else:
            raise ValueError(f"Unsupported source type: {source_kind}")
        required = set(cols.values()) | {chr_field}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing source fields: {sorted(required - set(reader.fieldnames or []))}")
        with output.open("xb") as raw_out:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw_out, mtime=0) as compressed:
                with __import__("io").TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    for row in reader:
                        chrom = str(row.get(chr_field, "")).removeprefix("chr")
                        if chrom.isdigit() and 1 <= int(chrom) <= 22:
                            chromosomes[chrom] = chromosomes.get(chrom, 0) + 1
                        snp = str(row.get(cols["SNP"], "")).strip()
                        a1, a2 = (str(row.get(cols[k], "")).strip().upper() for k in ("A1", "A2"))
                        try:
                            beta, se, n = (float(row[cols[k]]) for k in ("BETA", "SE", "N"))
                            z = beta / se
                        except (KeyError, TypeError, ValueError, ZeroDivisionError):
                            continue
                        if (not snp or a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"}
                                or a1 == a2 or not all(map(math.isfinite, (z, n))) or se <= 0 or n <= 0):
                            continue
                        p = math.erfc(abs(z) / math.sqrt(2.0))
                        if p == 0:
                            p = sys.float_info.min
                        writer.writerow({"SNP": snp, "A1": a1, "A2": a2,
                                         "Z": format(z, ".17g"), "P": format(p, ".17g"), "N": format(n, ".17g")})
                        rows += 1
    return {"source": str(source), "source_sha256": digest(source), "ldsc_input": str(output),
            "ldsc_input_sha256": digest(output), "input_rows_written": rows,
            "autosomal_rows_seen": sum(chromosomes.values()), "autosomal_chromosome_counts": chromosomes}


def run_logged(command: list[str], log: Path, cwd: Path) -> None:
    with log.open("xb") as stream:
        result = subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=False, text=True)
    if result.returncode:
        raise RuntimeError(f"LDSC command failed with exit code {result.returncode}; see {log}")


def parse_rg_log(path: Path) -> dict[str, float | str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "rg": r"(?m)^\s*Genetic Correlation:\s*([-+0-9.eE]+)\s*\(([-+0-9.eE]+)\)",
        "p": r"(?m)^\s*P:\s*([-+0-9.eE]+)",
        "cross_trait_intercept": r"(?m)^\s*Cross[- ]trait intercept:\s*([-+0-9.eE]+)\s*\(([-+0-9.eE]+)\)",
    }
    result: dict[str, float | str] = {"log_path": str(path)}
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            if name == "cross_trait_intercept":
                section = re.search(r"(?s)Genetic Covariance\s*-+\s*(.*?)(?:\n\s*Genetic Correlation|\Z)", text)
                match = re.search(r"(?m)^\s*Intercept:\s*([-+0-9.eE]+)\s*\(([-+0-9.eE]+)\)", section.group(1)) if section else None
            if not match:
                raise ValueError(f"Could not parse {name} from {path}")
        result[name] = float(match.group(1))
        if match.lastindex and match.lastindex > 1:
            result[f"{name}_se"] = float(match.group(2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pgc", type=Path, default=DEFAULT_PGC)
    parser.add_argument("--readme", type=Path, default=DEFAULT_README)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--reference-archive", type=Path, default=DEFAULT_REFERENCE_ARCHIVE)
    parser.add_argument("--insomnia", type=Path, default=None)
    parser.add_argument("--python", type=Path, default=ROOT / ".ldsc-env/bin/python")
    parser.add_argument("--ldsc-dir", type=Path, default=ROOT / "ldsc")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite an existing output directory: {output}")
    qc = ROOT / "brain6/qc/dense_source_qc_master.tsv"
    with qc.open(encoding="utf-8", newline="") as stream:
        insomnia_row = next((r for r in csv.DictReader(stream, delimiter="\t") if r["trait_id"] == INSOMNIA_TRAIT_ID), None)
    if insomnia_row is None:
        raise ValueError("Frozen Brain6 source QC manifest has no insomnia row")
    insomnia = args.insomnia or Path(insomnia_row["source_path"])
    if digest(insomnia) != insomnia_row["source_sha256"]:
        raise ValueError("Insomnia source hash does not match the frozen Brain6 source QC manifest")
    checks = ((args.pgc, "md5", PGC_MD5), (args.readme, "md5", README_MD5),
              (args.reference_archive, "md5", REFERENCE_MD5))
    for path, algorithm, expected in checks:
        if not path.is_file() or digest(path, algorithm) != expected:
            raise ValueError(f"Missing or checksum-mismatched input: {path}")
    needed = [args.reference / f"{chrom}.l2.ldscore.gz" for chrom in range(1, 23)]
    needed += [args.reference / f"{chrom}.l2.M_5_50" for chrom in range(1, 23)]
    needed += [args.reference / "w_hm3.snplist", args.ldsc_dir / "munge_sumstats.py", args.ldsc_dir / "ldsc.py", args.python]
    missing = [str(path) for path in needed if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing LDSC/reference inputs: " + ", ".join(missing))
    output.mkdir(parents=True)
    inputs = output / "inputs"
    munged_dir = output / "munge"
    ldsc_dir = output / "ldsc"
    for directory in (inputs, munged_dir, ldsc_dir):
        directory.mkdir()
    provenance: dict[str, object] = {
        "schema_version": 1,
        "analysis_id": "brain6_insomnia_pgc_mdd2025_rg_sensitivity_v1",
        "status": "RUNNING",
        "claim_limit": "Genome-wide LDSC external sensitivity only; not local LAVA, an independent replication claim, or a replacement for the frozen discovery estimate.",
        "data_use": "PGC data retained only in the external research workspace; raw or transformed PGC summary statistics must not be committed or redistributed.",
        "sources": {
            "insomnia": {"trait_id": "jansen_2019_insomnia_ukb", "source_path": str(insomnia),
                         "source_sha256": digest(insomnia), "source_qc_manifest": str(qc),
                         "source_qc_manifest_sha256": digest(qc), "study_id": insomnia_row["study_id"]},
            "mdd": {"study_id": "pgc_mdd2025_no23andme_noUKBB_eur_v3_49_24_11", "source_path": str(args.pgc),
                    "source_md5": digest(args.pgc, "md5"), "source_sha256": digest(args.pgc),
                    "readme_path": str(args.readme), "readme_md5": digest(args.readme, "md5"),
                    "genome_build": "GRCh37", "ancestry": "EUR", "cases": 357636, "controls": 1281936,
                    "effective_n": 967084, "cohorts": 74, "reported_variants": 7139872,
                    "cohort_overlap_caveat": "The source excludes UK Biobank. This reduces direct overlap with the UKB-only Jansen insomnia GWAS, but exact participant linkage and overlap through other recruitment sources are not established."}},
        "reference": {"url": "https://zenodo.org/records/14993076/files/eur_w_ld_chr.tgz?download=1",
                      "doi": "10.5281/zenodo.14993076", "archive_path": str(args.reference_archive),
                      "archive_md5": digest(args.reference_archive, "md5"), "reference_dir": str(args.reference),
                      "hm3_snplist_sha256": digest(args.reference / "w_hm3.snplist")},
        "method": "LDSC 1.0.x two-trait --rg using HapMap3 allele-merged summary statistics and 1000 Genomes Phase 3 EUR LD scores/weights; PGC per-variant effective N is used.",
        "software": {"python": str(args.python), "python_version": subprocess.check_output([str(args.python), "--version"], text=True).strip(),
                     "ldsc_source": str(args.ldsc_dir), "ldsc_git_commit": subprocess.check_output(["git", "-C", str(args.ldsc_dir), "rev-parse", "HEAD"], text=True).strip()},
        "outputs": {"external_output_dir": str(output)},
    }
    save_json(output / "provenance.json", provenance)
    try:
        converted: dict[str, dict[str, object]] = {}
        converted["insomnia"] = write_ldsc_input(insomnia, inputs / "insomnia.tsv.gz", "insomnia")
        converted["mdd2025"] = write_ldsc_input(args.pgc, inputs / "mdd2025.tsv.gz", "pgc_mdd2025")
        for label in ("insomnia", "mdd2025"):
            source = inputs / f"{label}.tsv.gz"
            prefix = munged_dir / label
            command = [str(args.python), str(args.ldsc_dir / "munge_sumstats.py"), "--sumstats", str(source),
                       "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P", "--signed-sumstats", "Z,0",
                       "--N-col", "N", "--merge-alleles", str(args.reference / "w_hm3.snplist"), "--out", str(prefix)]
            run_logged(command, munged_dir / f"{label}.stdout.log", ROOT)
        ld_prefix = str(args.reference) + "/"
        for label in ("insomnia", "mdd2025"):
            command = [str(args.python), str(args.ldsc_dir / "ldsc.py"), "--h2",
                       str(munged_dir / f"{label}.sumstats.gz"), "--ref-ld-chr", ld_prefix,
                       "--w-ld-chr", ld_prefix, "--out", str(ldsc_dir / f"h2_{label}")]
            run_logged(command, ldsc_dir / f"h2_{label}.stdout.log", ROOT)
        rg_prefix = ldsc_dir / "rg_insomnia__mdd2025"
        command = [str(args.python), str(args.ldsc_dir / "ldsc.py"), "--rg",
                   f"{munged_dir / 'insomnia.sumstats.gz'},{munged_dir / 'mdd2025.sumstats.gz'}",
                   "--ref-ld-chr", ld_prefix, "--w-ld-chr", ld_prefix, "--out", str(rg_prefix)]
        run_logged(command, ldsc_dir / "rg_insomnia__mdd2025.stdout.log", ROOT)
        result = parse_rg_log(rg_prefix.with_suffix(".log"))
        result["status"] = "PASS_LDSC_EXTERNAL_SENSITIVITY"
        result["source_rows"] = converted
        result["munged_sumstats_sha256"] = {
            label: digest(munged_dir / f"{label}.sumstats.gz") for label in ("insomnia", "mdd2025")}
        result["rg_log_sha256"] = digest(rg_prefix.with_suffix(".log"))
        provenance["result"] = result
        provenance["status"] = "PASS_LDSC_EXTERNAL_SENSITIVITY"
        save_json(output / "provenance.json", provenance)
        print(json.dumps(result, indent=2, sort_keys=True))
    except BaseException as exc:
        provenance["status"] = "FAILED_PRESERVED_FOR_AUDIT"
        provenance["error"] = f"{type(exc).__name__}: {exc}"
        save_json(output / "provenance.json", provenance)
        raise


if __name__ == "__main__":
    main()
