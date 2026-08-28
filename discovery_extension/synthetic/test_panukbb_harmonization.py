#!/usr/bin/env python3
"""Isolated smoke test for extension harmonization; writes no real result."""

from __future__ import annotations

import csv
import gzip
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def write_gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="sleep_gwas_extension_synthetic_") as temporary:
        root = Path(temporary)
        hm3 = root / "hm3.tsv.gz"
        variant = root / "variants.tsv.gz"
        source = root / "synthetic-continuous.tsv.bgz"
        binary_source = root / "synthetic-binary.tsv.bgz"
        panel = root / "panel.tsv"
        reference = root / "reference.tsv.gz"
        provenance = root / "reference.json"
        output = root / "harmonized.tsv.gz"
        binary_output = root / "harmonized-binary.tsv.gz"
        qc = root / "qc.tsv"
        binary_qc = root / "qc-binary.tsv"

        hm3_rows = [
            {"SNP": "rs1", "CHR": 1, "BP": 100, "A1": "G", "A2": "A"},
            {"SNP": "rs2", "CHR": 1, "BP": 200, "A1": "T", "A2": "C"},
            {"SNP": "rs3", "CHR": 6, "BP": 26000000, "A1": "C", "A2": "A"},
            {"SNP": "rs4", "CHR": 1, "BP": 300, "A1": "T", "A2": "A"},
            {"SNP": "rs5", "CHR": 1, "BP": 400, "A1": "G", "A2": "A"},
            {"SNP": "rs6", "CHR": 1, "BP": 500, "A1": "G", "A2": "A"},
        ]
        write_gzip_tsv(hm3, ["SNP", "CHR", "BP", "A1", "A2"], hm3_rows)
        variant_rows = [
            {"chrom": row["CHR"], "pos": row["BP"], "ref": row["A2"], "alt": row["A1"], "rsid": row["SNP"], "info": 0.8 if row["SNP"] == "rs5" else 0.99}
            for row in hm3_rows
        ]
        write_gzip_tsv(variant, ["chrom", "pos", "ref", "alt", "rsid", "info"], variant_rows)

        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/09_build_panukbb_hm3_reference.py"),
                "--variant-manifest", str(variant), "--hm3-map", str(hm3),
                "--out", str(reference), "--provenance-out", str(provenance),
            ],
            check=True,
        )

        source_rows = []
        for row in variant_rows:
            rsid = row["rsid"]
            source_rows.append({
                "chr": row["chrom"], "pos": row["pos"], "ref": row["ref"], "alt": row["alt"],
                "af_EUR": 0.005 if rsid == "rs6" else 0.2,
                "beta_EUR": 0.1, "se_EUR": 0.02, "neglog10_pval_EUR": 6,
                "low_confidence_EUR": "true" if rsid == "rs2" else "false",
            })
        source_rows.append({
            "chr": 1, "pos": 999, "ref": "A", "alt": "C", "af_EUR": 0.2,
            "beta_EUR": 0.1, "se_EUR": 0.02, "neglog10_pval_EUR": 6,
            "low_confidence_EUR": "false",
        })
        write_gzip_tsv(
            source,
            ["chr", "pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"],
            source_rows,
        )
        with panel.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, delimiter="\t",
                fieldnames=["extension_trait_id", "source_filename", "binary_or_continuous", "sample_size", "cases", "controls"],
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow({
                "extension_trait_id": "synthetic_trait", "source_filename": source.name,
                "binary_or_continuous": "continuous", "sample_size": 100000,
                "cases": "NA", "controls": "NA",
            })
            writer.writerow({
                "extension_trait_id": "synthetic_binary", "source_filename": binary_source.name,
                "binary_or_continuous": "binary", "sample_size": 4000,
                "cases": 1000, "controls": 3000,
            })

        write_gzip_tsv(
            binary_source,
            ["chr", "pos", "ref", "alt", "af_cases_EUR", "af_controls_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"],
            [{
                "chr": 1, "pos": 100, "ref": "A", "alt": "G",
                "af_cases_EUR": 0.3, "af_controls_EUR": 0.1,
                "beta_EUR": 0.2, "se_EUR": 0.04, "neglog10_pval_EUR": 5,
                "low_confidence_EUR": "false",
            }],
        )

        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/10_harmonize_panukbb.py"),
                "--trait-id", "synthetic_trait", "--source", str(source),
                "--panel", str(panel), "--reference", str(reference),
                "--out", str(output), "--qc-out", str(qc),
            ],
            check=True,
        )
        with gzip.open(output, "rt", newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(rows) != 1 or rows[0]["SNP"] != "rs1" or rows[0]["A1"] != "G":
            raise SystemExit(f"ERROR: unexpected synthetic survivors: {rows}")
        qc_rows = {row["metric"]: row["count"] for row in csv.DictReader(qc.open(), delimiter="\t")}
        expected = {
            "low_confidence_EUR": "1",
            "extended_MHC": "1",
            "strand_ambiguous": "1",
            "info_below_0_9_or_above_1": "1",
            "maf_below_0_01_or_invalid_frequency": "1",
            "not_exact_pinned_hapmap3_identity": "1",
            "output_rows": "1",
        }
        for metric, count in expected.items():
            if qc_rows.get(metric) != count:
                raise SystemExit(f"ERROR: synthetic QC {metric}={qc_rows.get(metric)} expected {count}")

        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/10_harmonize_panukbb.py"),
                "--trait-id", "synthetic_binary", "--source", str(binary_source),
                "--panel", str(panel), "--reference", str(reference),
                "--out", str(binary_output), "--qc-out", str(binary_qc),
            ],
            check=True,
        )
        with gzip.open(binary_output, "rt", newline="", encoding="utf-8") as handle:
            binary_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(binary_rows) != 1 or float(binary_rows[0]["FRQ"]) != 0.15 or float(binary_rows[0]["N"]) != 3000:
            raise SystemExit(f"ERROR: unexpected binary synthetic output: {binary_rows}")
    print("PANUKBB_HARMONIZATION_SYNTHETIC_OK continuous_survivors=1 binary_survivors=1 isolated=true")


if __name__ == "__main__":
    main()
