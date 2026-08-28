#!/usr/bin/env python3
"""Build the pre-analysis Pan-UKB EUR phenotype universe from official metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterable


KEY = ["trait_type", "phenocode", "pheno_sex", "coding", "modifier"]
QUANTITATIVE = {"continuous", "biomarkers"}
SOURCE_NAME = "Pan-UK Biobank 2025 EUR"
SOURCE_PMID = "40968291"
SOURCE_DOI = "10.1038/s41588-025-02335-7"
PHENOTYPE_MANIFEST_URL = (
    "https://pan-ukb-us-east-1.s3.amazonaws.com/"
    "sumstats_release/phenotype_manifest.tsv.bgz"
)
H2_MANIFEST_URL = (
    "https://pan-ukb-us-east-1.s3.amazonaws.com/"
    "sumstats_release/h2_manifest.tsv.bgz"
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(column, "") for column in KEY)


def number(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed


def integer_text(value: float) -> str:
    return str(int(value))


def stable_trait_id(row: dict[str, str]) -> str:
    pieces = [row[column] or "na" for column in KEY]
    normalized = [
        "".join(char.lower() if char.isalnum() else "_" for char in piece).strip("_")
        or "na"
        for piece in pieces
    ]
    return "panukbb_" + "__".join(normalized)


def phenotype_name(row: dict[str, str]) -> str:
    description = row.get("description", "")
    if description and description != "NA":
        return description
    phenocode = row.get("phenocode", "")
    if phenocode and phenocode != "NA":
        return phenocode.replace("_", " ")
    return row.get("category", "") or "UNRESOLVED"


def phenotype_definition(row: dict[str, str]) -> str:
    parts = [phenotype_name(row)]
    for field in ("coding_description", "description_more"):
        value = row.get(field, "")
        if value and value != "NA" and value not in parts:
            parts.append(value)
    return " | ".join(parts)


def eligible(
    phenotype: dict[str, str],
    heritability: dict[str, str],
    min_binary_cases: int,
    min_binary_controls: int,
    min_quantitative_n: int,
    min_h2_z: float,
) -> tuple[bool, str]:
    if phenotype.get("pheno_sex") != "both_sexes":
        return False, "not_both_sexes"
    if phenotype.get("phenotype_qc_EUR") != "PASS":
        return False, f"PanUKBB_EUR_QC_{phenotype.get('phenotype_qc_EUR', 'MISSING')}"
    h2_z = number(heritability.get("estimates.final.h2_z", ""))
    if h2_z is None or h2_z < min_h2_z:
        return False, "PanUKBB_final_h2_z_below_threshold"
    cases_or_n = number(phenotype.get("n_cases_EUR", ""))
    controls = number(phenotype.get("n_controls_EUR", ""))
    if cases_or_n is None:
        return False, "missing_EUR_sample_size"
    if phenotype["trait_type"] in QUANTITATIVE:
        if cases_or_n < min_quantitative_n:
            return False, "quantitative_n_below_threshold"
    elif cases_or_n < min_binary_cases or controls is None or controls < min_binary_controls:
        return False, "binary_case_control_count_below_threshold"
    return True, "PASS"


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phenotype-manifest",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/phenotype_manifest.tsv.bgz"),
    )
    parser.add_argument(
        "--h2-manifest",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/h2_manifest.tsv.bgz"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("discovery_extension/results/panukbb_eur_eligible_universe.tsv"),
    )
    parser.add_argument(
        "--run-metadata",
        type=Path,
        default=Path("discovery_extension/provenance/panukbb/universe_build.json"),
    )
    parser.add_argument("--min-binary-cases", type=int, default=1000)
    parser.add_argument("--min-binary-controls", type=int, default=1000)
    parser.add_argument("--min-quantitative-n", type=int, default=10000)
    parser.add_argument("--min-h2-z", type=float, default=4.0)
    args = parser.parse_args()

    for path in (args.phenotype_manifest, args.h2_manifest):
        if not path.is_file():
            fail(f"required official metadata file is missing: {path}")

    phenotypes = read_gzip_tsv(args.phenotype_manifest)
    h2_rows = [row for row in read_gzip_tsv(args.h2_manifest) if row.get("pop") == "EUR"]
    h2_by_key: dict[tuple[str, ...], dict[str, str]] = {}
    for row in h2_rows:
        row_key = key(row)
        if row_key in h2_by_key:
            fail(f"duplicate EUR h2 key: {row_key}")
        h2_by_key[row_key] = row

    output: list[dict[str, str]] = []
    exclusions: dict[str, int] = {}
    for phenotype in phenotypes:
        heritability = h2_by_key.get(key(phenotype))
        if heritability is None:
            exclusions["missing_EUR_h2_row"] = exclusions.get("missing_EUR_h2_row", 0) + 1
            continue
        keep, reason = eligible(
            phenotype,
            heritability,
            args.min_binary_cases,
            args.min_binary_controls,
            args.min_quantitative_n,
            args.min_h2_z,
        )
        if not keep:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue

        quantitative = phenotype["trait_type"] in QUANTITATIVE
        cases_or_n = number(phenotype["n_cases_EUR"])
        controls = number(phenotype.get("n_controls_EUR", ""))
        assert cases_or_n is not None
        if quantitative:
            sample_size = cases_or_n
            cases = "NA"
            controls_text = "NA"
            phenotype_type = "continuous"
        else:
            assert controls is not None
            sample_size = cases_or_n + controls
            cases = integer_text(cases_or_n)
            controls_text = integer_text(controls)
            phenotype_type = "binary"

        output.append(
            {
                "extension_trait_id": stable_trait_id(phenotype),
                "phenotype_name": phenotype_name(phenotype),
                "phenotype_category": phenotype["category"],
                "phenotype_definition": phenotype_definition(phenotype),
                "source": SOURCE_NAME,
                "study_accession": "PANUKBB_2025_RELEASE",
                "PMID": SOURCE_PMID,
                "DOI": SOURCE_DOI,
                "ancestry": "EUR",
                "sample_size": integer_text(sample_size),
                "cases": cases,
                "controls": controls_text,
                "build": "GRCh37",
                "source_url": (
                    "https://pan-ukb-us-east-1.s3.amazonaws.com/"
                    f"sumstats_flat_files/{phenotype['filename']}"
                ),
                "checksum": f"md5:{phenotype['md5_hex']}",
                "source_tabix_url": (
                    "https://pan-ukb-us-east-1.s3.amazonaws.com/"
                    f"sumstats_flat_files_tabix/{phenotype['filename_tabix']}"
                ),
                "source_tabix_checksum": f"md5:{phenotype['md5_hex_tabix']}",
                "source_tabix_size_bytes": phenotype["size_in_bytes_tabix"],
                "binary_or_continuous": phenotype_type,
                "available_beta": "YES:beta_EUR",
                "available_se": "YES:se_EUR",
                "available_effect_allele": "YES:alt",
                "available_other_allele": "YES:ref",
                "available_frequency": (
                    "YES:af_EUR" if quantitative else "YES:af_cases_EUR+af_controls_EUR"
                ),
                "available_info": "YES:separate_PanUKBB_variant_manifest",
                "h2": heritability.get("estimates.ldsc.h2_observed", "NA"),
                "h2_se": heritability.get("estimates.ldsc.h2_observed_se", "NA"),
                "h2_z": heritability.get("estimates.ldsc.h2_z", "NA"),
                "LDSC_intercept": heritability.get("estimates.ldsc.intercept", "NA"),
                "PanUKBB_QC_status": phenotype["phenotype_qc_EUR"],
                "PanUKBB_final_h2_z": heritability.get("estimates.final.h2_z", "NA"),
                "PanUKBB_LDSC_ratio": heritability.get("estimates.ldsc.ratio", "NA"),
                "in_max_independent_set": phenotype.get("in_max_independent_set", "false"),
                "trait_type": phenotype["trait_type"],
                "phenocode": phenotype["phenocode"],
                "pheno_sex": phenotype["pheno_sex"],
                "coding": phenotype.get("coding", ""),
                "modifier": phenotype.get("modifier", ""),
                "filename": phenotype["filename"],
                "file_size_bytes": phenotype["size_in_bytes"],
                "filename_tabix": phenotype["filename_tabix"],
            }
        )

    output.sort(key=lambda row: (row["phenotype_category"], row["phenotype_name"], row["extension_trait_id"]))
    if len({row["extension_trait_id"] for row in output}) != len(output):
        fail("generated extension_trait_id values are not unique")
    fields = list(output[0]) if output else []
    write_tsv(args.out, fields, output)

    metadata = {
        "schema_version": "1.0.0",
        "source_urls": [PHENOTYPE_MANIFEST_URL, H2_MANIFEST_URL],
        "phenotype_manifest_sha256": sha256(args.phenotype_manifest),
        "h2_manifest_sha256": sha256(args.h2_manifest),
        "source_phenotype_rows": len(phenotypes),
        "source_eur_h2_rows": len(h2_rows),
        "eligible_rows": len(output),
        "criteria": {
            "ancestry": "EUR",
            "pheno_sex": "both_sexes",
            "PanUKBB_phenotype_qc_EUR": "PASS",
            "PanUKBB_final_h2_z_minimum": args.min_h2_z,
            "minimum_binary_cases": args.min_binary_cases,
            "minimum_binary_controls": args.min_binary_controls,
            "minimum_quantitative_n": args.min_quantitative_n,
            "note": (
                "The Pan-UKB final h2 Z/QC fields are pre-analysis screening evidence. "
                "Every locked trait must still pass the extension's own standard LDSC h2 gate."
            ),
        },
        "exclusion_counts": dict(sorted(exclusions.items())),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
    }
    args.run_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.run_metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {args.out} eligible={len(output)} sha256={metadata['output_sha256']}")


if __name__ == "__main__":
    main()
