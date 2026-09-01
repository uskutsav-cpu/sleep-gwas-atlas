#!/usr/bin/env python3
"""Build and verify the Track B dense full-summary-statistics QC gate."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
PANEL = Path("config/analysis_panel.tsv")
SCHEMAS = Path("config/gwas_schemas.tsv")
READINESS = Path("results/tables/pleiotropy_input_readiness.tsv")
OUT = Path("results/track_b/03_dense_input_qc.tsv")
LOCK = Path("results/track_b/03_dense_input_qc.lock.json")

TRAITS = {
    "snoring": ("A", "TRAIT1"),
    "parental_lifespan": ("A", "TRAIT2"),
    "insomnia": ("B;CONTROL", "TRAIT1"),
    "adhd": ("B", "TRAIT2"),
    "frailty": ("CONTROL", "TRAIT2"),
}

FIELDS = [
    "pair_ids", "trait_id", "pair_role", "source_gwas_id", "source",
    "publication", "ancestry", "source_build", "analysis_build",
    "dense_file", "dense_file_sha256", "compressed_size_bytes", "rows_in",
    "rows_out", "pct_retained", "observed_schema", "schema_status",
    "rsid_status", "effect_allele_status", "other_allele_status",
    "effect_status", "SE_status", "P_status", "EAF_status", "N_status",
    "case_control_counts", "INFO_status", "liftover_status",
    "duplicate_variant_handling", "multiallelic_handling",
    "strand_ambiguity_handling", "MAF_handling", "MHC_handling",
    "sample_size_convention", "fine_mapping_readiness", "overall_QC",
    "major_limitation",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_qc(path: Path) -> tuple[dict[str, str], dict[str, tuple[str, str]]]:
    metadata: dict[str, str] = {}
    steps: dict[str, tuple[str, str]] = {}
    in_steps = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if fields == ["step", "dropped", "remaining"]:
            in_steps = True
            continue
        if in_steps and len(fields) == 3:
            steps[fields[0]] = (fields[1], fields[2])
        elif len(fields) == 2:
            metadata[fields[0]] = fields[1]
    return metadata, steps


def find_step(steps: dict[str, tuple[str, str]], prefix: str) -> tuple[str, str]:
    matches = [value for key, value in steps.items() if key.startswith(prefix)]
    if len(matches) != 1:
        raise SystemExit(f"ERROR: expected one QC step starting {prefix!r}, found {len(matches)}")
    return matches[0]


def one(index: dict[str, dict[str, str]], schema_index: dict[str, dict[str, str]],
        readiness_index: dict[str, dict[str, str]], trait: str) -> dict[str, str]:
    panel = index[trait]
    schema = schema_index[panel["source_id"]]
    readiness = readiness_index[trait]
    path = Path(readiness["harmonized_file"])
    qc_path = path.with_name(f"{trait}.qc.txt")
    if readiness["input_status"] != "READY_FULL_SUMSTATS" or not path.is_file() or not qc_path.is_file():
        raise SystemExit(f"ERROR: dense source not ready: {trait}")
    metadata, steps = parse_qc(qc_path)
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        observed = next(csv.reader(handle, delimiter="\t"))
    required = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
    if observed != required:
        raise SystemExit(f"ERROR: dense schema drifted for {trait}: {observed}")
    if metadata.get("configured_build") != "hg19" or metadata.get("output_build") != "hg19":
        raise SystemExit(f"ERROR: Track B analysis build is not hg19 for {trait}")
    duplicates = find_step(steps, "duplicate SNP ID")
    multiallelic = find_step(steps, "non-SNP / indel / multi-base allele")
    strand = find_step(steps, "strand-ambiguous")
    maf = find_step(steps, "invalid FRQ or MAF")
    mhc = find_step(steps, "MHC chr6")
    info_steps = [(key, value) for key, value in steps.items() if key.startswith("INFO ")]
    if len(info_steps) != 1:
        raise SystemExit(f"ERROR: INFO audit missing for {trait}")
    info_key, info_value = info_steps[0]
    info_absent = "absent" in info_key
    fine_status = (
        "PASS_DENSE_WITH_SOURCE_LEVEL_INFO_LIMITATION"
        if info_absent else "PASS_DENSE_HARMONIZED"
    )
    case_control = (
        f"{panel['ncase']}/{panel['ncontrol']}"
        if panel["type"] == "binary" else "NA/NA"
    )
    sample_step = next(
        (key for key in steps if key.startswith("sample-size mode:")), None
    )
    if sample_step is None:
        raise SystemExit(f"ERROR: sample-size convention absent for {trait}")
    limitation = (
        "Release lacks row-wise INFO; source-level upstream INFO QC is documented but cannot support variant-level INFO sensitivity. "
        if info_absent else ""
    )
    limitation += (
        "Harmonized dense file excludes ambiguous SNPs, MAF<=0.01, MHC, duplicates, and non-SNP/multibase alleles; it is not a HapMap3-only file but is filtered rather than an unmodified raw release."
    )
    return {
        "pair_ids": TRAITS[trait][0],
        "trait_id": trait,
        "pair_role": TRAITS[trait][1],
        "source_gwas_id": panel["source_id"],
        "source": panel["source_note"],
        "publication": f"PMID:{panel['pmid']};DOI:{panel['doi']}",
        "ancestry": panel["ancestry"],
        "source_build": panel["build"],
        "analysis_build": metadata["output_build"],
        "dense_file": str(path),
        "dense_file_sha256": sha256(path),
        "compressed_size_bytes": str(path.stat().st_size),
        "rows_in": metadata["rows_in"],
        "rows_out": metadata["rows_out"],
        "pct_retained": metadata["pct_retained"],
        "observed_schema": ",".join(observed),
        "schema_status": schema["schema_status"],
        "rsid_status": "REQUIRED_RSIDS;MISSING_OR_NON_RSID_REMOVED",
        "effect_allele_status": f"A1_CANONICAL_EFFECT_ALLELE;source={schema['effect_allele']}",
        "other_allele_status": f"A2_CANONICAL_OTHER_ALLELE;source={schema['other_allele']}",
        "effect_status": f"BETA_SIGNED;source={schema['effect']}:{schema['effect_convention']}",
        "SE_status": f"PRESENT_POSITIVE;source={schema['standard_error']}",
        "P_status": f"PRESENT_IN_0_TO_1;source={schema['p_value']}",
        "EAF_status": f"FRQ_EFFECT_ALLELE_FREQUENCY;source={schema['eaf']}",
        "N_status": f"PRESENT;source={schema['sample_size']}",
        "case_control_counts": case_control,
        "INFO_status": (
            f"ABSENT_IN_RELEASE_SOURCE_LEVEL_ONLY;dropped={info_value[0]}"
            if info_absent else f"ROW_WISE_FILTERED_TO_0.6_TO_1;dropped={info_value[0]}"
        ),
        "liftover_status": "NOT_REQUIRED_SOURCE_AND_OUTPUT_GRCH37_HG19",
        "duplicate_variant_handling": f"REMOVED;dropped={duplicates[0]}",
        "multiallelic_handling": f"NON_SNP_INDEL_MULTIBASE_REMOVED;dropped={multiallelic[0]}",
        "strand_ambiguity_handling": f"AT_CG_REMOVED;dropped={strand[0]}",
        "MAF_handling": f"MAF_GT_0.01;dropped={maf[0]}",
        "MHC_handling": f"CHR6_25_34MB_REMOVED;dropped={mhc[0]}",
        "sample_size_convention": sample_step.removeprefix("sample-size mode: "),
        "fine_mapping_readiness": fine_status,
        "overall_QC": fine_status,
        "major_limitation": limitation,
    }


def tsv_text(rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def build() -> tuple[str, str]:
    pair_rows = read_tsv(PAIR_MANIFEST)
    if {(r["pair_id"], r["sleep_trait"], r["external_trait"]) for r in pair_rows} != {
        ("A", "snoring", "parental_lifespan"),
        ("B", "insomnia", "adhd"),
        ("CONTROL", "insomnia", "frailty"),
    }:
        raise SystemExit("ERROR: pair manifest drifted")
    panel = {row["trait_id"]: row for row in read_tsv(PANEL)}
    schemas = {row["source_id"]: row for row in read_tsv(SCHEMAS)}
    readiness = {row["trait_id"]: row for row in read_tsv(READINESS)}
    output_rows = [one(panel, schemas, readiness, trait) for trait in TRAITS]
    table = tsv_text(output_rows)
    lock = {
        "schema_version": 1,
        "pair_manifest_sha256": sha256(PAIR_MANIFEST),
        "readiness_table_sha256": sha256(READINESS),
        "dense_qc_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "trait_count": len(output_rows),
        "unique_dense_trait_family": list(TRAITS),
        "all_are_dense_not_hapmap3_only": True,
        "all_analysis_builds": ["hg19"],
        "rowwise_info_limitations": [
            row["trait_id"] for row in output_rows
            if row["INFO_status"].startswith("ABSENT_IN_RELEASE")
        ],
        "fine_mapping_gate_policy": "PASS_WITH_LIMITATION_IS_NOT_MECHANISTIC_EVIDENCE;LOCUS_LEVEL_LD_AND_VARIANT_QC_REMAIN_REQUIRED",
    }
    return table, json.dumps(lock, indent=2, sort_keys=True) + "\n"


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build()
    paths = [OUT, LOCK]
    if args.verify:
        for path, value in zip(paths, expected):
            if not path.is_file() or path.read_text(encoding="utf-8") != value:
                raise SystemExit(f"ERROR: dense Track B QC missing or drifted: {path}")
        print("verified Track B dense input QC: 5 traits across 3 frozen pairs")
        return
    for path, value in zip(paths, expected):
        atomic_text(path, value)
    print("wrote Track B dense input QC: 5 traits across 3 frozen pairs")


if __name__ == "__main__":
    main()
