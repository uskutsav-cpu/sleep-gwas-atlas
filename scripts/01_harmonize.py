#!/usr/bin/env python3
"""Harmonize one EUR GRCh37/hg19 GWAS summary-statistics file for LDSC.

This is the Phase 0 gatekeeper. It deliberately stops when it cannot verify a
non-negotiable assumption rather than producing a plausible-looking output.
It never infers ancestry or silently changes genome build. A registered hg38
source may opt into the checksum-pinned UCSC hg38-to-hg19 chain; other sources
requiring identity completion may opt into the pinned, provenance-checked
GRCh37 HapMap3 mapper. Both paths fail closed on missing or ambiguous matches.

Output schema (tab-separated, gzipped)::

    SNP  CHR  BP  A1  A2  FRQ  BETA  SE  P  N

``A1`` is the source effect allele and ``BETA`` is on the source effect scale
(log-odds for OR inputs). ``munge_sumstats.py --merge-alleles`` performs the
final HapMap3 allele compatibility check.

Usage::

    python3 scripts/01_harmonize.py --trait insomnia \
      --infile data/raw/insomnia.txt.gz --outdir data/harmonized
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys

import numpy as np
import pandas as pd

from variant_map import (
    STRATEGIES as VARIANT_MAP_STRATEGIES,
    VariantMapError,
    coordinate_match,
    load_variant_map,
    rsid_match,
)
from liftover_chain import LiftoverError, load_chain, reverse_complement


# Every threshold is either taken from CDG3 or required to make LDSC input
# auditable. Do not relax one to make a file pass.
INFO_MIN = 0.6  # [CDG3]
MAF_MIN = 0.01  # [CDG3]
N_EFF_MIN_FRACTION = 0.50  # [CDG3]
MHC_CHR = 6
MHC_START = 25_000_000
MHC_END = 34_000_000
TARGET_BUILD = "hg19"
AMBIGUOUS_ALLELES = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
VALID_ALLELES = {"A", "C", "G", "T"}
RSID = re.compile(r"rs[0-9]+$", re.IGNORECASE)

# Extend this table when a source has a documented alternate header. Do not
# silently guess a column based on its position.
ALIASES = {
    "SNP": ["snp", "rsid", "rs_id", "rsids", "markername", "variant_id", "id", "marker"],
    "CHR": ["chr", "chrom", "chromosome", "#chrom", "hg19chr"],
    "BP": [
        "bp", "pos", "position", "base_pair_location", "bp_hg19", "pos_hg19",
        # Graham et al. 2021 GLGC European releases explicitly label this
        # build-37 coordinate POS_b37 in their official README and headers.
        "pos_b37",
    ],
    "A1": [
        "a1", "effect_allele", "effec_allele", "ea", "allele1",
        "tested_allele", "alt",
    ],
    "A2": ["a2", "a0", "other_allele", "nea", "allele0", "allele2", "non_effect_allele", "ref"],
    "FRQ": [
        "frq", "freq", "eaf", "effect_allele_frequency", "maf", "a1freq",
        "freq1", "freq_tested_allele_in_hrs", "pooled_alt_af", "fcon",
        "frq_u_186843", "af_alt",
    ],
    "BETA": [
        "beta", "effect", "b", "log_odds", "logor", "effect_size",
        # Documented headers in the public UKB sleep releases registered in
        # config/public_gwas_sources.tsv (Dashti et al. 2019).
        "beta_sleepduration", "beta_shortsleep", "beta_longsleep",
        # Timmers et al. 2019 parental-survival release: a1 effect allele,
        # a0 reference allele, beta1 log-hazard protection ratio.
        "beta1",
    ],
    "OR": ["or", "odds_ratio", "oddsratio"],
    "SE": [
        "se", "standard_error", "stderr", "sebeta", "log_odds_se", "logor_se",
        # Howard et al. 2019 public MDD release: SE of LogOR for A1.
        "stderrlogor",
        "se_sleepduration", "se_shortsleep", "se_longsleep",
    ],
    "P": [
        "p", "pval", "pvalue", "p_value", "p_bolt_lmm", "p-value", "p.value",
        "p_sleepduration", "p_shortsleep", "p_longsleep",
    ],
    "N": ["n", "n_total", "samplesize", "n_complete_samples"],
    "N_EFF": ["neff", "n_eff", "effective_sample_size"],
    # Deelen 2019 defines its literal Effective_N as
    # 2/(1/N_cases + 1/N_controls), i.e. half the conventional LDSC N_eff.
    # PGC bipolar uses the semantically identical literal NEFFDIV2.
    "N_EFF_HALF": ["neffdiv2", "effective_n"],
    "NCASE": ["ncase", "n_cas", "n_cases", "cases", "ncas", "nca"],
    "NCONTROL": ["ncontrol", "n_con", "n_controls", "controls", "ncon", "nco"],
    "INFO": [
        "info", "imputation_info", "rsq", "r2", "imp_quality", "impinfo",
        # Schumacher et al. 2018 PRACTICAL EUR meta-analysis (GCST006085).
        "oncoarray_imputation_r2",
    ],
    "BUILD": ["build", "genome_build", "assembly"],
}


def log(message):
    print(f"  {message}", flush=True)


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def normalise_build(value):
    value = str(value).strip().lower()
    aliases = {
        "hg19": "hg19", "grch37": "hg19", "grch37/hg19": "hg19",
        "b37": "hg19", "hg38": "hg38", "grch38": "hg38",
        "grch38/hg38": "hg38", "b38": "hg38",
    }
    return aliases.get(value, value)


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def effective_n(ncase, ncontrol):
    """Effective N, 4 / (1/ncase + 1/ncontrol), used by CDG3."""
    return 4.0 / (1.0 / float(ncase) + 1.0 / float(ncontrol))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sniff_separator(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        metadata_lines = 0
        for header in handle:
            if header.startswith("##"):
                metadata_lines += 1
                continue
            break
        else:
            fail("source contains metadata but no tabular header")
    if header.count("\t") >= 2:
        return "\t", "tab", metadata_lines
    if header.count(",") >= 2:
        return ",", "comma", metadata_lines
    if len(header.split()) >= 3:
        return r"\s+", "whitespace", metadata_lines
    fail("could not determine a delimiter from the header")


def map_columns(columns):
    lower = {str(column).lower().strip(): column for column in columns}
    mapped = {}
    for standard, aliases in ALIASES.items():
        for alias in aliases:
            if alias in lower:
                mapped[standard] = lower[alias]
                break
    return mapped


def apply_schema_sample_size_mapping(config_path, trait, header_columns, mapped):
    """Apply an explicit source-schema effective-N conversion, if registered."""
    schema_path = os.path.join(os.path.dirname(os.path.abspath(config_path)), "gwas_schemas.tsv")
    if not os.path.isfile(schema_path):
        return mapped, "not supplied"
    schema = pd.read_csv(schema_path, sep="\t", dtype=str)
    rows = schema.loc[schema["trait_id"] == trait]
    if len(rows) > 1:
        fail(f"expected at most one source-schema row for trait '{trait}', found {len(rows)}")
    if len(rows) == 0:
        return mapped, "not supplied"
    raw_mapping = rows.iloc[0].get("sample_size", "")
    mapping = "" if pd.isna(raw_mapping) else str(raw_mapping).strip()
    match = re.fullmatch(r"DERIVED_N_EFF=2\*([A-Za-z0-9_.-]+)", mapping)
    if not match:
        return mapped, mapping or "not supplied"

    literal = match.group(1)
    matches = [column for column in header_columns if str(column).strip().casefold() == literal.casefold()]
    if len(matches) != 1:
        fail(
            f"schema requires {mapping}, but the source header contains "
            f"{len(matches)} columns named {literal!r}"
        )
    actual = matches[0]
    mapped = {standard: original for standard, original in mapped.items() if original != actual}
    mapped["N_EFF_HALF"] = actual
    return mapped, mapping


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True)
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--infile", required=True)
    parser.add_argument("--outdir", default="data/harmonized")
    parser.add_argument(
        "--source-build",
        help="Build declared by the downloaded file. Must agree with config; use this when the file header has no build column.",
    )
    parser.add_argument(
        "--variant-map",
        help="Pinned GRCh37 HapMap3 map; requires its .provenance.json companion.",
    )
    parser.add_argument(
        "--variant-map-strategy",
        choices=sorted(VARIANT_MAP_STRATEGIES),
        help="BY_COORD_ALLELES assigns rsIDs; BY_RSID_ALLELES assigns GRCh37 coordinates.",
    )
    parser.add_argument(
        "--expected-variant-map-sha256",
        help="Registered SHA-256 for the variant map (required by production entry points).",
    )
    parser.add_argument(
        "--expected-variant-map-bytes",
        type=int,
        help="Registered byte count for the variant map (required by production entry points).",
    )
    parser.add_argument(
        "--liftover-chain",
        help="Pinned UCSC hg38-to-hg19 chain for explicit point liftover.",
    )
    parser.add_argument(
        "--expected-liftover-chain-sha256",
        help="Registered SHA-256 for the liftover chain.",
    )
    parser.add_argument(
        "--expected-liftover-chain-bytes",
        type=int,
        help="Registered byte count for the liftover chain.",
    )
    parser.add_argument(
        "--prefilter-provenance",
        help="Provenance JSON from a registered bounded-memory source prefilter.",
    )
    args = parser.parse_args()

    if bool(args.variant_map) != bool(args.variant_map_strategy):
        fail("--variant-map and --variant-map-strategy must be supplied together")
    if (args.expected_variant_map_sha256 or args.expected_variant_map_bytes is not None) and not args.variant_map:
        fail("registered variant-map hash/bytes require --variant-map")
    if (
        args.expected_liftover_chain_sha256
        or args.expected_liftover_chain_bytes is not None
    ) and not args.liftover_chain:
        fail("registered liftover-chain hash/bytes require --liftover-chain")

    if not os.path.isfile(args.infile):
        fail(f"input file not found: {args.infile}")

    prefilter_provenance = None
    if args.prefilter_provenance:
        if not os.path.isfile(args.prefilter_provenance):
            fail(f"prefilter provenance not found: {args.prefilter_provenance}")
        with open(args.prefilter_provenance, encoding="utf-8") as handle:
            prefilter_provenance = json.load(handle)
        required_prefilter_fields = {
            "strategy", "input_path", "input_bytes", "input_sha256",
            "allowlist_path", "allowlist_count", "allowlist_sha256",
            "source_rows", "retained_rows", "output_path", "output_bytes",
            "output_sha256",
        }
        missing_prefilter_fields = required_prefilter_fields.difference(
            prefilter_provenance
        )
        if missing_prefilter_fields:
            fail(
                "prefilter provenance is missing fields: "
                f"{sorted(missing_prefilter_fields)}"
            )
        if prefilter_provenance["strategy"] != "HAPMAP3_RSID_ALLOWLIST":
            fail("unsupported prefilter provenance strategy")
        if int(prefilter_provenance["output_bytes"]) != os.path.getsize(args.infile):
            fail("prefilter output byte count does not match its provenance")
        if prefilter_provenance["output_sha256"].lower() != sha256(args.infile):
            fail("prefilter output SHA-256 does not match its provenance")
        if int(prefilter_provenance["retained_rows"]) > int(
            prefilter_provenance["source_rows"]
        ):
            fail("prefilter retained row count exceeds source row count")

    config = pd.read_csv(args.config, sep="\t", dtype=str)
    trait_rows = config.loc[config["trait_id"] == args.trait]
    if len(trait_rows) != 1:
        fail(f"expected exactly one config row for trait '{args.trait}', found {len(trait_rows)}")
    metadata = trait_rows.iloc[0]

    configured_build = normalise_build(metadata.get("build", ""))
    rsid_mapping_sets_build = args.variant_map_strategy == "BY_RSID_ALLELES"
    liftover_sets_build = bool(args.liftover_chain)
    valid_build_decision = (
        (liftover_sets_build and configured_build == "hg38")
        or (not liftover_sets_build and configured_build == TARGET_BUILD)
        or (
            not liftover_sets_build
            and rsid_mapping_sets_build
            and configured_build in {"", "unresolved"}
        )
    )
    if not valid_build_decision:
        fail(
            f"trait '{args.trait}' has build {metadata.get('build')!r}, which does not "
            "match the explicit hg19, coordinate-free-rsID, or hg38-liftover path."
        )
    if liftover_sets_build:
        if not args.source_build or normalise_build(args.source_build) != "hg38":
            fail("--liftover-chain requires an explicit --source-build hg38/GRCh38")
    elif args.source_build and normalise_build(args.source_build) != TARGET_BUILD:
        fail(f"--source-build {args.source_build!r} is not hg19/GRCh37; refusing silent liftover")

    separator, separator_name, metadata_lines = sniff_separator(args.infile)
    # Let pandas retain native numeric columns. Coercing every one of ~10M
    # rows to Python strings creates a multi-gigabyte working set before the
    # explicit validation below; SNP and allele fields are still converted to
    # normalized strings, and every numerical field is still passed through
    # pd.to_numeric before any QC decision. Chunked dtype inference also keeps
    # peak memory bounded for the largest public releases; explicit coercion
    # below, rather than inferred dtype, remains authoritative for every QC
    # decision.
    read_kwargs = {"low_memory": True, "comment": None, "skiprows": metadata_lines}
    if separator == r"\s+":
        read_kwargs.pop("low_memory")
        read_kwargs.update({"sep": separator, "engine": "python"})
    else:
        read_kwargs["sep"] = separator
    # Validate the source schema before parsing millions of records. This is
    # deliberately separate from the full read so a missing documented alias
    # fails in seconds rather than after allocating the complete table.
    header = pd.read_csv(args.infile, nrows=0, **read_kwargs)
    columns = map_columns(header.columns)
    columns, sample_size_schema_mapping = apply_schema_sample_size_mapping(
        args.config, args.trait, header.columns, columns
    )
    required = ["A1", "A2", "SE", "P"]
    if args.variant_map_strategy == "BY_COORD_ALLELES":
        if "CHR" in columns and "BP" in columns:
            required += ["CHR", "BP"]
        else:
            required += ["SNP"]
    elif args.variant_map_strategy == "BY_RSID_ALLELES":
        required += ["SNP"]
    else:
        required += ["SNP", "CHR", "BP"]
    if args.liftover_chain:
        required = list(dict.fromkeys(required + ["SNP", "CHR", "BP"]))
    missing = [column for column in required if column not in columns]
    if missing:
        fail(
            f"could not map required columns {missing}. File columns: {list(header.columns)}. "
            "Add a documented alias to ALIASES rather than renaming by position."
        )
    if "BETA" not in columns and "OR" not in columns:
        fail("neither BETA nor OR is present; an effect size is required")

    # A source can expose both a beta and a redundant odds ratio. BETA is the
    # documented preferred effect below, so do not allocate the unused OR
    # column. More generally, loading only registered mapped columns prevents
    # wide source-specific annotation fields from exhausting memory.
    if "BETA" in columns:
        columns.pop("OR", None)
    selected_source_columns = list(dict.fromkeys(columns.values()))
    read_kwargs["usecols"] = selected_source_columns
    raw = pd.read_csv(args.infile, **read_kwargs)
    n_input = len(raw)
    if n_input == 0:
        fail("input file has a header but no data rows")
    log(f"delimiter detected: {separator_name}; read {n_input:,} rows and {len(raw.columns)} columns")

    log(f"column map: {columns}")

    if "BUILD" in columns:
        file_builds = {
            normalise_build(value) for value in raw[columns["BUILD"]].dropna().unique()
            if str(value).strip()
        }
        expected_file_build = "hg38" if liftover_sets_build else TARGET_BUILD
        if file_builds and file_builds != {expected_file_build}:
            fail(
                f"file build values {sorted(file_builds)} do not all equal "
                f"{expected_file_build}; inspect the source before proceeding"
            )

    # All selected source columns map one-to-one to canonical fields. Renaming
    # in place avoids holding a second full dataframe during harmonization.
    source_to_standard = {original: standard for standard, original in columns.items()}
    if len(source_to_standard) != len(columns):
        fail("multiple canonical fields map to one source column; inspect the registered aliases")
    raw.rename(columns=source_to_standard, inplace=True)
    out = raw

    coordinate_label_parsed = False
    if args.variant_map_strategy == "BY_COORD_ALLELES" and (
        "CHR" not in out or "BP" not in out
    ):
        if "CHR" in out or "BP" in out:
            fail("coordinate mapping found only one of CHR/BP; inspect the source schema")
        labels = out["SNP"].astype("string").str.strip()
        parsed = labels.str.extract(r"^(?:chr)?([0-9]+)[:_]([0-9]+)(?:[:_]|$)")
        out["CHR"] = parsed[0]
        out["BP"] = parsed[1]
        coordinate_label_parsed = True
        log("parsed CHR/BP from the documented coordinate-based variant label")

    if "BETA" in out:
        out["BETA"] = pd.to_numeric(out["BETA"], errors="coerce")
    else:
        odds_ratio = pd.to_numeric(out["OR"], errors="coerce")
        out["BETA"] = np.log(odds_ratio.where(odds_ratio > 0))
        log("converted OR to log(OR)")

    for column in ["CHR", "BP", "FRQ", "SE", "P", "N", "N_EFF", "N_EFF_HALF", "NCASE", "NCONTROL", "INFO"]:
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    for column in ["A1", "A2"]:
        out[column] = out[column].astype("string").str.strip()
    if "SNP" in out:
        out["SNP"] = out["SNP"].astype("string").str.strip()
    else:
        out["SNP"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out["SNP"] = out["SNP"].str.lower()
    out["A1"] = out["A1"].str.upper()
    out["A2"] = out["A2"].str.upper()
    if "CHR" not in out:
        out["CHR"] = np.nan
    if "BP" not in out:
        out["BP"] = np.nan
    out["CHR"] = out["CHR"].astype("string").str.replace("chr", "", case=False, regex=False)
    out["CHR"] = pd.to_numeric(out["CHR"], errors="coerce")

    steps = []
    if coordinate_label_parsed:
        steps.append(("parsed CHR/BP from coordinate-based variant label", 0, len(out)))

    def drop(mask, reason):
        nonlocal out
        mask = mask.fillna(True) if isinstance(mask, pd.Series) else mask
        removed = int(mask.sum())
        if removed:
            out = out.loc[~mask].copy()
        steps.append((reason, removed, len(out)))

    variant_map_provenance = None
    liftover_provenance = None
    if args.liftover_chain:
        try:
            chain_index, liftover_provenance = load_chain(
                args.liftover_chain,
                args.expected_liftover_chain_sha256,
                args.expected_liftover_chain_bytes,
            )
        except LiftoverError as error:
            fail(str(error))
        statuses = []
        target_chromosomes = []
        target_positions = []
        target_a1 = []
        target_a2 = []
        reverse_strand_count = 0
        for chromosome, position, a1, a2 in zip(
            out["CHR"], out["BP"], out["A1"], out["A2"]
        ):
            status, mapped = chain_index.map_point(chromosome, position)
            statuses.append(status)
            if mapped is None:
                target_chromosomes.append(None)
                target_positions.append(None)
                target_a1.append(a1)
                target_a2.append(a2)
                continue
            target_chromosome, target_position, strand = mapped
            target_chromosomes.append(target_chromosome)
            target_positions.append(target_position)
            if strand == "-":
                target_a1.append(reverse_complement(a1))
                target_a2.append(reverse_complement(a2))
                reverse_strand_count += 1
            else:
                target_a1.append(a1)
                target_a2.append(a2)
        out["CHR"] = pd.Series(target_chromosomes, index=out.index, dtype="float64")
        out["BP"] = pd.Series(target_positions, index=out.index, dtype="float64")
        out["A1"] = pd.Series(target_a1, index=out.index, dtype="string")
        out["A2"] = pd.Series(target_a2, index=out.index, dtype="string")
        status_series = pd.Series(statuses, index=out.index, dtype="string")
        for status, reason in [
            ("invalid_source_coordinate", "invalid source coordinate before hg38-to-hg19 liftover"),
            ("unmapped", "not mapped by pinned UCSC hg38-to-hg19 chain"),
            ("ambiguous", "multiple target points in pinned UCSC hg38-to-hg19 chain"),
            ("non_autosomal_target", "liftover target is not an autosome"),
        ]:
            drop(status_series.loc[out.index] == status, reason)
        steps.append(
            (
                f"reverse-strand liftover mappings with allele complements: "
                f"{reverse_strand_count}",
                0,
                len(out),
            )
        )
        steps.append(
            (
                f"coordinates lifted hg38 to hg19; "
                f"chain_sha256={liftover_provenance['chain_sha256']}",
                0,
                len(out),
            )
        )

    if args.variant_map:
        try:
            variant_index, variant_map_provenance = load_variant_map(
                args.variant_map,
                args.variant_map_strategy,
                args.expected_variant_map_sha256,
                args.expected_variant_map_bytes,
            )
        except VariantMapError as error:
            fail(str(error))
        statuses = []
        mapped_snps = []
        mapped_chromosomes = []
        mapped_positions = []
        if args.variant_map_strategy == "BY_COORD_ALLELES":
            for chromosome, position, a1, a2 in zip(
                out["CHR"], out["BP"], out["A1"], out["A2"]
            ):
                status, snp = coordinate_match(
                    variant_index, chromosome, position, a1, a2
                )
                statuses.append(status)
                mapped_snps.append(snp)
            out["SNP"] = pd.Series(mapped_snps, index=out.index, dtype="string")
        else:
            for snp, chromosome, position, a1, a2 in zip(
                out["SNP"], out["CHR"], out["BP"], out["A1"], out["A2"]
            ):
                status, mapped = rsid_match(
                    variant_index, snp, chromosome, position, a1, a2
                )
                statuses.append(status)
                if mapped is None:
                    mapped_snps.append(None)
                    mapped_chromosomes.append(None)
                    mapped_positions.append(None)
                else:
                    mapped_snps.append(mapped[0])
                    mapped_chromosomes.append(mapped[1])
                    mapped_positions.append(mapped[2])
            out["SNP"] = pd.Series(mapped_snps, index=out.index, dtype="string")
            out["CHR"] = pd.Series(mapped_chromosomes, index=out.index, dtype="float64")
            out["BP"] = pd.Series(mapped_positions, index=out.index, dtype="float64")
        status_series = pd.Series(statuses, index=out.index, dtype="string")
        for status, reason in [
            ("unmatched_reference", "not present in pinned GRCh37 EUR HapMap3 map"),
            ("allele_conflict", "alleles conflict with pinned GRCh37 HapMap3 map"),
            ("coordinate_conflict", "source coordinate conflicts with pinned GRCh37 HapMap3 map"),
            ("ambiguous_reference", "coordinate/alleles map ambiguously to multiple rsIDs"),
        ]:
            drop(status_series.loc[out.index] == status, reason)
        steps.append(
            (
                f"variant identity assigned by {args.variant_map_strategy}; "
                f"map_sha256={variant_map_provenance['map_sha256']}",
                0,
                len(out),
            )
        )

    drop(out["SNP"].isna() | ~out["SNP"].map(lambda value: bool(RSID.fullmatch(str(value)))), "missing or non-rsID SNP (dbSNP resolution required)")
    drop(out["A1"].isna() | out["A2"].isna() | ~out["A1"].isin(VALID_ALLELES) | ~out["A2"].isin(VALID_ALLELES), "non-SNP / indel / multi-base allele")
    drop(out["A1"] == out["A2"], "A1 equals A2")
    ambiguous = pd.Series([(a1, a2) in AMBIGUOUS_ALLELES for a1, a2 in zip(out["A1"], out["A2"])], index=out.index)
    drop(ambiguous, "strand-ambiguous A/T or C/G [CDG3]")
    drop(out["P"].isna() | (out["P"] <= 0) | (out["P"] > 1), "P missing or outside (0, 1]")
    drop(out["BETA"].isna() | ~np.isfinite(out["BETA"]), "BETA missing or non-finite")
    drop(out["SE"].isna() | (out["SE"] <= 0) | ~np.isfinite(out["SE"]), "SE missing, non-positive, or non-finite")
    drop(out["CHR"].isna() | (out["CHR"] % 1 != 0) | ~out["CHR"].between(1, 22), "missing or non-autosomal CHR (EUR LDSC panel is autosomal)")
    drop(out["BP"].isna() | (out["BP"] <= 0) | (out["BP"] % 1 != 0), "missing or invalid BP")

    if "INFO" in out:
        drop(out["INFO"].notna() & ((out["INFO"] <= INFO_MIN) | (out["INFO"] > 1)), f"INFO outside ({INFO_MIN}, 1] [CDG3]")
    else:
        steps.append(("INFO column absent - source-level INFO QC must be documented", 0, len(out)))
    if "FRQ" in out:
        maf = out["FRQ"].where(out["FRQ"] <= 0.5, 1 - out["FRQ"])
        drop(out["FRQ"].notna() & ((out["FRQ"] < 0) | (out["FRQ"] > 1) | (maf <= MAF_MIN)), f"invalid FRQ or MAF <= {MAF_MIN} [CDG3]")
    else:
        steps.append(("FRQ column absent - source-level MAF QC must be documented", 0, len(out)))

    drop((out["CHR"] == MHC_CHR) & out["BP"].between(MHC_START, MHC_END), f"MHC chr6:{MHC_START}-{MHC_END} [CDG3]")
    drop(out["SNP"].duplicated(keep="first"), "duplicate SNP ID")

    binary = metadata["type"] == "binary"
    metadata_ncase = finite_number(metadata.get("ncase"))
    metadata_ncontrol = finite_number(metadata.get("ncontrol"))
    metadata_total = finite_number(metadata.get("n_total"))
    n_mode = ""
    n_reference = None
    n_reference_label = ""
    if binary and "N_EFF" in out:
        out["N"] = out["N_EFF"]
        n_mode = "per-SNP N_eff from source"
        n_reference = finite_number(out["N"].max())
        n_reference_label = "source maximum"
    elif binary and "N_EFF_HALF" in out:
        out["N"] = 2.0 * out["N_EFF_HALF"]
        n_mode = f"per-SNP N_eff derived as 2 * source {columns['N_EFF_HALF']}"
        n_reference = finite_number(out["N"].max())
        n_reference_label = "source maximum"
    elif binary and "NCASE" in out and "NCONTROL" in out:
        valid_counts = (out["NCASE"] > 0) & (out["NCONTROL"] > 0)
        out["N"] = np.where(valid_counts, 4.0 / (1.0 / out["NCASE"] + 1.0 / out["NCONTROL"]), np.nan)
        n_mode = "per-SNP N_eff derived from NCASE/NCONTROL"
        n_reference = finite_number(out["N"].max())
        n_reference_label = "source maximum"
    elif binary and metadata_ncase and metadata_ncontrol:
        n_reference = effective_n(metadata_ncase, metadata_ncontrol)
        n_reference_label = "configured effective N"
        out["N"] = n_reference
        n_mode = "constant N_eff derived from config ncase/ncontrol"
    elif not binary and "N" in out and out["N"].notna().any():
        n_reference = finite_number(out["N"].max())
        n_reference_label = "source maximum"
        n_mode = "per-SNP N from source"
    elif not binary and metadata_total:
        out["N"] = metadata_total
        n_reference = metadata_total
        n_reference_label = "configured total N"
        n_mode = "constant N from config n_total"
    else:
        fail(
            "could not establish an LDSC sample size. Binary traits need N_EFF, NCASE/NCONTROL, "
            "or cited config ncase/ncontrol; continuous traits need N or cited config n_total."
        )
    if n_reference is None or n_reference <= 0:
        fail("could not establish a positive sample-size reference for per-SNP QC")
    drop(
        out["N"].isna()
        | (out["N"] <= 0)
        | (out["N"] < N_EFF_MIN_FRACTION * n_reference),
        f"sample size below {N_EFF_MIN_FRACTION:.0%} of {n_reference_label} "
        f"({n_reference:,.1f}) [CDG3]",
    )
    steps.append((f"sample-size mode: {n_mode}", 0, len(out)))

    for column in ["FRQ"]:
        if column not in out:
            out[column] = np.nan
    out["CHR"] = out["CHR"].astype(int)
    out["BP"] = out["BP"].astype(int)
    out = out[["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]]

    os.makedirs(args.outdir, exist_ok=True)
    output_path = os.path.join(args.outdir, f"{args.trait}.harmonized.tsv.gz")
    report_path = os.path.join(args.outdir, f"{args.trait}.qc.txt")
    out.to_csv(output_path, sep="\t", index=False, na_rep="NA", compression="gzip")

    with open(report_path, "w", encoding="utf-8") as report:
        report.write(f"trait\t{args.trait}\n")
        report.write(f"source\t{metadata.get('source_note', '')}\n")
        report.write(f"configured_build\t{metadata.get('build', '')}\n")
        report.write(f"output_build\t{TARGET_BUILD}\n")
        report.write(f"source_build_argument\t{args.source_build or 'not supplied'}\n")
        report.write(f"sample_size_schema_mapping\t{sample_size_schema_mapping}\n")
        report.write(f"variant_map\t{os.path.abspath(args.variant_map) if args.variant_map else 'not supplied'}\n")
        report.write(f"variant_map_strategy\t{args.variant_map_strategy or 'not supplied'}\n")
        report.write(
            f"variant_map_sha256\t"
            f"{variant_map_provenance['map_sha256'] if variant_map_provenance else 'not supplied'}\n"
        )
        report.write(
            f"liftover_chain\t"
            f"{os.path.abspath(args.liftover_chain) if args.liftover_chain else 'not supplied'}\n"
        )
        report.write(
            f"liftover_chain_sha256\t"
            f"{liftover_provenance['chain_sha256'] if liftover_provenance else 'not supplied'}\n"
        )
        report.write(f"infile\t{os.path.abspath(args.infile)}\n")
        report.write(f"infile_sha256\t{sha256(args.infile)}\n")
        report.write(
            f"prefilter_strategy\t"
            f"{prefilter_provenance['strategy'] if prefilter_provenance else 'not supplied'}\n"
        )
        if prefilter_provenance:
            report.write(
                f"prefilter_source_rows\t{prefilter_provenance['source_rows']}\n"
            )
            report.write(
                f"prefilter_retained_rows\t{prefilter_provenance['retained_rows']}\n"
            )
            report.write(
                f"prefilter_source_sha256\t{prefilter_provenance['input_sha256']}\n"
            )
            report.write(
                f"prefilter_allowlist_sha256\t"
                f"{prefilter_provenance['allowlist_sha256']}\n"
            )
        report.write(f"rows_in\t{n_input}\n")
        report.write("\nstep\tdropped\tremaining\n")
        for reason, removed, remaining in steps:
            report.write(f"{reason}\t{removed}\t{remaining}\n")
        report.write(f"\nrows_out\t{len(out)}\n")
        report.write(f"pct_retained\t{100 * len(out) / n_input:.2f}\n")
        if prefilter_provenance:
            report.write(
                f"pct_retained_from_source\t"
                f"{100 * len(out) / int(prefilter_provenance['source_rows']):.2f}\n"
            )

    log(f"kept {len(out):,} / {n_input:,} ({100 * len(out) / n_input:.1f}%)")
    log(f"wrote {output_path}")
    log(f"wrote {report_path}")


if __name__ == "__main__":
    main()
