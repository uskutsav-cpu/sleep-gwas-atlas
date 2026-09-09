#!/usr/bin/env python3
"""Materialize the public BCAC 2020 overall-breast-cancer meta-analysis.

The official data dictionary defines ``SNP.iCOGs`` and ``SNP.Onco`` as
``rsid:position:non-effect-allele:effect-allele`` when an rsID is available;
otherwise they are coordinate labels.  This source-specific materializer keeps
only rows with an explicit rsID tag that agrees with the documented GRCh37
``var_name`` and meta-analysis alleles.  It never infers an rsID from a
coordinate-only row.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

from panel_guard import require_locked_traits
import re
import tempfile


ARCHIVE_BYTES = 4_764_280_363
ARCHIVE_SHA256 = "4f21dfc33d52753d496cd33756beed6e78e31fd6849fbd5755c406ad2c64e9c8"
HEADER = [
    "var_name", "Effect.Gwas", "Baseline.Gwas", "Freq.Gwas", "beta.Gwas",
    "SE.Gwas", "P.value.Gwas", "SNP.iCOGs", "chr.iCOGs", "Position.iCOGs",
    "Effect.iCOGs", "Baseline.iCOGs", "EAFcontrols.iCOGs", "EAFcases.iCOGs",
    "r2.iCOGs", "NumCalled.iCOGs", "beta.iCOGs", "SE.iCOGs", "SE_LRT.iCOGs",
    "chi2.iCOGs", "LRT.iCOGs", "P1df_risk_chi.iCOGs", "P1df_risk_LRT.iCOGs",
    "SNP.Onco", "chr.Onco", "Position.Onco", "Effect.Onco", "Baseline.Onco",
    "EAFcontrols.Onco", "EAFcases.Onco", "r2.Onco", "NumCalled.Onco",
    "beta.Onco", "SE.Onco", "SE_LRT.Onco", "chi2.Onco", "LRT.Onco",
    "P1df_risk_chi.Onco", "P1df_risk_LRT.Onco", "Effect.Meta", "Baseline.Meta",
    "Beta.meta", "var.meta", "sdE.meta", "chi.meta", "p.meta", "dir.meta",
]
REQUIRED_ROW_FIELDS = [
    "var_name", "SNP.iCOGs", "SNP.Onco", "Effect.Meta", "Baseline.Meta",
    "Beta.meta", "sdE.meta", "p.meta",
]
AUTOSOMES = {str(chromosome) for chromosome in range(1, 23)}
RSID_TAG = re.compile(r"(rs[0-9]+):(\d+):([ACGT]):([ACGT])$", re.IGNORECASE)
VARIANT = re.compile(r"([0-9]+|X|Y):?[_]([0-9]+)[_]([ACGT])[_]([ACGT])$", re.IGNORECASE)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_archive(source: Path) -> None:
    if not source.is_file():
        fail(f"source file not found: {source}")
    if source.stat().st_size != ARCHIVE_BYTES:
        fail(f"source byte count is {source.stat().st_size}, expected {ARCHIVE_BYTES}")
    if sha256(source) != ARCHIVE_SHA256:
        fail("source SHA-256 does not match the registered BCAC release")


def parse_variant(value: str, line_number: int) -> tuple[str, str, str, str] | None:
    match = VARIANT.fullmatch(value.strip())
    if match is None:
        return None
    chromosome, position, non_effect, effect = match.groups()
    return chromosome.upper(), position, non_effect.upper(), effect.upper()


def explicit_rsid(
    icogs: str, onco: str, chromosome: str, position: str,
    non_effect: str, effect: str, line_number: int,
) -> str | None:
    valid = []
    for label, value in (("SNP.iCOGs", icogs), ("SNP.Onco", onco)):
        match = RSID_TAG.fullmatch(value.strip())
        if match is None:
            continue
        rsid, tag_position, tag_non_effect, tag_effect = match.groups()
        candidate = (tag_position, tag_non_effect.upper(), tag_effect.upper())
        if candidate != (position, non_effect, effect):
            fail(
                f"{label} tag conflicts with var_name/meta alleles at line {line_number}: "
                f"{value!r} versus {chromosome}_{position}_{non_effect}_{effect}"
            )
        valid.append(rsid.lower())
    if not valid:
        return None
    if len(set(valid)) != 1:
        fail(f"conflicting explicit rsIDs at line {line_number}: {valid}")
    return valid[0]


def injected_header_fragment(row: dict[str, str | list[str] | None]) -> bool:
    """Recognize the one incomplete record followed by a repeated source header.

    The verified BCAC archive contains its normal header at line 1 and one
    additional header injected after 16 fields of a record at line 1,018.  The
    fragment lacks all remaining iCOGS/Onco/meta values and therefore cannot
    yield an association statistic.  Check its fixed adjacent header tokens so
    a genuinely malformed association cannot be silently discarded.
    """
    return (
        isinstance(row.get("beta.iCOGs"), str)
        and str(row["beta.iCOGs"]).endswith("var_name")
        and row.get("SE.iCOGs") == "Effect.Gwas"
        and row.get("P1df_risk_LRT.iCOGs") == "P.value.Gwas"
        and row.get("Effect.Meta") == "SNP.Onco"
        and row.get("p.meta") == "EAFcases.Onco"
    )


def audit_source(source: Path, output: Path, progress_every: int) -> None:
    """Count source-row states without creating a materialized GWAS file."""
    counts = {
        "source_rows": 0,
        "materializable_explicit_rsid_rows": 0,
        "coordinate_only_rows": 0,
        "nonautosomal_rows": 0,
        "injected_header_fragments": 0,
        "incomplete_rows": 0,
        "malformed_var_name_rows": 0,
        "meta_allele_conflict_rows": 0,
        "rsid_tag_conflict_rows": 0,
        "conflicting_dual_rsid_rows": 0,
    }
    examples: dict[str, list[dict[str, object]]] = {
        key: [] for key in counts if key not in {"source_rows", "materializable_explicit_rsid_rows"}
    }

    def record(kind: str, line_number: int, row: dict[str, str | list[str] | None]) -> None:
        counts[kind] += 1
        if len(examples[kind]) < 3:
            examples[kind].append({
                "line": line_number,
                "var_name": row.get("var_name"),
                "SNP.iCOGs": row.get("SNP.iCOGs"),
                "SNP.Onco": row.get("SNP.Onco"),
            })

    with source.open("rt", encoding="utf-8", newline="") as source_handle:
        reader = csv.DictReader(source_handle, delimiter=" ", skipinitialspace=True)
        if reader.fieldnames != HEADER:
            fail(f"unexpected BCAC header: {reader.fieldnames}")
        for line_number, row in enumerate(reader, start=2):
            counts["source_rows"] += 1
            if counts["source_rows"] % progress_every == 0:
                print(f"  audited {counts['source_rows']:,} source rows", flush=True)
            if injected_header_fragment(row):
                record("injected_header_fragments", line_number, row)
                continue
            if any(row.get(field) is None for field in REQUIRED_ROW_FIELDS):
                record("incomplete_rows", line_number, row)
                continue
            parsed = parse_variant(str(row["var_name"]), line_number)
            if parsed is None:
                record("malformed_var_name_rows", line_number, row)
                continue
            chromosome, position, non_effect, effect = parsed
            if chromosome not in AUTOSOMES:
                record("nonautosomal_rows", line_number, row)
                continue
            if str(row["Effect.Meta"]).strip().upper() != effect or str(row["Baseline.Meta"]).strip().upper() != non_effect:
                record("meta_allele_conflict_rows", line_number, row)
                continue
            rsids = []
            tag_conflict = False
            for value in (str(row["SNP.iCOGs"]), str(row["SNP.Onco"])):
                match = RSID_TAG.fullmatch(value.strip())
                if match is None:
                    continue
                rsid, tag_position, tag_non_effect, tag_effect = match.groups()
                if (tag_position, tag_non_effect.upper(), tag_effect.upper()) != (position, non_effect, effect):
                    tag_conflict = True
                    break
                rsids.append(rsid.lower())
            if tag_conflict:
                record("rsid_tag_conflict_rows", line_number, row)
            elif not rsids:
                record("coordinate_only_rows", line_number, row)
            elif len(set(rsids)) != 1:
                record("conflicting_dual_rsid_rows", line_number, row)
            else:
                counts["materializable_explicit_rsid_rows"] += 1

    payload = {"archive_sha256": ARCHIVE_SHA256, "counts": counts, "examples": examples}
    temporary = output.with_name(f".{output.name}.partial")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, output)
    print(f"Wrote source audit: {output}")
    print(json.dumps(counts, sort_keys=True))


def materialize(source: Path, output: Path, progress_every: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".partial", dir=output.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    total = kept = non_rsid = malformed = injected_header = incomplete = nonautosomal = 0
    try:
        with source.open("rt", encoding="utf-8", newline="") as source_handle, gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as output_handle:
            reader = csv.DictReader(source_handle, delimiter=" ", skipinitialspace=True)
            if reader.fieldnames != HEADER:
                fail(f"unexpected BCAC header: {reader.fieldnames}")
            writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"])
            for line_number, row in enumerate(reader, start=2):
                total += 1
                if total % progress_every == 0:
                    print(
                        f"  processed {total:,} source rows; retained {kept:,} explicit rsID rows",
                        flush=True,
                    )
                if injected_header_fragment(row):
                    injected_header += 1
                    continue
                if any(row.get(field) is None for field in REQUIRED_ROW_FIELDS):
                    # The source occasionally contains a partial constituent
                    # record with no meta-analysis fields. It cannot yield a
                    # beta/SE/P triplet and is counted rather than padded.
                    incomplete += 1
                    continue
                parsed = parse_variant(row["var_name"], line_number)
                if parsed is None:
                    malformed += 1
                    continue
                chromosome, position, non_effect, effect = parsed
                if chromosome not in AUTOSOMES:
                    nonautosomal += 1
                    continue
                if row["Effect.Meta"].strip().upper() != effect or row["Baseline.Meta"].strip().upper() != non_effect:
                    fail(
                        f"meta effect/non-effect alleles conflict with var_name at line {line_number}: "
                        f"{row['var_name']!r}"
                    )
                rsid = explicit_rsid(
                    row["SNP.iCOGs"], row["SNP.Onco"], chromosome, position,
                    non_effect, effect, line_number,
                )
                if rsid is None:
                    non_rsid += 1
                    continue
                writer.writerow([
                    rsid, chromosome, position, effect, non_effect,
                    row["Beta.meta"], row["sdE.meta"], row["p.meta"],
                ])
                kept += 1
        with gzip.open(temporary, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise

    print(f"Materialized {output}")
    print(f"  source rows: {total:,}")
    print(f"  retained explicit rsID rows: {kept:,}")
    print(f"  skipped coordinate-only rows: {non_rsid:,}")
    print(f"  skipped non-autosomal rows: {nonautosomal:,}")
    print(f"  skipped injected source-header fragments: {injected_header:,}")
    print(f"  skipped incomplete source rows without a full meta record: {incomplete:,}")
    print(f"  skipped malformed var_name rows: {malformed:,}")
    print(f"  output SHA-256: {sha256(output)}")


def main() -> None:
    require_locked_traits({"breast_cancer": "bcac_2020_overall_breast_cancer"})
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--audit-out", type=Path)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()
    if sum((args.verify_only, args.out is not None, args.audit_out is not None)) != 1:
        fail("supply exactly one of --verify-only, --out, or --audit-out")
    if args.progress_every <= 0:
        fail("--progress-every must be positive")
    verify_archive(args.source)
    print("Verified BCAC 2020 overall-breast-cancer source archive")
    if args.audit_out is not None:
        audit_source(args.source, args.audit_out, args.progress_every)
    elif not args.verify_only:
        materialize(args.source, args.out, args.progress_every)


if __name__ == "__main__":
    main()
