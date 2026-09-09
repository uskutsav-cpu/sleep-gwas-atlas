#!/usr/bin/env python3
"""Materialize the public Phelan et al. 2017 OCAC ovarian-GWAS release.

The original EBI archive has one tab-delimited table for each chromosome.  Its
small official data dictionary says that these tables describe analyses in
Europeans and that ``overall`` is the all-invasive histotype.  Although the
source labels the effect column ``overall_OR``, the same dictionary explicitly
defines it as a *log odds ratio*; this extractor therefore emits it as BETA,
not an odds ratio to be logged again.

Only source-supplied, explicit rsID labels and biallelic A/C/G/T SNPs are
retained.  It does not infer rsIDs, change build, or combine one of the
histotype-specific columns with the all-invasive phenotype.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
from pathlib import Path
import re
import tempfile
import zipfile

from panel_guard import require_locked_traits


ARCHIVE_BYTES = 3_774_309_689
MEMBERS = [f"Summary_chr{chromosome}.txt" for chromosome in range(1, 24)]
AUTOSOME_MEMBERS = MEMBERS[:22]
REQUIRED_COLUMNS = [
    "1000G_SNPname", "Chromosome", "Position", "Effect", "Baseline", "EAF",
    "R2_oncoarray", "overall_OR", "overall_SE", "overall_pvalue",
]
RSID_TAG = re.compile(r"^(rs[0-9]+):(\d+):([ACGT]+):([ACGT]+)$", re.IGNORECASE)
SNP_ALLELES = {"A", "C", "G", "T"}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_archive(source: Path, expected_sha256: str) -> None:
    if not source.is_file():
        fail(f"source file not found: {source}")
    if source.stat().st_size != ARCHIVE_BYTES:
        fail(
            f"source byte count is {source.stat().st_size}, expected {ARCHIVE_BYTES}"
        )
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        fail("expected SHA-256 must be a 64-character lowercase hexadecimal digest")
    if sha256(source) != expected_sha256:
        fail("source SHA-256 does not match the registered Phelan release")
    try:
        with zipfile.ZipFile(source) as archive:
            names = archive.namelist()
            missing = set(MEMBERS).difference(names)
            unexpected = {
                name for name in names
                if name.startswith("Summary_chr") and name not in MEMBERS
            }
            if missing or unexpected:
                fail(
                    "unexpected chromosome members: "
                    f"missing={sorted(missing)}, unexpected={sorted(unexpected)}"
                )
    except zipfile.BadZipFile as error:
        fail(f"invalid ZIP container: {error}")
    print("Verified Phelan 2017 OCAC ovarian-cancer source archive")


def materialize(source: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".partial", dir=output.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    total = kept = non_rsid = non_snp = malformed = tag_conflict = tag_complemented = 0
    try:
        with zipfile.ZipFile(source) as archive, gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as output_handle:
            writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "INFO",
            ])
            # The source's 23rd file is chromosome X.  It is verified as part
            # of the release container but deliberately never emitted to the
            # autosomal LDSC input.
            for chromosome, member in enumerate(AUTOSOME_MEMBERS, start=1):
                with archive.open(member, "r") as binary_handle, io.TextIOWrapper(
                    binary_handle, encoding="utf-8", newline=""
                ) as source_handle:
                    header = source_handle.readline()
                    delimiter = "\t" if header.count("\t") >= header.count(",") else ","
                    fieldnames = next(csv.reader([header], delimiter=delimiter), None)
                    if fieldnames is None or any(
                        column not in fieldnames for column in REQUIRED_COLUMNS
                    ):
                        fail(f"unexpected header in {member}: {fieldnames}")
                    reader = csv.DictReader(
                        source_handle, fieldnames=fieldnames, delimiter=delimiter
                    )
                    for line_number, row in enumerate(reader, start=2):
                        total += 1
                        if total % 1_000_000 == 0:
                            print(
                                f"  processed {total:,} rows; retained {kept:,} rsID SNPs",
                                flush=True,
                            )
                        row_fields = (
                            "1000G_SNPname", "Chromosome", "Position", "Effect", "Baseline",
                            "EAF", "R2_oncoarray", "overall_OR", "overall_SE", "overall_pvalue",
                        )
                        if any(row.get(column) is None or not row[column].strip() for column in row_fields):
                            malformed += 1
                            continue
                        if row["Chromosome"].strip() != str(chromosome):
                            fail(
                                f"chromosome/member disagreement at {member}:{line_number}: "
                                f"{row['Chromosome']!r}"
                            )
                        effect = row["Effect"].strip().upper()
                        baseline = row["Baseline"].strip().upper()
                        if not row["Position"].strip().isdigit():
                            malformed += 1
                            continue
                        if effect not in SNP_ALLELES or baseline not in SNP_ALLELES or effect == baseline:
                            non_snp += 1
                            continue
                        tag = RSID_TAG.fullmatch(row["1000G_SNPname"].strip())
                        if tag is None:
                            non_rsid += 1
                            continue
                        rsid, tagged_position, tagged_allele_one, tagged_allele_two = tag.groups()
                        tagged_alleles = {tagged_allele_one.upper(), tagged_allele_two.upper()}
                        row_alleles = {baseline, effect}
                        tag_is_complement = {
                            COMPLEMENT[allele] for allele in tagged_alleles
                        } == row_alleles
                        if tagged_position != row["Position"].strip() or (
                            tagged_alleles != row_alleles and not tag_is_complement
                        ):
                            # This public source occasionally has a 1000G tag
                            # incompatible with its reported association row.
                            # That row cannot safely provide an rsID, so retain
                            # neither the tag nor an inferred replacement.
                            tag_conflict += 1
                            continue
                        if tag_is_complement and tagged_alleles != row_alleles:
                            tag_complemented += 1
                        writer.writerow([
                            rsid.lower(), chromosome, row["Position"].strip(), effect, baseline,
                            row["EAF"].strip(), row["overall_OR"].strip(),
                            row["overall_SE"].strip(), row["overall_pvalue"].strip(),
                            row["R2_oncoarray"].strip(),
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
    print(f"  retained source-explicit rsID SNPs: {kept:,}")
    print(f"  skipped non-rsID rows: {non_rsid:,}")
    print(f"  skipped non-SNP rows: {non_snp:,}")
    print(f"  skipped incomplete rows: {malformed:,}")
    print(f"  skipped conflicting rsID-tag rows: {tag_conflict:,}")
    print(f"  source rsID tags on complementary strand: {tag_complemented:,}")
    print(f"  output SHA-256: {sha256(output)}")


def main() -> None:
    require_locked_traits({"ovarian_cancer": "phelan_2017_ovarian_cancer"})
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if args.verify_only == (args.out is not None):
        fail("supply exactly one of --verify-only or --out")
    verify_archive(args.source, args.expected_sha256)
    if not args.verify_only:
        materialize(args.source, args.out)


if __name__ == "__main__":
    main()
