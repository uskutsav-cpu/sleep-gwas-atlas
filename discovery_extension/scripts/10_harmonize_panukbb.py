#!/usr/bin/env python3
"""Stream one locked Pan-UKB trait into the canonical extension LDSC schema."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import re
import tempfile
from collections import Counter
from pathlib import Path


VALID = {"A", "C", "G", "T"}
AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}
RSID = re.compile(r"rs[0-9]+")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: str, label: str, counts: Counter[str]) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        counts[f"invalid_{label}"] += 1
        return None
    if not math.isfinite(number):
        counts[f"invalid_{label}"] += 1
        return None
    return number


def load_reference(path: Path) -> dict[tuple[int, int, str, str], dict[str, str]]:
    reference: dict[tuple[int, int, str, str], dict[str, str]] = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "CHR", "BP", "REF", "ALT", "INFO"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"ERROR: reference lacks fields: {sorted(missing)}")
        for row in reader:
            key = (int(row["CHR"]), int(row["BP"]), row["REF"], row["ALT"])
            if key in reference:
                raise SystemExit(f"ERROR: duplicate exact key in variant reference: {key}")
            reference[key] = row
    return reference


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-id", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument(
        "--panel", type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--reference", type=Path,
        default=Path("discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz"),
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--qc-out", type=Path)
    args = parser.parse_args()

    panel = {row["extension_trait_id"]: row for row in read_tsv(args.panel)}
    if args.trait_id not in panel:
        raise SystemExit(f"ERROR: trait is not in the locked panel: {args.trait_id}")
    trait = panel[args.trait_id]
    if args.source.name != trait["source_filename"]:
        raise SystemExit("ERROR: source basename differs from the locked source filename")
    for path in (args.source, args.reference):
        if not path.is_file():
            raise SystemExit(f"ERROR: required input is missing: {path}")
    out = args.out or Path(f"discovery_extension/data/harmonized/{args.trait_id}.txt.gz")
    qc_out = args.qc_out or Path(f"discovery_extension/results/qc/harmonization/{args.trait_id}.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    qc_out.parent.mkdir(parents=True, exist_ok=True)

    reference = load_reference(args.reference)
    counts: Counter[str] = Counter()
    seen_rsids: set[str] = set()
    binary = trait["binary_or_continuous"] == "binary"
    cases = int(trait["cases"]) if binary else None
    controls = int(trait["controls"]) if binary else None
    sample_n = 4.0 / (1.0 / cases + 1.0 / controls) if binary else float(trait["sample_size"])

    handle_fd, temporary_name = tempfile.mkstemp(prefix=out.name + ".", suffix=".tmp", dir=out.parent)
    os.close(handle_fd)
    temporary = Path(temporary_name)
    try:
        with gzip.open(args.source, "rt", newline="", encoding="utf-8") as source, gzip.open(
            temporary, "wt", newline="", encoding="utf-8"
        ) as destination:
            reader = csv.DictReader(source, delimiter="\t")
            common = {"chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"}
            frequency = {"af_cases_EUR", "af_controls_EUR"} if binary else {"af_EUR"}
            missing = (common | frequency) - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(f"ERROR: phenotype source lacks fields: {sorted(missing)}")
            writer = csv.DictWriter(
                destination, delimiter="\t",
                fieldnames=["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N", "INFO"],
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                counts["input_rows"] += 1
                if row["low_confidence_EUR"].strip().lower() not in {"false", "0"}:
                    counts["low_confidence_EUR"] += 1
                    continue
                try:
                    chromosome, position = int(row["chr"]), int(row["pos"])
                except ValueError:
                    counts["invalid_coordinate"] += 1
                    continue
                ref, alt = row["ref"].upper(), row["alt"].upper()
                identity = reference.get((chromosome, position, ref, alt))
                if identity is None:
                    counts["not_exact_pinned_hapmap3_identity"] += 1
                    continue
                if chromosome < 1 or chromosome > 22:
                    counts["non_autosomal"] += 1
                    continue
                if ref not in VALID or alt not in VALID or ref == alt:
                    counts["non_snp_or_invalid_alleles"] += 1
                    continue
                if frozenset((ref, alt)) in AMBIGUOUS:
                    counts["strand_ambiguous"] += 1
                    continue
                if chromosome == 6 and 25_000_000 <= position <= 34_000_000:
                    counts["extended_MHC"] += 1
                    continue
                beta = finite(row["beta_EUR"], "beta", counts)
                se = finite(row["se_EUR"], "se", counts)
                neglog = finite(row["neglog10_pval_EUR"], "neglog10_p", counts)
                info = finite(identity["INFO"], "info", counts)
                if None in (beta, se, neglog, info):
                    continue
                assert beta is not None and se is not None and neglog is not None and info is not None
                if se <= 0 or neglog < 0:
                    counts["invalid_se_or_p"] += 1
                    continue
                if not 0.9 <= info <= 1:
                    counts["info_below_0_9_or_above_1"] += 1
                    continue
                if binary:
                    af_cases = finite(row["af_cases_EUR"], "af_cases", counts)
                    af_controls = finite(row["af_controls_EUR"], "af_controls", counts)
                    if af_cases is None or af_controls is None:
                        continue
                    frequency_value = (af_cases * cases + af_controls * controls) / (cases + controls)
                else:
                    frequency_value = finite(row["af_EUR"], "frequency", counts)
                    if frequency_value is None:
                        continue
                if not 0 <= frequency_value <= 1 or min(frequency_value, 1 - frequency_value) < 0.01:
                    counts["maf_below_0_01_or_invalid_frequency"] += 1
                    continue
                rsid = identity["SNP"]
                if not RSID.fullmatch(rsid):
                    raise SystemExit(f"ERROR: reference yielded invalid rsID: {rsid}")
                if rsid in seen_rsids:
                    raise SystemExit(f"ERROR: duplicate rsID in harmonized source: {rsid}")
                seen_rsids.add(rsid)
                p_value = max(10.0 ** (-neglog), 1e-300)
                writer.writerow({
                    "SNP": rsid, "CHR": chromosome, "BP": position,
                    "A1": alt, "A2": ref, "FRQ": f"{frequency_value:.12g}",
                    "BETA": f"{beta:.12g}", "SE": f"{se:.12g}",
                    "P": f"{p_value:.12g}", "N": f"{sample_n:.12g}",
                    "INFO": f"{info:.12g}",
                })
                counts["output_rows"] += 1
        temporary.replace(out)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise

    with qc_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["metric", "count", "notes"], lineterminator="\n")
        writer.writeheader()
        for metric, count in sorted(counts.items()):
            writer.writerow({"metric": metric, "count": count, "notes": "prospective locked filter"})
        writer.writerow({"metric": "output_sha256", "count": sha256(out), "notes": str(out)})
        writer.writerow({"metric": "reference_sha256", "count": sha256(args.reference), "notes": str(args.reference)})
        writer.writerow({"metric": "source_sha256", "count": sha256(args.source), "notes": str(args.source)})

    if counts["output_rows"] == 0:
        raise SystemExit("ERROR: no variants survived harmonization")
    print(
        f"PANUKBB_HARMONIZED trait={args.trait_id} rows={counts['output_rows']} "
        f"sha256={sha256(out)}"
    )


if __name__ == "__main__":
    main()
