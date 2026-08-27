#!/usr/bin/env python3
"""Download and materialize a registered OpenGWAS GWAS-VCF release.

OpenGWAS exposes stable dataset pages but generates short-lived download URLs.
The downloader therefore POSTs the registered dataset ID to the public trial
endpoint, selects the exact VCF object (never the index), and verifies the
complete registered byte count and SHA-256 before replacing the destination.

The materializer is intentionally strict about GWAS-VCF semantics.  It emits
the atlas canonical columns with ALT as the effect allele, converts LP from
``-log10(P)`` to P, and refuses a build, sample, FORMAT, or allele mismatch.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation, localcontext
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import ssl
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from panel_guard import require_locked_traits


DATASET_ID = "ieu-b-38"
SOURCE_ID = "evangelou_2018_sbp_opengwas"
TOTAL_VARIANTS = 7_088_083
SAMPLE_REQUIRED = {
    "ID": DATASET_ID,
    "TotalVariants": str(TOTAL_VARIANTS),
    "VariantsNotRead": "0",
    "HarmonisedVariants": str(TOTAL_VARIANTS),
    "VariantsNotHarmonised": "0",
    "StudyType": "Continuous",
}
FORMAT_REQUIRED = {"ES", "SE", "LP", "AF", "SS", "ID"}
VCF_HEADER = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    "FORMAT", DATASET_ID,
]
RSID = re.compile(r"rs[0-9]+$", re.IGNORECASE)
AUTOSOMES = {str(chromosome) for chromosome in range(1, 23)}
ALLELES = {"A", "C", "G", "T"}
P_FLOOR = Decimal("1e-300")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ssl_context() -> ssl.SSLContext:
    context = ssl.create_default_context()
    try:
        import certifi  # type: ignore[import-not-found]
    except ImportError:
        return context
    return ssl.create_default_context(cafile=certifi.where())


def json_urls(value: object) -> list[str]:
    """Collect URL strings from the endpoint without assuming response shape."""
    if isinstance(value, str):
        return [value] if value.startswith("https://") else []
    if isinstance(value, list):
        return [url for item in value for url in json_urls(item)]
    if isinstance(value, dict):
        return [url for item in value.values() for url in json_urls(item)]
    return []


def select_vcf_url(payload: object, dataset_id: str) -> str:
    expected_suffix = f"/{dataset_id}/{dataset_id}.vcf.gz"
    candidates = {
        url for url in json_urls(payload)
        if urlparse(url).scheme == "https"
        and urlparse(url).path.endswith(expected_suffix)
    }
    if len(candidates) != 1:
        fail(
            f"expected exactly one HTTPS VCF URL ending {expected_suffix!r}; "
            f"found {len(candidates)}"
        )
    return candidates.pop()


def request_download_url(endpoint: str, dataset_id: str) -> str:
    parsed = urlparse(endpoint)
    if parsed.scheme != "https" or not parsed.netloc:
        fail("OpenGWAS link endpoint must be an absolute HTTPS URL")
    body = json.dumps({"id": dataset_id, "type": "vcf"}).encode("utf-8")
    request = Request(
        endpoint,
        data=body,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "sleep-gwas-atlas/1.0",
        },
        method="POST",
    )
    with urlopen(request, timeout=60, context=ssl_context()) as response:
        if response.status != 200:
            fail(f"OpenGWAS link endpoint returned HTTP {response.status}")
        try:
            payload = json.load(response)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            fail(f"OpenGWAS link endpoint did not return valid JSON: {error}")
    return select_vcf_url(payload, dataset_id)


def verify_file(path: Path, expected_bytes: int, expected_sha256: str) -> None:
    if not path.is_file():
        fail(f"downloaded file not found: {path}")
    if path.stat().st_size != expected_bytes:
        fail(
            f"downloaded byte count is {path.stat().st_size}, expected {expected_bytes}"
        )
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        fail("expected SHA-256 must be 64 lowercase hexadecimal characters")
    if sha256(path) != expected_sha256:
        fail("downloaded SHA-256 does not match the registered OpenGWAS release")


def download(
    endpoint: str, dataset_id: str, output: Path,
    expected_bytes: int, expected_sha256: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file():
        try:
            verify_file(output, expected_bytes, expected_sha256)
        except SystemExit:
            pass
        else:
            print(f"Registered OpenGWAS archive already present: {output}")
            return

    vcf_url = request_download_url(endpoint, dataset_id)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".partial", dir=output.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        request = Request(vcf_url, headers={"User-Agent": "sleep-gwas-atlas/1.0"})
        digest = hashlib.sha256()
        downloaded = 0
        with urlopen(request, timeout=120, context=ssl_context()) as response, temporary.open("wb") as handle:
            if response.status != 200:
                fail(f"OpenGWAS VCF download returned HTTP {response.status}")
            for block in iter(lambda: response.read(1024 * 1024), b""):
                handle.write(block)
                digest.update(block)
                downloaded += len(block)
        if downloaded != expected_bytes:
            fail(f"downloaded {downloaded} bytes, expected {expected_bytes}")
        if digest.hexdigest() != expected_sha256:
            fail("downloaded SHA-256 does not match the registered OpenGWAS release")
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    print(f"Downloaded and verified OpenGWAS dataset {dataset_id}: {output}")


def angle_metadata(line: str, prefix: str) -> dict[str, str]:
    if not line.startswith(prefix) or not line.endswith(">"):
        fail(f"malformed VCF metadata line: {line[:120]!r}")
    body = line[len(prefix):-1]
    values: dict[str, str] = {}
    for item in body.split(","):
        if "=" not in item:
            fail(f"malformed VCF metadata item: {item!r}")
        key, value = item.split("=", 1)
        values[key] = value.strip('"')
    return values


def finite_decimal(value: str, label: str, line_number: int) -> Decimal:
    try:
        number = Decimal(value)
    except InvalidOperation:
        fail(f"invalid {label} at VCF line {line_number}: {value!r}")
    if not number.is_finite():
        fail(f"non-finite {label} at VCF line {line_number}: {value!r}")
    return number


def p_from_lp(value: str, line_number: int) -> tuple[str, bool]:
    lp = finite_decimal(value, "LP", line_number)
    if lp < 0:
        fail(f"negative LP at VCF line {line_number}: {value!r}")
    if lp > 300:
        return "1.000000000000000E-300", True
    with localcontext() as context:
        context.prec = 50
        p_value = Decimal(10) ** (-lp)
    return f"{p_value:.15E}", False


def parse_info(value: str, line_number: int) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in value.split(";"):
        if "=" not in item:
            continue
        key, field_value = item.split("=", 1)
        values[key] = field_value
    if "AF" not in values:
        fail(f"INFO/AF is missing at VCF line {line_number}")
    return values


def materialize_vcf(
    source: Path, output: Path, dataset_id: str = DATASET_ID,
    expected_variants: int = TOTAL_VARIANTS,
) -> dict[str, int]:
    if not source.is_file():
        fail(f"source file not found: {source}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".partial", dir=output.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    counts = {"source_rows": 0, "p_values_floored_at_1e-300": 0}
    sample_metadata: dict[str, str] | None = None
    format_ids: set[str] = set()
    build_declarations = 0
    header_seen = False
    try:
        with gzip.open(source, "rt", encoding="utf-8", newline="") as source_handle, temporary.open("wb") as binary_output:
            with gzip.GzipFile(
                filename="", mode="wb", compresslevel=1, mtime=0,
                fileobj=binary_output,
            ) as compressed_output, io.TextIOWrapper(
                compressed_output, encoding="utf-8", newline=""
            ) as text_output:
                writer = csv.writer(text_output, delimiter="\t", lineterminator="\n")
                writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"])
                for line_number, raw_line in enumerate(source_handle, start=1):
                    line = raw_line.rstrip("\r\n")
                    if line.startswith("##SAMPLE=<"):
                        if sample_metadata is not None:
                            fail("VCF declares more than one SAMPLE metadata record")
                        sample_metadata = angle_metadata(line, "##SAMPLE=<")
                        continue
                    if line.startswith("##FORMAT=<ID="):
                        match = re.match(r"##FORMAT=<ID=([^,>]+)", line)
                        if match is None:
                            fail(f"malformed FORMAT metadata line: {line[:120]!r}")
                        format_ids.add(match.group(1))
                        continue
                    if line.startswith("##contig=<"):
                        contig = angle_metadata(line, "##contig=<")
                        assembly = contig.get("assembly", "").lower()
                        if assembly not in {"hg19", "grch37"}:
                            fail(f"contig metadata does not declare HG19/GRCh37: {line[:160]!r}")
                        build_declarations += 1
                        continue
                    if line.startswith("##"):
                        continue
                    fields = line.split("\t")
                    if not header_seen:
                        if fields != VCF_HEADER[:-1] + [dataset_id]:
                            fail(f"unexpected GWAS-VCF header: {fields}")
                        if sample_metadata is None:
                            fail("VCF header appeared before SAMPLE metadata")
                        for key, expected in {
                            **SAMPLE_REQUIRED,
                            "ID": dataset_id,
                            "TotalVariants": str(expected_variants),
                            "HarmonisedVariants": str(expected_variants),
                        }.items():
                            if sample_metadata.get(key) != expected:
                                fail(
                                    f"SAMPLE/{key} is {sample_metadata.get(key)!r}, "
                                    f"expected {expected!r}"
                                )
                        missing_formats = FORMAT_REQUIRED.difference(format_ids)
                        if missing_formats:
                            fail(f"VCF metadata is missing FORMAT definitions: {sorted(missing_formats)}")
                        if build_declarations == 0:
                            fail("VCF has no HG19/GRCh37 contig declaration")
                        header_seen = True
                        continue
                    if len(fields) != 10:
                        fail(f"expected 10 VCF fields at line {line_number}, got {len(fields)}")
                    chrom, position, variant_id, ref, alt, _qual, row_filter, info, format_text, sample_text = fields
                    if chrom not in AUTOSOMES or not position.isdigit() or int(position) <= 0:
                        fail(f"invalid autosomal coordinate at VCF line {line_number}")
                    if not RSID.fullmatch(variant_id):
                        fail(f"non-rsID VCF ID at line {line_number}: {variant_id!r}")
                    ref = ref.upper()
                    alt = alt.upper()
                    if ref not in ALLELES or alt not in ALLELES or ref == alt:
                        fail(f"non-biallelic-SNP alleles at VCF line {line_number}: {ref}/{alt}")
                    if row_filter != "PASS":
                        fail(f"non-PASS variant at VCF line {line_number}: {row_filter!r}")
                    keys = format_text.split(":")
                    values = sample_text.split(":")
                    if len(keys) != len(values) or len(keys) != len(set(keys)):
                        fail(f"malformed FORMAT/sample fields at VCF line {line_number}")
                    sample = dict(zip(keys, values))
                    missing = FORMAT_REQUIRED.difference(sample)
                    if missing:
                        fail(f"missing FORMAT fields at VCF line {line_number}: {sorted(missing)}")
                    if sample["ID"].lower() != variant_id.lower():
                        fail(f"VCF ID and FORMAT/ID disagree at line {line_number}")
                    effect = finite_decimal(sample["ES"], "ES", line_number)
                    standard_error = finite_decimal(sample["SE"], "SE", line_number)
                    frequency = finite_decimal(sample["AF"], "AF", line_number)
                    sample_size = finite_decimal(sample["SS"], "SS", line_number)
                    if (
                        standard_error <= 0 or frequency < 0 or frequency > 1
                        or sample_size <= 0 or sample_size != sample_size.to_integral_value()
                    ):
                        fail(f"invalid SE, AF, or SS at VCF line {line_number}")
                    info_frequency = finite_decimal(parse_info(info, line_number)["AF"], "INFO/AF", line_number)
                    if info_frequency != frequency:
                        fail(f"INFO/AF and FORMAT/AF disagree at VCF line {line_number}")
                    p_value, floored = p_from_lp(sample["LP"], line_number)
                    counts["source_rows"] += 1
                    counts["p_values_floored_at_1e-300"] += int(floored)
                    writer.writerow([
                        variant_id.lower(), chrom, position, alt, ref, sample["AF"],
                        str(effect), str(standard_error), p_value, str(sample_size),
                    ])
                    if counts["source_rows"] % 1_000_000 == 0:
                        print(f"  materialized {counts['source_rows']:,} GWAS-VCF rows", flush=True)
        if not header_seen:
            fail("VCF header was not found")
        if counts["source_rows"] != expected_variants:
            fail(
                f"materialized {counts['source_rows']:,} variants, "
                f"expected {expected_variants:,}"
            )
        with gzip.open(temporary, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    print(
        f"Materialized {counts['source_rows']:,} variants; "
        f"floored {counts['p_values_floored_at_1e-300']:,} P values at 1e-300"
    )
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--download", action="store_true")
    mode.add_argument("--materialize", action="store_true")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--endpoint")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--expected-bytes", type=int)
    parser.add_argument("--expected-sha256")
    args = parser.parse_args()

    if args.dataset_id != DATASET_ID:
        fail(f"unsupported dataset ID: {args.dataset_id!r}")
    require_locked_traits({"sbp": SOURCE_ID})
    if args.download:
        if not args.endpoint or args.expected_bytes is None or not args.expected_sha256:
            fail("--download requires --endpoint, --expected-bytes, and --expected-sha256")
        if args.source is not None:
            fail("--source is not valid with --download")
        download(
            args.endpoint, args.dataset_id, args.out,
            args.expected_bytes, args.expected_sha256,
        )
    else:
        if args.source is None:
            fail("--materialize requires --source")
        if args.endpoint or args.expected_bytes is not None or args.expected_sha256:
            fail("download-only arguments are not valid with --materialize")
        materialize_vcf(args.source, args.out, args.dataset_id)


if __name__ == "__main__":
    main()
