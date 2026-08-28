#!/usr/bin/env python3
"""Stream an exact source gzip to the registered HapMap3 rsID family.

This is a bounded-memory pre-harmonization operation for very large releases.
Selected rows are copied byte-for-byte; the only new bytes are a deterministic
gzip wrapper. A provenance JSON binds the source, allowlist, output, and exact
row counts. Malformed physical rows fail closed.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys


RSID = re.compile(br"rs[0-9]+$", re.IGNORECASE)


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_allowlist(path):
    variants = set()
    with open(path, "rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            token = line.rstrip(b"\r\n").split(b"\t", 1)[0].strip().lower()
            if line_number == 1 and token == b"snp":
                continue
            if not RSID.fullmatch(token):
                fail(f"invalid rsID in allowlist {path} at line {line_number}")
            variants.add(token)
    if not variants:
        fail(f"allowlist contains no rsIDs: {path}")
    return variants


def field_at(line, index):
    start = 0
    for _ in range(index):
        start = line.find(b"\t", start)
        if start < 0:
            return None
        start += 1
    end = line.find(b"\t", start)
    if end < 0:
        end = len(line.rstrip(b"\r\n"))
    return line[start:end].strip().lower()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--snp-column", required=True)
    parser.add_argument("--allowlist", required=True)
    parser.add_argument("--expected-input-bytes", type=int, required=True)
    parser.add_argument("--expected-input-sha256", required=True)
    args = parser.parse_args()

    for path in [args.input, args.allowlist]:
        if not os.path.isfile(path):
            fail(f"required input not found: {path}")
    input_bytes = os.path.getsize(args.input)
    if input_bytes != args.expected_input_bytes:
        fail(
            f"input byte count mismatch: expected {args.expected_input_bytes}, "
            f"found {input_bytes}"
        )
    input_sha = sha256(args.input)
    if input_sha.lower() != args.expected_input_sha256.lower():
        fail(
            f"input SHA-256 mismatch: expected {args.expected_input_sha256}, "
            f"found {input_sha}"
        )

    allowlist = load_allowlist(args.allowlist)
    allowlist_sha = sha256(args.allowlist)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.provenance)), exist_ok=True)
    if os.path.isfile(args.output) and os.path.isfile(args.provenance):
        try:
            with open(args.provenance, encoding="utf-8") as handle:
                existing = json.load(handle)
            reusable = (
                existing.get("strategy") == "HAPMAP3_RSID_ALLOWLIST"
                and existing.get("input_sha256", "").lower() == input_sha.lower()
                and int(existing.get("input_bytes", -1)) == input_bytes
                and existing.get("allowlist_sha256", "").lower()
                == allowlist_sha.lower()
                and existing.get("snp_column") == args.snp_column
                and int(existing.get("output_bytes", -1))
                == os.path.getsize(args.output)
                and existing.get("output_sha256", "").lower()
                == sha256(args.output).lower()
            )
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            reusable = False
        if reusable:
            print(
                f"Reusing verified prefilter with {existing['retained_rows']:,} / "
                f"{existing['source_rows']:,} rows; output SHA-256 "
                f"{existing['output_sha256']}"
            )
            return
    temp_output = f"{args.output}.tmp.{os.getpid()}"
    temp_provenance = f"{args.provenance}.tmp.{os.getpid()}"
    source_rows = 0
    retained_rows = 0

    try:
        with gzip.open(args.input, "rb") as source:
            header = source.readline()
            if not header:
                fail("source gzip is empty")
            header_fields = header.rstrip(b"\r\n").split(b"\t")
            matches = [
                index for index, field in enumerate(header_fields)
                if field.decode("utf-8").strip() == args.snp_column
            ]
            if len(matches) != 1:
                fail(
                    f"expected exactly one source column {args.snp_column!r}; "
                    f"found {len(matches)} in {header_fields!r}"
                )
            snp_index = matches[0]
            expected_fields = len(header_fields)
            with open(temp_output, "wb") as raw_output:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw_output, mtime=0
                ) as output:
                    output.write(header)
                    for line_number, line in enumerate(source, start=2):
                        source_rows += 1
                        if line.count(b"\t") != expected_fields - 1:
                            fail(
                                f"malformed source row {line_number}: expected "
                                f"{expected_fields} tab-separated fields"
                            )
                        snp = field_at(line, snp_index)
                        if snp in allowlist:
                            output.write(line)
                            retained_rows += 1
                        if source_rows % 5_000_000 == 0:
                            print(
                                f"Scanned {source_rows:,} rows; retained "
                                f"{retained_rows:,}",
                                flush=True,
                            )
        if source_rows == 0:
            fail("source contains a header but no data rows")
        if retained_rows == 0:
            fail("no source rows matched the registered HapMap3 allowlist")
        os.replace(temp_output, args.output)
        output_bytes = os.path.getsize(args.output)
        output_sha = sha256(args.output)
        provenance = {
            "strategy": "HAPMAP3_RSID_ALLOWLIST",
            "input_path": os.path.abspath(args.input),
            "input_bytes": input_bytes,
            "input_sha256": input_sha,
            "snp_column": args.snp_column,
            "allowlist_path": os.path.abspath(args.allowlist),
            "allowlist_count": len(allowlist),
            "allowlist_sha256": allowlist_sha,
            "source_rows": source_rows,
            "retained_rows": retained_rows,
            "output_path": os.path.abspath(args.output),
            "output_bytes": output_bytes,
            "output_sha256": output_sha,
        }
        with open(temp_provenance, "w", encoding="utf-8") as handle:
            json.dump(provenance, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temp_provenance, args.provenance)
    finally:
        for path in [temp_output, temp_provenance]:
            if os.path.exists(path):
                os.unlink(path)

    print(
        f"Prefiltered {retained_rows:,} / {source_rows:,} rows; "
        f"output SHA-256 {output_sha}"
    )


if __name__ == "__main__":
    main()
