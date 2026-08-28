#!/usr/bin/env python3
"""Stream the pinned EUR archive into the locked allele-matched HapMap3 PLINK subset."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def append_suffix(path: Path, suffix: str) -> Path:
    return Path(str(path) + suffix)


def safe_relative(root: Path, value: object, label: str) -> Path:
    relative = Path(str(value))
    if relative.is_absolute() or ".." in relative.parts:
        fail(f"unsafe {label} path")
    return root / relative


def require_pin(path: Path, expected_bytes: object, expected_hash: object, label: str) -> None:
    if (
        not path.is_file() or path.stat().st_size != expected_bytes
        or sha256(path) != expected_hash
    ):
        fail(f"{label} is absent or differs from its exact byte/SHA-256 pin")


def load_instrument_universe(path: Path, expected_rows: int) -> dict[str, frozenset[str]]:
    targets: dict[str, frozenset[str]] = {}
    with path.open(encoding="ascii", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["SNP", "A1", "A2"]:
            fail("HapMap3 instrument universe has an unexpected schema")
        for line_number, row in enumerate(reader, start=2):
            snp, a1, a2 = row["SNP"], row["A1"].upper(), row["A2"].upper()
            if (
                not snp.startswith("rs") or not snp[2:].isdigit() or snp in targets
                or a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"}
                or a1 == a2
            ):
                fail(f"invalid or duplicate HapMap3 instrument at line {line_number}")
            targets[snp] = frozenset((a1, a2))
    if len(targets) + 1 != expected_rows:
        fail("HapMap3 instrument-universe row count differs from its pin")
    return targets


def expected_reference_members(manifest: dict[str, object]) -> dict[str, dict[str, object]]:
    expected = {row["name"]: row for row in manifest.get("reference_members", [])}
    required = {"README", "g1000_eur.bed", "g1000_eur.bim", "g1000_eur.fam", "g1000_eur.synonyms"}
    if set(expected) != required:
        fail("EUR reference member family differs from its pin")
    return expected


def stream_bim(
    archive: zipfile.ZipFile,
    member: zipfile.ZipInfo,
    expected: dict[str, object],
    targets: dict[str, frozenset[str]],
    output: Path,
    expected_reference_variants: int,
    minimum_selected_variants: int = 100000,
) -> tuple[list[int], dict[str, int], str]:
    digest = hashlib.sha256()
    selected_indices: list[int] = []
    observed_targets: set[str] = set()
    selected_targets: set[str] = set()
    reference_rows = 0
    with archive.open(member) as source, output.open("wb") as sink:
        for index, raw_line in enumerate(source):
            digest.update(raw_line)
            values = raw_line.split()
            if len(values) != 6:
                fail(f"EUR BIM row {index + 1} does not have six fields")
            try:
                snp = values[1].decode("ascii")
                alleles = frozenset((values[4].decode("ascii").upper(), values[5].decode("ascii").upper()))
            except UnicodeError as exc:
                fail(f"EUR BIM contains a non-ASCII identifier or allele: {exc}")
            if snp in targets:
                observed_targets.add(snp)
                if alleles == targets[snp]:
                    if snp in selected_targets:
                        fail(f"EUR BIM contains a duplicate exact target: {snp}")
                    selected_targets.add(snp)
                    selected_indices.append(index)
                    sink.write(raw_line)
            reference_rows += 1
    if member.file_size != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
        fail("streamed EUR BIM differs from its exact member pin")
    if reference_rows != expected_reference_variants:
        fail("EUR BIM variant count differs from the BED-derived exact count")
    counters = {
        "reference_variants": reference_rows,
        "target_variants": len(targets),
        "reference_id_observed": len(observed_targets),
        "allele_matched_variants": len(selected_targets),
        "allele_mismatch_variants": len(observed_targets - selected_targets),
        "reference_absent_variants": len(targets) - len(observed_targets),
    }
    if len(selected_indices) < minimum_selected_variants:
        fail(f"fewer than {minimum_selected_variants:,} exact allele-matched HapMap3 variants remain")
    return selected_indices, counters, digest.hexdigest()


def stream_bed(
    archive: zipfile.ZipFile,
    member: zipfile.ZipInfo,
    expected: dict[str, object],
    selected_indices: list[int],
    output: Path,
    record_bytes: int,
    expected_reference_variants: int,
) -> str:
    digest = hashlib.sha256()
    selected_cursor = 0
    base_index = 0
    records_per_chunk = 65536
    with archive.open(member) as source, output.open("wb") as sink:
        header = source.read(3)
        digest.update(header)
        if header != b"\x6c\x1b\x01":
            fail("EUR BED is not a PLINK 1 SNP-major binary file")
        sink.write(header)
        while True:
            block = source.read(record_bytes * records_per_chunk)
            if not block:
                break
            digest.update(block)
            if len(block) % record_bytes:
                fail("EUR BED stream ends with a partial variant record")
            record_count = len(block) // record_bytes
            end_index = base_index + record_count
            while selected_cursor < len(selected_indices) and selected_indices[selected_cursor] < end_index:
                index = selected_indices[selected_cursor]
                if index < base_index:
                    fail("selected BED indices are not strictly increasing")
                offset = (index - base_index) * record_bytes
                sink.write(block[offset:offset + record_bytes])
                selected_cursor += 1
            base_index = end_index
    if member.file_size != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
        fail("streamed EUR BED differs from its exact member pin")
    if base_index != expected_reference_variants or selected_cursor != len(selected_indices):
        fail("EUR BED record count or selected-index coverage differs from the BIM")
    if output.stat().st_size != 3 + len(selected_indices) * record_bytes:
        fail("materialized causal BED has an unexpected byte count")
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--out-prefix", default="work/causal_reference/g1000_eur_hm3")
    parser.add_argument("--acknowledge-streamed-reference", action="store_true")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    causal = policy["causal_inference"]
    runtime_path = safe_relative(root, causal["component_manifest"], "causal runtime manifest")
    if not runtime_path.is_file() or sha256(runtime_path) != causal["component_manifest_sha256"]:
        fail("causal runtime manifest is absent or differs from the policy pin")
    runtime = json.loads(runtime_path.read_text(encoding="utf-8"))
    reference = runtime["ld_reference"]
    source_manifest_path = safe_relative(root, reference["source_manifest"], "reference source manifest")
    if not source_manifest_path.is_file() or sha256(source_manifest_path) != reference["source_manifest_sha256"]:
        fail("EUR reference source manifest differs from the causal runtime pin")
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    members = expected_reference_members(source_manifest)
    record_bytes = int(reference["bed_record_bytes"])
    if (members["g1000_eur.bed"]["bytes"] - 3) % record_bytes:
        fail("EUR BED member pin is not an integer SNP-major record family")
    expected_reference_variants = (members["g1000_eur.bed"]["bytes"] - 3) // record_bytes
    if expected_reference_variants != 22665064:
        fail("EUR BED-derived reference variant count differs from its release pin")
    archive_path = safe_relative(root, reference["archive_path"], "EUR archive")
    require_pin(archive_path, reference["archive_bytes"], reference["archive_sha256"], "EUR reference archive")
    instrument_path = safe_relative(root, reference["instrument_universe_path"], "instrument universe")
    require_pin(
        instrument_path, reference["instrument_universe_bytes"],
        reference["instrument_universe_sha256"], "HapMap3 instrument universe",
    )
    plink_component = next(
        row for row in runtime["components"]
        if row["component_id"] == "PLINK_1_9_STABLE_UNIVERSAL_BINARY"
    )
    plink = safe_relative(root, plink_component["path"], "PLINK binary")
    require_pin(plink, plink_component["bytes"], plink_component["sha256"], "PLINK binary")
    output_prefix = safe_relative(root, args.out_prefix, "causal reference output")
    final_paths = [
        output_prefix.with_suffix(suffix) for suffix in (".bed", ".bim", ".fam", ".provenance.json")
    ] + [
        append_suffix(output_prefix.with_name(output_prefix.name + ".validation"), suffix)
        for suffix in (".frq.counts", ".log", ".nosex")
    ]
    if any(path.exists() for path in final_paths):
        fail("causal reference output is immutable and one or more target files already exist")
    if not args.execute:
        print(
            "CAUSAL_REFERENCE_READY_TO_STREAM "
            f"archive_bytes={archive_path.stat().st_size} uncompressed_bytes={source_manifest['reference_archive_uncompressed_bytes']} "
            f"targets={reference['instrument_universe_rows_including_header'] - 1} out={output_prefix}"
        )
        return 0
    if not args.acknowledge_streamed_reference:
        fail("streaming 3.60 GB from the compressed EUR reference requires --acknowledge-streamed-reference")
    free_bytes = shutil.disk_usage(output_prefix.parent if output_prefix.parent.exists() else root).free
    if free_bytes < 400 * 1024 * 1024:
        fail("at least 400 MiB free disk is required for the disk-bounded causal reference subset")

    targets = load_instrument_universe(
        instrument_path, int(reference["instrument_universe_rows_including_header"]),
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    temporary_directory = Path(tempfile.mkdtemp(prefix=".causal-reference-", dir=output_prefix.parent))
    temporary_prefix = temporary_directory / output_prefix.name
    try:
        with zipfile.ZipFile(archive_path) as archive:
            infos = {info.filename: info for info in archive.infolist() if not info.is_dir()}
            if set(infos) != set(members):
                fail("EUR archive member family differs from its pin")
            for name, info in infos.items():
                if info.file_size != members[name]["bytes"]:
                    fail(f"EUR archive member byte count differs from its pin: {name}")
            selected_indices, counters, bim_member_hash = stream_bim(
                archive, infos["g1000_eur.bim"], members["g1000_eur.bim"], targets,
                temporary_prefix.with_suffix(".bim"), expected_reference_variants,
            )
            with archive.open(infos["g1000_eur.fam"]) as handle:
                fam = handle.read()
            if (
                len(fam) != members["g1000_eur.fam"]["bytes"]
                or hashlib.sha256(fam).hexdigest() != members["g1000_eur.fam"]["sha256"]
                or len(fam.splitlines()) != reference["sample_count"]
            ):
                fail("EUR FAM differs from its exact member or sample-count pin")
            temporary_prefix.with_suffix(".fam").write_bytes(fam)
            bed_member_hash = stream_bed(
                archive, infos["g1000_eur.bed"], members["g1000_eur.bed"], selected_indices,
                temporary_prefix.with_suffix(".bed"), record_bytes, expected_reference_variants,
            )

        validation_prefix = temporary_prefix.with_name(temporary_prefix.name + ".validation")
        validation = subprocess.run(
            [
                str(plink), "--bfile", str(temporary_prefix), "--freq", "counts",
                "--allow-no-sex", "--out", str(validation_prefix),
            ],
            capture_output=True, text=True, check=False,
        )
        if validation.returncode:
            fail("PLINK failed to read and frequency-check the streamed causal subset")
        frequency_path = append_suffix(validation_prefix, ".frq.counts")
        log_path = append_suffix(validation_prefix, ".log")
        if not frequency_path.is_file() or not log_path.is_file():
            fail("PLINK causal-subset validation did not produce frequency and log artifacts")
        with frequency_path.open(encoding="utf-8") as handle:
            frequency_rows = sum(1 for _ in handle) - 1
        if frequency_rows != len(selected_indices):
            fail("PLINK frequency validation row count differs from the causal subset")

        artifacts: dict[str, dict[str, object]] = {}
        generated = [
            temporary_prefix.with_suffix(".bed"), temporary_prefix.with_suffix(".bim"),
            temporary_prefix.with_suffix(".fam"), frequency_path, log_path,
        ]
        nosex_path = append_suffix(validation_prefix, ".nosex")
        if nosex_path.is_file():
            generated.append(nosex_path)
        for path in generated:
            artifacts[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
        provenance = {
            "schema_version": "sleep-atlas-causal-reference.1",
            "policy_sha256": sha256(policy_path),
            "causal_runtime_manifest_sha256": sha256(runtime_path),
            "source_manifest_sha256": sha256(source_manifest_path),
            "source_archive": str(archive_path.relative_to(root)),
            "source_archive_sha256": sha256(archive_path),
            "instrument_universe": str(instrument_path.relative_to(root)),
            "instrument_universe_sha256": sha256(instrument_path),
            "selection_rule": reference["instrument_universe_rule"],
            "sample_count": reference["sample_count"],
            "bed_record_bytes": reference["bed_record_bytes"],
            "source_bim_sha256": bim_member_hash,
            "source_bed_sha256": bed_member_hash,
            "counters": counters,
            "plink_version": validation.stdout.splitlines()[0] if validation.stdout else "UNKNOWN",
            "plink_returncode": validation.returncode,
            "frequency_rows": frequency_rows,
            "artifacts": artifacts,
            "script_sha256": sha256(Path(__file__).resolve()),
        }
        provenance_path = temporary_prefix.with_suffix(".provenance.json")
        provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        generated.append(provenance_path)
        for temporary in generated:
            destination = output_prefix.parent / temporary.name
            temporary.replace(destination)
    finally:
        shutil.rmtree(temporary_directory, ignore_errors=True)

    print(
        "CAUSAL_REFERENCE_READY "
        f"variants={counters['allele_matched_variants']} samples={reference['sample_count']} "
        f"bed_bytes={output_prefix.with_suffix('.bed').stat().st_size} out={output_prefix}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
