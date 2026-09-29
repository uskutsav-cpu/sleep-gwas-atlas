#!/usr/bin/env python3
"""Validate and publish complete PLACO+ v3 pair collations without clobbering."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "extensions/brain6"
sys.path.insert(0, str(EXTENSION))
from brain6.artifacts import verify_artifact  # noqa: E402
from brain6.io import exclusive_lock  # noqa: E402

FIELDS = ["SNP", "CHR", "BP", "Z1", "Z2", "P_PLACO", "Q_WITHIN_PAIR",
          "Q_BONFERRONI_ACROSS_PAIRS", "status", "headline"]
MASTER_FIELDS = ["pair_id", "sleep_trait", "brain_disorder", "stage_status",
                 "total_rows", "tested_rows", "excluded_extreme_z", "numerical_failures",
                 "failure_rate", "across_track_p_threshold", "headline_variants",
                 "independent_loci_status", "output_path", "output_sha256",
                 "source_artifact_fingerprint", "source_receipt_sha256", "family_lock_sha256"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_variant_table(path: Path, expected_rows: int, threshold: float) -> dict:
    """Check schema, row identity, status accounting, P/q domains, and headline flags."""
    n = tested = excluded = failures = headline = 0
    snps: set[str] = set()
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != FIELDS:
            raise ValueError(f"Unexpected PLACO result columns: {reader.fieldnames}")
        for row in reader:
            n += 1
            snp = row["SNP"]
            if not snp or snp in snps:
                raise ValueError(f"Missing or duplicate SNP identifier: {snp!r}")
            snps.add(snp)
            status = row["status"]
            if status == "TESTED":
                tested += 1
                values = [float(row[k]) for k in ("P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS")]
                if not all(math.isfinite(v) and 0 <= v <= 1 for v in values):
                    raise ValueError(f"Invalid P/q value at {snp}")
                is_headline = values[0] < threshold
                if row["headline"] != str(is_headline):
                    raise ValueError(f"Headline flag disagrees with frozen threshold at {snp}")
                headline += is_headline
            elif status == "EXCLUDED_EXTREME_Z":
                excluded += 1
                if any(row[k] != "NA" for k in ("P_PLACO", "Q_WITHIN_PAIR", "Q_BONFERRONI_ACROSS_PAIRS")) or row["headline"] != "False":
                    raise ValueError(f"Excluded row has inferential values: {snp}")
            elif status == "NUMERICAL_FAILURE":
                failures += 1
                if row["headline"] != "False":
                    raise ValueError(f"Failed row cannot be headline: {snp}")
            else:
                raise ValueError(f"Unknown PLACO row status: {status}")
    if n != expected_rows:
        raise ValueError(f"Expected {expected_rows} PLACO rows, found {n}")
    return {"total_rows": n, "tested_rows": tested, "excluded_extreme_z": excluded,
            "numerical_failures": failures, "headline_variants": headline}


def atomic_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    source_hash = sha256(source)
    if target.exists():
        if sha256(target) != source_hash:
            raise FileExistsError(f"Refusing to replace different published output: {target}")
        return
    fd, name = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        shutil.copyfile(source, temporary)
        if sha256(temporary) != source_hash:
            raise IOError(f"Copy verification failed: {source}")
        os.chmod(temporary, 0o644)
        try:
            os.link(temporary, target)
        except FileExistsError:
            if sha256(target) != source_hash:
                raise FileExistsError(f"Refusing to replace different published output: {target}")
        except OSError as exc:
            # Match the extension's no-clobber fallback for filesystems without hard links.
            import errno
            unsupported = {errno.EPERM, errno.EXDEV, errno.ENOSYS}
            unsupported.update(getattr(errno, name) for name in ("EOPNOTSUPP", "ENOTSUP") if hasattr(errno, name))
            if exc.errno not in unsupported:
                raise
            fd = os.open(target, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            os.close(fd)
            try:
                os.rename(temporary, target)
            except BaseException:
                target.unlink(missing_ok=True)
                raise
    finally:
        temporary.unlink(missing_ok=True)


def publish(pair: str, source_root: Path) -> dict:
    family_path = ROOT / "brain6/config/placo_family_v3/family_lock.json"
    family = json.loads(family_path.read_text())
    if pair not in family["pairs"]:
        raise ValueError(f"Pair is not in the frozen PLACO v3 family: {pair}")
    spec = family["pairs"][pair]
    pair_root = source_root / pair
    chunk_dirs = sorted(pair_root.glob("chunk_[0-9]*"))
    expected_chunks = int(spec["n_chunks"])
    if len(chunk_dirs) != expected_chunks:
        raise ValueError(f"Expected {expected_chunks} chunk artifacts, found {len(chunk_dirs)}")
    for index, chunk in enumerate(chunk_dirs):
        if chunk.name != f"chunk_{index:06d}":
            raise ValueError(f"Chunk inventory has a gap or unexpected name at {chunk}")
        verify_artifact(chunk)

    collated = pair_root / "collated"
    receipt = verify_artifact(collated)
    if receipt.get("stage") != "collate_placo" or receipt.get("scientific_status") != "PASS":
        raise ValueError("Collated PLACO artifact did not pass its frozen pair-level QC")
    status = json.loads((collated / "status.json").read_text())
    if status.get("status") != "PASS":
        raise ValueError(f"PLACO pair QC failed: {status.get('status')}")
    expected_rows = int(spec["expected_rows"])
    if int(receipt["parameters"]["expected_rows"]) != expected_rows:
        raise ValueError("Collation denominator differs from the frozen pair plan")
    threshold = float(family["headline_threshold"])
    observed = validate_variant_table(collated / "results.tsv.gz", expected_rows, threshold)
    for key, value in observed.items():
        if int(status[key if key != "tested_rows" else "tested"]) != value:
            raise ValueError(f"Collation status disagrees with results for {key}")
    if float(status["failure_rate"]) > float(family["max_pair_numerical_failure_rate"]):
        raise ValueError("PLACO numerical failure rate exceeds the frozen ceiling")

    out = ROOT / "brain6/results/placo" / pair
    output = out / "variants.tsv.gz"
    atomic_copy(collated / "results.tsv.gz", output)
    atomic_copy(collated / "receipt.json", out / "source_artifact_receipt.json")
    atomic_copy(collated / "status.json", out / "status.json")
    row = {
        "pair_id": pair, "sleep_trait": pair.split("__", 1)[0],
        "brain_disorder": pair.split("__", 1)[1],
        "stage_status": "COMPLETE_QC_PASS_LOCUS_DEDUP_PENDING",
        "total_rows": observed["total_rows"], "tested_rows": observed["tested_rows"],
        "excluded_extreme_z": observed["excluded_extreme_z"],
        "numerical_failures": observed["numerical_failures"],
        "failure_rate": status["failure_rate"], "across_track_p_threshold": threshold,
        "headline_variants": observed["headline_variants"],
        "independent_loci_status": status["independent_loci"], "output_path": str(output.relative_to(ROOT)),
        "output_sha256": sha256(output), "source_artifact_fingerprint": receipt["fingerprint"],
        "source_receipt_sha256": sha256(collated / "receipt.json"),
        "family_lock_sha256": sha256(family_path),
    }
    master = ROOT / "brain6/results/placo/placo_master.tsv"
    master.parent.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(master.with_suffix(".lock")):
        existing = {}
        if master.exists():
            with master.open(newline="", encoding="utf-8") as stream:
                for old in csv.DictReader(stream, delimiter="\t"):
                    existing[old["pair_id"]] = old
        if pair in existing and existing[pair]["output_sha256"] != row["output_sha256"]:
            raise FileExistsError(f"Master already contains a different result for {pair}")
        existing[pair] = {key: str(row[key]) for key in MASTER_FIELDS}
        fd, name = tempfile.mkstemp(prefix=".placo_master.", dir=master.parent)
        os.close(fd)
        temp_master = Path(name)
        try:
            with temp_master.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=MASTER_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(existing[key] for key in sorted(existing))
            os.chmod(temp_master, 0o644)
            os.replace(temp_master, master)
        finally:
            temp_master.unlink(missing_ok=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", required=True, help="Frozen pair ID, such as longsleep__scz")
    parser.add_argument("--source-root", type=Path,
                        default=Path("/Volumes/Extreme SSD/brain6-work/placo-family-v3-results"))
    args = parser.parse_args()
    print(json.dumps(publish(args.pair, args.source_root), indent=2))


if __name__ == "__main__":
    main()
