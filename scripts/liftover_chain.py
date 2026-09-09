"""Fail-closed point liftover using a checksum-pinned UCSC chain file."""
from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
import gzip
import hashlib
from pathlib import Path
import re


AUTOSOME = re.compile(r"(?:chr)?([1-9]|1[0-9]|2[0-2])$", re.IGNORECASE)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class LiftoverError(ValueError):
    pass


def file_hash(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def autosome(value) -> int | None:
    match = AUTOSOME.fullmatch(str(value).strip())
    return int(match.group(1)) if match else None


def reverse_complement(allele: str) -> str:
    return str(allele).translate(COMPLEMENT)[::-1]


class ChainIndex:
    def __init__(self, intervals_by_chromosome):
        self.intervals = {}
        for chromosome, intervals in intervals_by_chromosome.items():
            intervals.sort(key=lambda interval: (interval[0], interval[1], interval[6]))
            starts = [interval[0] for interval in intervals]
            prefix_max_end = []
            maximum = 0
            for interval in intervals:
                maximum = max(maximum, interval[1])
                prefix_max_end.append(maximum)
            self.intervals[chromosome] = (starts, prefix_max_end, intervals)

    def map_point(self, chromosome, position):
        """Map a 1-based point; reject absent, non-autosomal, or ambiguous hits."""
        source_chromosome = autosome(chromosome)
        try:
            source_position = int(position)
        except (TypeError, ValueError, OverflowError):
            return "invalid_source_coordinate", None
        if source_chromosome is None or source_position <= 0:
            return "invalid_source_coordinate", None
        indexed = self.intervals.get(source_chromosome)
        if indexed is None:
            return "unmapped", None
        coordinate = source_position - 1
        starts, prefix_max_end, intervals = indexed
        cursor = bisect_right(starts, coordinate) - 1
        hits = set()
        while cursor >= 0 and prefix_max_end[cursor] > coordinate:
            (
                source_start,
                source_end,
                target_chromosome,
                target_start,
                target_size,
                target_strand,
                _chain_id,
            ) = intervals[cursor]
            if source_start <= coordinate < source_end:
                offset = coordinate - source_start
                if target_strand == "+":
                    target_coordinate = target_start + offset
                else:
                    target_coordinate = target_size - 1 - (target_start + offset)
                hits.add((target_chromosome, target_coordinate + 1, target_strand))
            cursor -= 1
        if not hits:
            return "unmapped", None
        if len(hits) != 1:
            return "ambiguous", None
        target = hits.pop()
        if target[0] is None:
            return "non_autosomal_target", None
        return "mapped", target


def load_chain(
    chain_path: str | Path,
    expected_sha256: str | None = None,
    expected_bytes: int | None = None,
):
    path = Path(chain_path)
    if not path.is_file():
        raise LiftoverError(f"liftover chain not found: {path}")
    actual_bytes = path.stat().st_size
    if expected_bytes is not None and actual_bytes != expected_bytes:
        raise LiftoverError(
            f"liftover chain bytes differ from the registered chain: "
            f"expected {expected_bytes}, got {actual_bytes}"
        )
    actual_sha256 = file_hash(path)
    if expected_sha256 is not None and actual_sha256 != expected_sha256:
        raise LiftoverError(
            f"liftover chain SHA-256 differs from the registered chain: "
            f"expected {expected_sha256}, got {actual_sha256}"
        )

    intervals_by_chromosome = defaultdict(list)
    chain_count = 0
    block_count = 0
    current = None
    source_cursor = target_cursor = None
    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.strip().split()
                if not fields:
                    current = None
                    source_cursor = target_cursor = None
                    continue
                if fields[0] == "chain":
                    if len(fields) != 13:
                        raise LiftoverError(
                            f"invalid chain header at line {line_number}: {fields}"
                        )
                    try:
                        (
                            _chain,
                            _score,
                            source_name,
                            _source_size,
                            source_strand,
                            source_start,
                            _source_end,
                            target_name,
                            target_size,
                            target_strand,
                            target_start,
                            _target_end,
                            chain_id,
                        ) = fields
                        source_cursor = int(source_start)
                        target_cursor = int(target_start)
                        target_size = int(target_size)
                        chain_id = int(chain_id)
                    except ValueError as error:
                        raise LiftoverError(
                            f"non-numeric chain header field at line {line_number}"
                        ) from error
                    if source_strand != "+" or target_strand not in {"+", "-"}:
                        raise LiftoverError(
                            f"unsupported chain orientation at line {line_number}"
                        )
                    current = (
                        autosome(source_name), autosome(target_name), target_size,
                        target_strand, chain_id,
                    )
                    chain_count += 1
                    continue
                if current is None or len(fields) not in {1, 3}:
                    raise LiftoverError(f"invalid chain block at line {line_number}")
                try:
                    size = int(fields[0])
                    source_gap = int(fields[1]) if len(fields) == 3 else 0
                    target_gap = int(fields[2]) if len(fields) == 3 else 0
                except ValueError as error:
                    raise LiftoverError(
                        f"non-numeric chain block at line {line_number}"
                    ) from error
                if size <= 0 or source_gap < 0 or target_gap < 0:
                    raise LiftoverError(f"invalid chain block size/gap at line {line_number}")
                source_chromosome, target_chromosome, target_size, target_strand, chain_id = current
                if source_chromosome is not None:
                    intervals_by_chromosome[source_chromosome].append(
                        (
                            source_cursor,
                            source_cursor + size,
                            target_chromosome,
                            target_cursor,
                            target_size,
                            target_strand,
                            chain_id,
                        )
                    )
                source_cursor += size + source_gap
                target_cursor += size + target_gap
                block_count += 1
    except (OSError, UnicodeDecodeError) as error:
        raise LiftoverError(f"could not read liftover chain: {error}") from error
    if not chain_count or not block_count or not intervals_by_chromosome:
        raise LiftoverError("liftover chain contains no autosomal mapping blocks")
    return ChainIndex(intervals_by_chromosome), {
        "chain_path": str(path),
        "chain_bytes": actual_bytes,
        "chain_sha256": actual_sha256,
        "chain_count": chain_count,
        "block_count": block_count,
    }
