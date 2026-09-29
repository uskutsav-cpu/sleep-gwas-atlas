#!/usr/bin/env python3
"""Verify source-card identity and bytes for every trait in the frozen pair lock."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "extensions" / "brain6"))
from brain6.io import json_hash, sha256, write_tsv


FIELDS = ["pair_lock_sha256", "selection_sha256", "trait_id", "source_card", "source_card_sha256", "dense_path",
          "dense_bytes", "expected_sha256", "observed_sha256", "status"]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-lock", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=Path,
                        help="Reviewed trait-to-source-card JSON registry")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    lock = json.loads(args.pair_lock.read_text(encoding="utf-8"))
    registry = json.loads(args.sources.read_text(encoding="utf-8"))
    lock_digest = lock.get("lock_sha256")
    lock_body = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if lock.get("synthetic", False) or json_hash(lock_body) != lock_digest:
        raise SystemExit("Pair lock is synthetic or its checksum is invalid")
    root = Path(__file__).resolve().parents[2]
    selection = root / "brain6/config/deep_tracks_v1.tsv"
    selection_sha = sha256(selection)
    expected_selection_sha = selection.with_suffix(selection.suffix + ".sha256").read_text(
        encoding="utf-8").split()[0]
    if selection_sha != expected_selection_sha:
        raise SystemExit("Frozen deep-track selection checksum is invalid")
    traits = sorted({pair[key] for pair in lock["pairs"]
                     for key in ("sleep_trait", "disease_trait")})
    missing = set(traits) - set(registry)
    if missing:
        raise SystemExit(f"Source registry is missing locked traits: {sorted(missing)}")

    rows = []
    for trait in traits:
        card_path = Path(registry[trait]).resolve(strict=True)
        card = json.loads(card_path.read_text(encoding="utf-8"))
        if card.get("trait_id") != trait:
            raise SystemExit(f"Source-card trait mismatch for {trait}: {card.get('trait_id')!r}")
        if card.get("synthetic", False) or not card.get("access_permitted", False):
            raise SystemExit(f"Refusing synthetic or unauthorized source card for {trait}")
        dense_path = Path(card["path"]).resolve(strict=True)
        expected = str(card["sha256"]).lower()
        observed = digest(dense_path)
        status = "VERIFIED" if observed == expected else "HASH_MISMATCH"
        rows.append({"pair_lock_sha256": lock_digest, "selection_sha256": selection_sha,
                     "trait_id": trait, "source_card": str(card_path),
                     "source_card_sha256": digest(card_path),
                     "dense_path": str(dense_path), "dense_bytes": dense_path.stat().st_size,
                     "expected_sha256": expected, "observed_sha256": observed,
                     "status": status})

    write_tsv(args.out, FIELDS, rows)
    if any(row["status"] != "VERIFIED" for row in rows):
        raise SystemExit("One or more dense source files do not match their reviewed card")
    print(f"Verified {len(rows)} locked trait inputs; manifest={args.out}")


if __name__ == "__main__":
    main()
