#!/usr/bin/env python3
"""Build pair-specific LAVA covariance inputs from the locked atlas LDSC intercepts."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions" / "brain6"))
from brain6.io import atomic_text, json_hash, sha256, write_tsv


FIELDS = ["pair_lock_sha256", "atlas_map_sha256", "trait_metadata_sha256",
          "pair_id", "trait1", "trait2", "source_id1", "source_id2",
          "source_card_sha256_1", "source_card_sha256_2",
          "sleep_intercept", "disorder_intercept", "cross_trait_intercept",
          "cross_trait_intercept_se", "error_correlation", "matrix_path",
          "matrix_sha256", "status"]


def number(value: str, name: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise SystemExit(f"Non-finite {name}: {value}")
    return x


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-lock", type=Path,
                        default=ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json")
    parser.add_argument("--sources", type=Path,
                        default=ROOT / "extensions/brain6/work/overnight-v03/sources.reviewed.json")
    parser.add_argument("--atlas-map", type=Path,
                        default=ROOT / "brain6/results/global/brain6_72_locked.tsv")
    parser.add_argument("--traits", type=Path, default=ROOT / "results/atlas/traits.tsv")
    parser.add_argument("--out-dir", type=Path,
                        default=ROOT / "brain6/manifests/lava_overlap_v1")
    parser.add_argument("--version", default="v1", choices=["v1", "v2"])
    args = parser.parse_args()
    args.out_dir = args.out_dir.resolve()
    manifest_name = f"lava_overlap_family_{args.version}.tsv"
    if args.out_dir.exists():
        raise SystemExit("Refusing to overwrite an existing LAVA overlap family")

    lock = json.loads(args.pair_lock.read_text(encoding="utf-8"))
    lock_digest = lock.get("lock_sha256")
    if lock.get("synthetic", False) or json_hash(
            {key: value for key, value in lock.items() if key != "lock_sha256"}) != lock_digest:
        raise SystemExit("Pair lock is synthetic or its checksum is invalid")
    sources = json.loads(args.sources.read_text(encoding="utf-8"))
    traits = {row["trait_id"]: row for row in csv.DictReader(
        args.traits.open(encoding="utf-8", newline=""), delimiter="\t")}
    atlas_rows = list(csv.DictReader(args.atlas_map.open(encoding="utf-8", newline=""),
                                     delimiter="\t"))
    pair_lock_sha = lock_digest
    atlas_map_sha = sha256(args.atlas_map)
    trait_metadata_sha = sha256(args.traits)
    args.out_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=".lava-overlap-v1.", dir=args.out_dir.parent))

    rows = []
    try:
        for pair in lock["pairs"]:
            sleep, disorder = pair["sleep_trait"], pair["disease_trait"]
            candidates = [row for row in atlas_rows
                          if row["sleep_trait"] == sleep and row["external_trait"] == disorder]
            if len(candidates) != 1:
                raise SystemExit(f"Expected one locked atlas row for {sleep} × {disorder}")
            atlas = candidates[0]
            expected_sources = [(sleep, atlas["source_GWAS_ID_sleep"], "sleep"),
                                (disorder, atlas["source_GWAS_ID_disorder"], "disorder")]
            cards = []
            for trait_id, expected_id, side in expected_sources:
                if trait_id not in sources or trait_id not in traits:
                    raise SystemExit(f"Missing source/trait metadata for {trait_id}")
                card_path = Path(sources[trait_id]).resolve(strict=True)
                card = json.loads(card_path.read_text(encoding="utf-8"))
                metadata = traits[trait_id]
                if (card.get("study_id") != expected_id or metadata.get("source_id") != expected_id
                        or card.get("sha256") != metadata.get("harmonized_sha256")):
                    raise SystemExit(f"Source identity/hash mismatch for {side} trait {trait_id}")
                cards.append((card_path, card))

            h1 = number(atlas["sleep_h2_intercept"], "sleep LDSC intercept")
            h2 = number(atlas["external_h2_intercept"], "disorder LDSC intercept")
            cross = number(atlas["LDSC_cross_trait_intercept"], "cross-trait LDSC intercept")
            cross_se = number(atlas["LDSC_cross_trait_intercept_se"], "cross-trait intercept SE")
            if h1 <= 0 or h2 <= 0 or cross_se <= 0:
                raise SystemExit(f"Invalid LDSC intercept inputs for {pair['pair_id']}")
            if abs(cross) >= math.sqrt(h1 * h2):
                raise SystemExit(f"Non-positive-definite two-trait intercept matrix for {pair['pair_id']}")
            rho = cross / math.sqrt(h1 * h2)

            staged_matrix = temp_dir / f"{pair['pair_id']}.txt"
            with atomic_text(staged_matrix) as stream:
                stream.write(f" {sleep} {disorder}\n")
                stream.write(f"{sleep} {h1:.8f} {cross:.8f}\n")
                stream.write(f"{disorder} {cross:.8f} {h2:.8f}\n")
            rows.append({"pair_lock_sha256": pair_lock_sha,
                         "atlas_map_sha256": atlas_map_sha,
                         "trait_metadata_sha256": trait_metadata_sha,
                         "pair_id": pair["pair_id"], "trait1": sleep, "trait2": disorder,
                         "source_id1": cards[0][1]["study_id"],
                         "source_id2": cards[1][1]["study_id"],
                         "source_card_sha256_1": cards[0][1]["sha256"],
                         "source_card_sha256_2": cards[1][1]["sha256"],
                         "sleep_intercept": h1, "disorder_intercept": h2,
                         "cross_trait_intercept": cross,
                         "cross_trait_intercept_se": cross_se,
                         "error_correlation": rho,
                         "matrix_path": str((args.out_dir / staged_matrix.name).relative_to(ROOT)),
                         "matrix_sha256": sha256(staged_matrix),
                         "status": "PASS_LDSC_INTERCEPT_DERIVED"})
        if len(rows) != 5:
            raise SystemExit(f"Expected five primary pairs, found {len(rows)}")
        write_tsv(temp_dir / manifest_name, FIELDS, rows)
        os.rename(temp_dir, args.out_dir)
    except BaseException:
        # A partial directory is never published under the locked output name.
        raise
    print(f"Prepared {len(rows)} pair-specific LAVA overlap matrices.")
    print(f"Atlas map SHA256: {sha256(args.atlas_map)}")
    print(f"Trait metadata SHA256: {sha256(args.traits)}")
    print(f"Pair lock SHA256: {lock_digest}")
    print(f"Overlap manifest: {args.out_dir / manifest_name}")


if __name__ == "__main__":
    main()
