#!/usr/bin/env python3
"""Build a rank-based descriptive sensitivity for the locked disorder profiles."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import yaml
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "brain6/config/brain6_locked_family.yaml"
DEFAULT_MAP = ROOT / "brain6/results/global/brain6_72_locked.tsv"
DEFAULT_OUTPUT = ROOT / "brain6/results/global/disorder_profile_similarity_spearman.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_locked_profiles(config_path: Path, map_path: Path) -> tuple[list[str], list[str], dict[str, dict[str, float]]]:
    checksum_path = config_path.with_suffix(config_path.suffix + ".sha256")
    if not checksum_path.is_file():
        raise ValueError(f"locked family checksum is missing: {checksum_path}")
    expected_hash = checksum_path.read_text(encoding="utf-8").split()[0]
    if sha256(config_path) != expected_hash:
        raise ValueError("locked family configuration checksum mismatch")
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("family_name") != "original_locked_atlas_396" or config.get("recompute_brain6_fdr") is not False:
        raise ValueError("family config does not preserve the original locked atlas family")
    sleeps = config.get("sleep_traits")
    disorders = config.get("brain_disorders")
    if not isinstance(sleeps, list) or not isinstance(disorders, list) or len(set(sleeps)) != len(sleeps) or len(set(disorders)) != len(disorders):
        raise ValueError("family config must provide unique sleep traits and brain disorders")

    profiles = {disorder: {} for disorder in disorders}
    with map_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != len(sleeps) * len(disorders):
        raise ValueError("global map does not have the configured sleep × disorder dimensions")
    for row in rows:
        sleep, disorder = row["sleep_trait"], row["brain_disorder"]
        if sleep not in sleeps or disorder not in disorders or sleep in profiles[disorder]:
            raise ValueError(f"unexpected or duplicate profile cell: {sleep} × {disorder}")
        try:
            value = float(row["rg"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"non-numeric rg for {sleep} × {disorder}") from exc
        if not (-1 <= value <= 1):
            raise ValueError(f"rg outside [-1, 1] for {sleep} × {disorder}")
        profiles[disorder][sleep] = value
    if any(list(profile) != sleeps for profile in profiles.values()):
        # CSV row order is not significant; validate coverage and restore manifest order.
        if any(set(profile) != set(sleeps) for profile in profiles.values()):
            raise ValueError("one or more profiles are missing configured sleep traits")
    return sleeps, disorders, profiles


def render_table(config_path: Path, map_path: Path) -> str:
    sleeps, disorders, profiles = read_locked_profiles(config_path, map_path)
    lines = ["brain_disorder\t" + "\t".join(disorders)]
    for left in disorders:
        values = []
        for right in disorders:
            coefficient = spearmanr(
                [profiles[left][sleep] for sleep in sleeps],
                [profiles[right][sleep] for sleep in sleeps],
            ).statistic
            values.append(f"{coefficient:.12g}")
        lines.append(left + "\t" + "\t".join(values))
    return "\n".join(lines) + "\n"


def write_exclusive(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
        stream.write(content)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--map", dest="map_path", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    table = render_table(args.config, args.map_path)
    write_exclusive(args.output, table)
    provenance = {
        "status": "PASS",
        "method": "Spearman rank correlation across the locked sleep-trait rg profile",
        "interpretation": "descriptive sensitivity only; no inferential p-values or new discovery family",
        "config_sha256": sha256(args.config),
        "global_map_sha256": sha256(args.map_path),
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": sha256(args.output),
    }
    provenance_path = args.output.with_suffix(".provenance.json")
    write_exclusive(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), **provenance}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
