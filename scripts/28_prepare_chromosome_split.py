#!/usr/bin/env python3
"""Split locked munged GWAS and LDSC M totals into odd/even chromosomes."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
from pathlib import Path


def read_panel(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    traits = [row.get("trait_id", "") for row in rows]
    if len(traits) != 45 or len(set(traits)) != 45:
        raise SystemExit("ERROR: expected the locked 45-trait panel")
    return traits


def build_parity_map(reference_dir: Path) -> dict[str, str]:
    parity: dict[str, str] = {}
    for chromosome in range(1, 23):
        path = reference_dir / f"{chromosome}.l2.ldscore.gz"
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                value = "odd" if chromosome % 2 else "even"
                previous = parity.setdefault(row["SNP"], value)
                if previous != value:
                    raise SystemExit(f"ERROR: SNP maps to both chromosome parities: {row['SNP']}")
    if len(parity) < 1_000_000:
        raise SystemExit(f"ERROR: unexpectedly small HapMap3 parity map: {len(parity)}")
    return parity


def split_trait(source: Path, odd: Path, even: Path, parity: dict[str, str]) -> dict[str, int]:
    counts = {"source": 0, "odd": 0, "even": 0, "unmapped": 0}
    odd.parent.mkdir(parents=True, exist_ok=True)
    even.parent.mkdir(parents=True, exist_ok=True)
    odd_tmp = odd.with_suffix(odd.suffix + ".tmp")
    even_tmp = even.with_suffix(even.suffix + ".tmp")
    with gzip.open(source, "rt", encoding="utf-8", newline="") as src, \
            gzip.open(odd_tmp, "wt", encoding="utf-8", newline="") as odd_out, \
            gzip.open(even_tmp, "wt", encoding="utf-8", newline="") as even_out:
        header = src.readline()
        if not header.startswith("SNP\t"):
            raise SystemExit(f"ERROR: unexpected munged header: {source}")
        odd_out.write(header)
        even_out.write(header)
        for line in src:
            counts["source"] += 1
            snp = line.split("\t", 1)[0]
            target = parity.get(snp)
            if target == "odd":
                odd_out.write(line)
                counts["odd"] += 1
            elif target == "even":
                even_out.write(line)
                counts["even"] += 1
            else:
                counts["unmapped"] += 1
    mapped = counts["odd"] + counts["even"]
    if counts["source"] == 0 or mapped / counts["source"] < 0.90:
        odd_tmp.unlink(missing_ok=True)
        even_tmp.unlink(missing_ok=True)
        raise SystemExit(f"ERROR: empty or <90% reference-mapped rows in {source}: {counts}")
    os.replace(odd_tmp, odd)
    os.replace(even_tmp, even)
    return counts


def prepare_reference(source: Path, destination: Path, keep_parity: str) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for chromosome in range(1, 23):
        keep = ("odd" if chromosome % 2 else "even") == keep_parity
        ld_source = source / f"{chromosome}.l2.ldscore.gz"
        m_source = source / f"{chromosome}.l2.M_5_50"
        if not ld_source.is_file() or not m_source.is_file():
            raise SystemExit(f"ERROR: incomplete LDSC reference at chromosome {chromosome}")
        ld_target = destination / ld_source.name
        ld_target.unlink(missing_ok=True)
        ld_target.symlink_to(ld_source.resolve())
        m_target = destination / m_source.name
        if keep:
            m_target.unlink(missing_ok=True)
            m_target.symlink_to(m_source.resolve())
        else:
            values = m_source.read_text(encoding="utf-8").split()
            m_target.unlink(missing_ok=True)
            m_target.write_text("\t".join("0" for _ in values) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", default="config/analysis_panel.tsv")
    parser.add_argument("--munged-dir", default="data/munged")
    parser.add_argument("--out-dir", default="data/munged_chromosome_split")
    parser.add_argument("--ld-dir", default="ref/eur_w_ld_chr")
    parser.add_argument("--out-ld-prefix", default="ref/eur_w_ld_chr")
    parser.add_argument("--provenance", default="results/tables/chromosome_split_provenance.json")
    args = parser.parse_args()

    traits = read_panel(Path(args.panel))
    source_dir = Path(args.munged_dir)
    sources = [source_dir / f"{trait}.sumstats.gz" for trait in traits]
    missing = [str(path) for path in sources if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise SystemExit(f"ERROR: missing munged inputs: {missing}")
    free = shutil.disk_usage(Path(args.out_dir).parent).free
    source_bytes = sum(path.stat().st_size for path in sources)
    if free < source_bytes * 2:
        raise SystemExit(
            f"ERROR: chromosome split requires at least {source_bytes * 2} free bytes; found {free}"
        )

    parity = build_parity_map(Path(args.ld_dir))
    output_root = Path(args.out_dir)
    ledger = {"panel_traits": len(traits), "ld_reference_snps": len(parity), "traits": {}}
    for index, (trait, source) in enumerate(zip(traits, sources), start=1):
        counts = split_trait(
            source,
            output_root / "odd" / source.name,
            output_root / "even" / source.name,
            parity,
        )
        ledger["traits"][trait] = counts
        print(f"[{index}/45] {trait}: odd={counts['odd']} even={counts['even']}")

    prepare_reference(Path(args.ld_dir), Path(f"{args.out_ld_prefix}_odd"), "odd")
    prepare_reference(Path(args.ld_dir), Path(f"{args.out_ld_prefix}_even"), "even")
    provenance = Path(args.provenance)
    provenance.parent.mkdir(parents=True, exist_ok=True)
    provenance.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Published odd/even chromosome split for all {len(traits)} locked traits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
