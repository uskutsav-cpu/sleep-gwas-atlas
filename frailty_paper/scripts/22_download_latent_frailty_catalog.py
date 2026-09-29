#!/usr/bin/env python3
"""Acquire accession-indexed latent-frailty Catalog files with official MD5 checks."""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
import yaml

FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch_checksum(session: requests.Session, base: str, filename: str) -> str:
    response = session.get(f"{base}/md5sum.txt", timeout=60)
    response.raise_for_status()
    match = re.search(rf"^([0-9a-fA-F]{{32}})\s+\*?{re.escape(filename)}\s*$", response.text, re.M)
    if not match:
        raise RuntimeError(f"official md5sum.txt has no entry for {filename} at {base}")
    return match.group(1).lower()


def acquire(session: requests.Session, base: str, filename: str, dest: Path, expected_md5: str) -> None:
    if dest.is_file() and md5(dest) == expected_md5:
        print(f"verified existing {dest.name}", flush=True)
        return
    if dest.exists():
        quarantine = dest.with_name(dest.name + ".unverified")
        if quarantine.exists():
            quarantine = dest.with_name(dest.name + f".unverified.{dest.stat().st_mtime_ns}")
        dest.replace(quarantine)
        print(f"retained mismatched prior file at {quarantine}", flush=True)
    part = dest.with_name(dest.name + ".part")
    url = f"{base}/{filename}"
    subprocess.run([
        "curl", "--fail", "--location", "--silent", "--show-error", "--retry", "5",
        "--retry-delay", "3", "--retry-all-errors", "--continue-at", "-",
        "--header", "Accept-Encoding: identity", "--output", str(part), url,
    ], check=True)
    actual = md5(part)
    if actual != expected_md5:
        quarantine = part.with_name(part.name + ".unverified")
        if quarantine.exists():
            quarantine = part.with_name(part.name + f".unverified.{part.stat().st_mtime_ns}")
        part.replace(quarantine)
        raise RuntimeError(f"MD5 mismatch for {filename}: expected {expected_md5}, got {actual}; retained {quarantine}")
    os.replace(part, dest)
    print(f"downloaded and MD5-verified {dest}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--accession", action="append", help="Download only this accession; repeatable")
    parser.add_argument("--jobs", type=int, default=3, help="Concurrent accession downloads (default: 3)")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    root = args.repo.resolve()
    index_path = root / "frailty_paper/manifests/latent_frailty_accession_index.tsv"
    with index_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    for field in ("genome_assembly", "sample_ancestry", "sample_size", "license"):
        if field not in fields:
            fields.append(field)
    selected = set(args.accession or [])
    if selected - {r["accession"] for r in rows}:
        raise SystemExit(f"unknown accession(s): {', '.join(sorted(selected - {r['accession'] for r in rows}))}")
    out_root = root / "frailty_paper/data/gwas/latent_frailty_catalog"
    out_root.mkdir(parents=True, exist_ok=True)
    selected_rows = [row for row in rows if not selected or row["accession"] in selected]

    def process(row: dict[str, str]) -> dict[str, str]:
        accession = row["accession"]
        base = f"{FTP_ROOT}/{accession}"
        tsv_name = f"{accession}.tsv"
        meta_name = f"{accession}.tsv-meta.yaml"
        folder = out_root / accession
        folder.mkdir(exist_ok=True)
        with requests.Session() as session:
            session.headers.update({"User-Agent": "sleep-gwas-atlas-frailty-paper/1.0 (reproducible public-data acquisition)"})
            for filename in (tsv_name, meta_name):
                expected = fetch_checksum(session, base, filename)
                acquire(session, base, filename, folder / filename, expected)
            # Retain the small official checksum listing for human audit.
            checksum_file = folder / "md5sum.txt"
            response = session.get(f"{base}/md5sum.txt", timeout=60)
            response.raise_for_status()
            checksum_file.write_text(response.text, encoding="ascii")
            metadata = yaml.safe_load((folder / meta_name).read_text(encoding="utf-8")) or {}
            trait = metadata.get("trait_description") or []
            if isinstance(trait, str):
                trait = [trait]
            trait = [str(value).strip() for value in trait if str(value).strip()]
            if trait:
                row["phenotype_label"] = "; ".join(trait)
                row["label_source"] = "Official GWAS Catalog FTP metadata YAML"
                row["label_status"] = "VERIFIED"
            row["genome_assembly"] = metadata.get("genome_assembly", "UNKNOWN")
            samples = metadata.get("samples") or []
            ancestries, sizes = [], []
            for sample in samples:
                if not isinstance(sample, dict):
                    continue
                ancestry = sample.get("sample_ancestry_category") or []
                if isinstance(ancestry, str):
                    ancestry = [ancestry]
                ancestries.extend(str(value) for value in ancestry if value)
                if sample.get("sample_size") is not None:
                    sizes.append(str(sample["sample_size"]))
            row["sample_ancestry"] = "; ".join(dict.fromkeys(ancestries)) or "UNKNOWN"
            row["sample_size"] = "; ".join(sizes) or "UNKNOWN"
            license_value = metadata.get("license") or "UNKNOWN"
            if isinstance(license_value, dict):
                license_value = license_value.get("license_uri") or license_value.get("name") or str(license_value)
            row["license"] = str(license_value)
            row["local_acquisition_status"] = "ACQUIRED_MD5_VERIFIED"
            row["local_qc_status"] = "NOT_RUN"
            row["notes"] = "Phenotype label, genome assembly, ancestry and sample size transcribed from the accession's official GWAS Catalog TSV metadata YAML; TSV and YAML MD5s verified against the official md5sum.txt."
        return row

    errors = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(process, row): row["accession"] for row in selected_rows}
        for future in as_completed(futures):
            try:
                result = future.result()
                for index, row in enumerate(rows):
                    if row["accession"] == result["accession"]:
                        rows[index] = result
                        break
            except Exception as exc:
                errors.append(f"{futures[future]}: {exc}")
    if errors:
        raise SystemExit("latent-frailty acquisition errors:\n" + "\n".join(errors))
    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"LATENT_FRAILTY_ACQUISITION_DONE records={sum(1 for r in rows if r.get('local_acquisition_status') == 'ACQUIRED_MD5_VERIFIED')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
