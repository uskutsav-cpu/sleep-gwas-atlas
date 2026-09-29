#!/usr/bin/env python3
"""Normalize licensed database exports without modifying their source files."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import unicodedata
from pathlib import Path

FIELDS = [
    "record_id", "record_type", "pmid", "doi", "normalized_title", "title",
    "abstract", "year", "journal", "query_id", "source_file", "source_database",
    "source_record_id", "authors", "import_status",
]
MANIFEST_FIELDS = ["source_file", "source_database", "format", "bytes", "sha256", "records_imported", "records_missing_title"]
DB_NAMES = {
    "embase": "Embase", "scopus": "Scopus", "webofscience": "Web of Science",
    "wos": "Web of Science", "web_of_science": "Web of Science",
    "psycinfo": "PsycINFO",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_title(value: str) -> str:
    value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode().lower()
    return " ".join(re.findall(r"[a-z0-9]+", value))


def database_from_filename(path: Path) -> str:
    stem = re.sub(r"[^a-z0-9]+", "_", path.stem.lower()).strip("_")
    for key in sorted(DB_NAMES, key=len, reverse=True):
        if stem == key or stem.startswith(key + "_") or stem.endswith("_" + key):
            return DB_NAMES[key]
    raise ValueError(f"filename must identify a supported database (Embase/Scopus/Web of Science/PsycINFO): {path.name}")


def first(values: dict[str, list[str]], *keys: str) -> str:
    for key in keys:
        for value in values.get(key, []):
            value = " ".join(value.split())
            if value:
                return value
    return ""


def normalized_record(database: str, path: Path, relative: str, source_id: str,
                      ordinal: int, title: str, abstract: str, year: str,
                      journal: str, doi: str, pmid: str, authors: str) -> dict[str, str]:
    title = " ".join(title.split())
    status = "READY" if title else "MISSING_TITLE"
    stable_id = source_id or sha256_bytes(f"{relative}\n{ordinal}\n{title}\n{doi}\n{pmid}".encode())[:20]
    return {
        "record_id": f"manual:{database}:{stable_id}:{ordinal}",
        "record_type": "manual_database_record", "pmid": pmid.strip(),
        "doi": doi.strip(), "normalized_title": normalize_title(title),
        "title": title, "abstract": " ".join(abstract.split()),
        "year": (re.search(r"(?:19|20)\d{2}", year) or [""])[0],
        "journal": " ".join(journal.split()), "query_id": "manual_" + re.sub(r"[^a-z0-9]+", "_", database.lower()).strip("_"),
        "source_file": relative, "source_database": database,
        "source_record_id": source_id or stable_id, "authors": " ".join(authors.split()),
        "import_status": status,
    }


def parse_ris(payload: str, database: str, path: Path, relative: str) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    fields: dict[str, list[str]] = {}
    last_tag = ""

    def flush() -> None:
        nonlocal fields, last_tag
        if not fields:
            return
        source_id = first(fields, "AN", "EID", "UT", "ID", "UR", "M1")
        records.append(normalized_record(
            database, path, relative, source_id, len(records) + 1,
            first(fields, "TI", "T1", "CT", "TT"), first(fields, "AB", "N2"),
            first(fields, "PY", "Y1", "DA"), first(fields, "JO", "JF", "T2", "JA", "PB"),
            first(fields, "DO", "DI"), first(fields, "PMID", "PM"),
            "; ".join(fields.get("AU", []) + fields.get("A1", [])),
        ))
        fields = {}
        last_tag = ""

    for line in payload.splitlines():
        match = re.match(r"^([A-Z0-9]{2,})\s+-\s?(.*)$", line)
        if match:
            tag, value = match.group(1), match.group(2).strip()
            if tag == "TY" and fields:
                flush()
            if tag == "ER":
                flush()
                continue
            fields.setdefault(tag, []).append(value)
            last_tag = tag
        elif line[:1].isspace() and last_tag:
            fields[last_tag][-1] += " " + line.strip()
    flush()
    return records


def canonical_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def find_column(row: dict[str, str], *candidates: str) -> str:
    by_name = {canonical_header(key): value or "" for key, value in row.items() if key}
    for candidate in candidates:
        value = by_name.get(canonical_header(candidate), "")
        if value.strip():
            return value.strip()
    return ""


def parse_csv(payload: bytes, database: str, path: Path, relative: str) -> list[dict[str, str]]:
    from io import StringIO

    text = payload.decode("utf-8-sig")
    reader = csv.DictReader(StringIO(text))
    if not reader.fieldnames:
        raise ValueError("CSV has no header row")
    records: list[dict[str, str]] = []
    for ordinal, row in enumerate(reader, 1):
        title = find_column(row, "Title", "Article Title", "Document Title")
        source_id = find_column(row, "EID", "Accession Number", "UT (Unique WOS ID)", "UT", "Record ID", "Document ID", "Embase ID")
        records.append(normalized_record(
            database, path, relative, source_id, ordinal, title,
            find_column(row, "Abstract", "Description"),
            find_column(row, "Year", "Publication Year", "Pub Year", "Date"),
            find_column(row, "Source title", "Source Title", "Publication Title", "Journal"),
            find_column(row, "DOI", "Digital Object Identifier"),
            find_column(row, "PMID", "PubMed ID", "Pubmed ID"),
            find_column(row, "Authors", "Author(s)", "Author", "Creators"),
        ))
    return records


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def import_exports(directory: Path, root: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    records: list[dict[str, str]] = []
    manifest: list[dict[str, str]] = []
    files = sorted(p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in {".ris", ".csv"})
    for path in files:
        database = database_from_filename(path)
        payload = path.read_bytes()
        relative = path.relative_to(root).as_posix()
        try:
            if path.suffix.lower() == ".ris":
                parsed = parse_ris(payload.decode("utf-8-sig"), database, path, relative)
                fmt = "RIS"
            else:
                parsed = parse_csv(payload, database, path, relative)
                fmt = "CSV"
        except (UnicodeDecodeError, csv.Error, ValueError) as exc:
            raise ValueError(f"failed to import {relative}: {exc}") from exc
        records.extend(parsed)
        manifest.append({
            "source_file": relative, "source_database": database, "format": fmt,
            "bytes": str(len(payload)), "sha256": sha256_bytes(payload),
            "records_imported": str(len(parsed)),
            "records_missing_title": str(sum(not row["title"] for row in parsed)),
        })
    return records, manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=".")
    args = parser.parse_args()
    root = Path(args.repo).resolve()
    directory = root / "frailty_paper/review/manual_exports"
    records, manifest = import_exports(directory, root)
    if not manifest:
        stale = directory / "manual_normalized_records.tsv"
        if stale.is_file():
            raise SystemExit(f"No raw RIS/CSV exports found but normalized import still exists: {stale}; restore the immutable exports or reconcile the import explicitly")
        print("NO_MANUAL_EXPORTS_FOUND; untouched source files and generated import outputs unchanged")
        return 0
    write_tsv(directory / "manual_normalized_records.tsv", FIELDS, records)
    write_tsv(directory / "manual_import_manifest.tsv", MANIFEST_FIELDS, manifest)
    print(f"Imported {len(records)} records from {len(manifest)} immutable source exports; missing titles={sum(not r['title'] for r in records)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
