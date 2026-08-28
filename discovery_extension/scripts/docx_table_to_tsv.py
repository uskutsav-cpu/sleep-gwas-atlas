#!/usr/bin/env python3
"""Extract one OOXML Word table to a plain TSV without modifying the source."""

from __future__ import annotations

import argparse
import csv
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--document", type=Path, required=True)
    parser.add_argument("--table", type=int, default=1, help="1-based table number")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if args.table < 1:
        raise SystemExit("ERROR: --table must be at least 1")
    with zipfile.ZipFile(args.document) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    tables = root.findall(f".//{W}tbl")
    if args.table > len(tables):
        raise SystemExit(f"ERROR: table {args.table} not found; document has {len(tables)}")

    rows: list[list[str]] = []
    for tr in tables[args.table - 1].findall(f"./{W}tr"):
        row: list[str] = []
        for tc in tr.findall(f"./{W}tc"):
            paragraphs = []
            for paragraph in tc.findall(f".//{W}p"):
                text = "".join(node.text or "" for node in paragraph.findall(f".//{W}t"))
                if text:
                    paragraphs.append(text)
            row.append(" ".join(paragraphs).strip())
        rows.append(row)

    width = max((len(row) for row in rows), default=0)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(row + [""] * (width - len(row)) for row in rows)
    print(f"WROTE {args.out} table={args.table} rows={len(rows)} columns={width}")


if __name__ == "__main__":
    main()
