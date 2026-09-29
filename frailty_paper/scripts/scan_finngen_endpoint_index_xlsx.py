#!/usr/bin/env python3
"""Scan all worksheet cell text in an XLSX, ignoring unreliable sheet dimensions."""

from __future__ import annotations

import argparse
import hashlib
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkg": "http://schemas.openxmlformats.org/package/2006/relationships",
}
TERMS = re.compile(r"frailty|hfrs|hospital\s+frailty", re.IGNORECASE)


def scan(path: Path) -> tuple[str, int, int, list[tuple[int, str]]]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in relationships.findall("pkg:Relationship", NS)
        }
        sheets = workbook.findall("main:sheets/main:sheet", NS)
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("main:si", NS):
                shared.append("".join(t.text or "" for t in item.findall(".//main:t", NS)))

        row_total = 0
        nonempty_rows = 0
        hits: list[tuple[int, str]] = []
        for sheet in sheets:
            rel_id = sheet.attrib[f"{{{NS['rel']}}}id"]
            target = targets[rel_id].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            root = ET.fromstring(archive.read(target))
            for row in root.findall(".//main:sheetData/main:row", NS):
                row_total += 1
                values: list[str] = []
                for cell in row.findall("main:c", NS):
                    kind = cell.attrib.get("t")
                    if kind == "s":
                        value = cell.findtext("main:v", default="", namespaces=NS)
                        text = shared[int(value)] if value else ""
                    elif kind == "inlineStr":
                        text = "".join(t.text or "" for t in cell.findall(".//main:t", NS))
                    else:
                        text = cell.findtext("main:v", default="", namespaces=NS)
                    if text:
                        values.append(text)
                if values:
                    nonempty_rows += 1
                    joined = " | ".join(values)
                    if TERMS.search(joined):
                        hits.append((int(row.attrib.get("r", "0")), joined))
    return digest, row_total, nonempty_rows, hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xlsx", type=Path)
    args = parser.parse_args()
    digest, rows, nonempty, hits = scan(args.xlsx)
    print(f"file={args.xlsx}")
    print(f"bytes={args.xlsx.stat().st_size}")
    print(f"sha256={digest}")
    print(f"worksheet_rows={rows}")
    print(f"nonempty_rows={nonempty}")
    print(f"matching_rows={len(hits)}")
    for row_number, content in hits:
        print(f"row_{row_number}={content}")


if __name__ == "__main__":
    main()
