#!/usr/bin/env python3
"""Extract one worksheet from an XLSX archive using only the Python stdlib.

This is intentionally read-only with respect to the source workbook. It exists
to make prior-study phenotype lists reviewable as plain, diffable TSV files
without introducing an Excel-library dependency into the analysis runtime.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path, PurePosixPath
import re
import xml.etree.ElementTree as ET
import zipfile


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL_REF = re.compile(r"^([A-Z]+)([0-9]+)$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def column_number(letters: str) -> int:
    value = 0
    for letter in letters:
        value = value * 26 + ord(letter) - ord("A") + 1
    return value


def archive_path(target: str) -> str:
    target = target.lstrip("/")
    if target.startswith("xl/"):
        return target
    return str(PurePosixPath("xl") / target)


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    result: list[str] = []
    for item in root.findall(f"{{{MAIN_NS}}}si"):
        result.append("".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t")))
    return result


def sheet_target(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        node.attrib["Id"]: node.attrib["Target"]
        for node in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")
    }
    for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet"):
        if sheet.attrib.get("name") == sheet_name:
            rel_id = sheet.attrib.get(f"{{{REL_NS}}}id")
            if not rel_id or rel_id not in targets:
                fail(f"worksheet relationship missing for {sheet_name}")
            return archive_path(targets[rel_id])
    available = [
        sheet.attrib.get("name", "")
        for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet")
    ]
    fail(f"worksheet {sheet_name!r} not found; available: {available}")
    raise AssertionError("unreachable")


def cell_value(cell: ET.Element, strings: list[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    value = cell.find(f"{{{MAIN_NS}}}v")
    raw = "" if value is None or value.text is None else value.text
    if cell_type == "s" and raw:
        return strings[int(raw)]
    if cell_type == "inlineStr":
        inline = cell.find(f"{{{MAIN_NS}}}is")
        if inline is None:
            return ""
        return "".join(node.text or "" for node in inline.iter(f"{{{MAIN_NS}}}t"))
    if cell_type == "b":
        return "true" if raw == "1" else "false"
    return raw


def extract(workbook: Path, sheet_name: str) -> list[list[str]]:
    with zipfile.ZipFile(workbook) as archive:
        strings = shared_strings(archive)
        target = sheet_target(archive, sheet_name)
        root = ET.fromstring(archive.read(target))

    sparse: dict[int, dict[int, str]] = {}
    max_column = 0
    for cell in root.findall(f".//{{{MAIN_NS}}}c"):
        reference = cell.attrib.get("r", "")
        match = CELL_REF.match(reference)
        if not match:
            continue
        column = column_number(match.group(1))
        row = int(match.group(2))
        sparse.setdefault(row, {})[column] = cell_value(cell, strings)
        max_column = max(max_column, column)
    if not sparse:
        return []

    rows: list[list[str]] = []
    for row_number in range(1, max(sparse) + 1):
        values = sparse.get(row_number, {})
        row = [values.get(column, "") for column in range(1, max_column + 1)]
        while row and row[-1] == "":
            row.pop()
        rows.append(row)
    while rows and not rows[-1]:
        rows.pop()
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if not args.workbook.is_file():
        fail(f"workbook not found: {args.workbook}")
    rows = extract(args.workbook, args.sheet)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    print(f"WROTE {args.out} sheet={args.sheet} rows={len(rows)}")


if __name__ == "__main__":
    main()

