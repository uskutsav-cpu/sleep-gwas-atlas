#!/usr/bin/env python3
"""Literal-token crosswalk of the paper's HFRS ICD-10 codes to FinnGen DF13."""

from __future__ import annotations

import argparse
import csv
import hashlib
import posixpath
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkg": "http://schemas.openxmlformats.org/package/2006/relationships",
}
HFRS_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs43587-025-00925-y/MediaObjects/43587_2025_925_MOESM1_ESM.xlsx"
)
FINNGEN_URL = (
    "https://www.finngen.fi/sites/default/files/inline-files/"
    "FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx"
)
EXPECTED_SHA256 = {
    "hfrs": "123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625",
    "finngen": "0004624f727b96b4a56c6fac75ff182eb8bc3daf78cdf690bdad276b31fca6c7",
}
ICD10_FIELDS = ("OUTPAT_ICD", "HD_ICD_10_ATC", "HD_ICD_10", "COD_ICD_10")
TOKEN_SPLIT = re.compile(r"[|,&;\s]+")


def _column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if not letters:
        raise ValueError(f"Invalid XLSX cell reference: {reference!r}")
    index = 0
    for char in letters.group(0):
        index = index * 26 + ord(char) - ord("A") + 1
    return index - 1


def read_xlsx_sheet(path: Path, sheet_name: str) -> tuple[list[list[str]], str, int]:
    """Read one worksheet via XML; do not trust the optional dimension element."""
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall("pkg:Relationship", NS)
        }
        sheets = workbook.findall("main:sheets/main:sheet", NS)
        selected = next((s for s in sheets if s.attrib.get("name") == sheet_name), None)
        if selected is None:
            raise ValueError(f"Worksheet {sheet_name!r} not found in {path}")
        target = targets[selected.attrib[f"{{{NS['rel']}}}id"]]
        if target.startswith("/"):
            worksheet_path = target.lstrip("/")
        else:
            worksheet_path = posixpath.normpath(posixpath.join("xl", target))

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(part.text or "" for part in item.findall(".//main:t", NS))
                for item in shared_root.findall("main:si", NS)
            ]

        root = ET.fromstring(archive.read(worksheet_path))
        data: list[list[str]] = []
        for row in root.findall(".//main:sheetData/main:row", NS):
            cells: dict[int, str] = {}
            for cell in row.findall("main:c", NS):
                kind = cell.attrib.get("t")
                if kind == "s":
                    value = cell.findtext("main:v", default="", namespaces=NS)
                    value = shared_strings[int(value)] if value else ""
                elif kind == "inlineStr":
                    value = "".join(p.text or "" for p in cell.findall(".//main:t", NS))
                else:
                    value = cell.findtext("main:v", default="", namespaces=NS)
                cells[_column_index(cell.attrib["r"])] = value
            if cells:
                row_values = [""] * (max(cells) + 1)
                for index, value in cells.items():
                    row_values[index] = value
                data.append(row_values)
            else:
                data.append([])
    return data, digest, path.stat().st_size


def _get(row: list[str], index: int) -> str:
    return row[index] if index < len(row) else ""


def crosswalk(hfrs_path: Path, finngen_path: Path, output: Path) -> dict[str, object]:
    hfrs_rows, hfrs_hash, hfrs_bytes = read_xlsx_sheet(hfrs_path, "ST18")
    finngen_rows, finngen_hash, finngen_bytes = read_xlsx_sheet(finngen_path, "Sheet 1")
    if hfrs_hash != EXPECTED_SHA256["hfrs"]:
        raise ValueError(f"Unexpected publisher workbook SHA-256: {hfrs_hash}")
    if finngen_hash != EXPECTED_SHA256["finngen"]:
        raise ValueError(f"Unexpected FinnGen workbook SHA-256: {finngen_hash}")

    headers = finngen_rows[0]
    required = ("TAGS", "NAME", "LONGNAME", *ICD10_FIELDS)
    missing = [name for name in required if name not in headers]
    if missing:
        raise ValueError(f"FinnGen endpoint sheet lacks columns: {missing}")
    positions = {name: headers.index(name) for name in required}
    endpoints: list[dict[str, str]] = []
    for row in finngen_rows[1:]:
        if not row:
            continue
        endpoint = {
            "tags": _get(row, positions["TAGS"]),
            "name": _get(row, positions["NAME"]),
            "longname": _get(row, positions["LONGNAME"]),
            "tokens": set(),
            "fields": {},
        }
        for field in ICD10_FIELDS:
            raw = _get(row, positions[field])
            endpoint["fields"][field] = raw  # type: ignore[assignment]
            endpoint["tokens"].update(t.strip() for t in TOKEN_SPLIT.split(raw) if t.strip())  # type: ignore[union-attr]
        if endpoint["name"]:
            endpoints.append(endpoint)

    code_rows: list[dict[str, str]] = []
    for row in hfrs_rows[2:]:
        if not row or not _get(row, 0).strip().isdigit():
            continue
        code = _get(row, 2).strip()
        description = _get(row, 1).replace("\xa0", " ").strip()
        raw_weight = _get(row, 3).strip()
        weight = format(float(raw_weight), ".15g") if raw_weight else ""
        matches = [e for e in endpoints if code in e["tokens"]]  # type: ignore[operator]
        code_rows.append(
            {
                "hfrs_row": _get(row, 0).strip(),
                "icd10_code": code,
                "description": description,
                "weight": weight,
                "literal_token_match_count": str(len(matches)),
                "matched_finngen_tags": "|".join(e["tags"] for e in matches)
                or "NO_LITERAL_EXACT_TOKEN",
                "matched_finngen_names": "|".join(e["name"] for e in matches)
                or "NO_LITERAL_EXACT_TOKEN",
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(code_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(code_rows)
    matched = sum(int(row["literal_token_match_count"]) > 0 for row in code_rows)
    return {
        "hfrs_source_url": HFRS_URL,
        "hfrs_bytes": hfrs_bytes,
        "hfrs_sha256": hfrs_hash,
        "finngen_source_url": FINNGEN_URL,
        "finngen_bytes": finngen_bytes,
        "finngen_sha256": finngen_hash,
        "finngen_endpoint_rows_with_names": len(endpoints),
        "hfrs_code_rows": len(code_rows),
        "hfrs_codes_with_literal_exact_token_match": matched,
        "hfrs_codes_without_literal_exact_token_match": len(code_rows) - matched,
        "definition": "Literal exact token matches only in OUTPAT_ICD, HD_ICD_10_ATC, HD_ICD_10, COD_ICD_10; encoded ranges/patterns and hierarchical expansion are not interpreted.",
        "output": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hfrs_supplement_xlsx", type=Path)
    parser.add_argument("finngen_df13_endpoint_xlsx", type=Path)
    parser.add_argument("output_tsv", type=Path)
    args = parser.parse_args()
    for key, value in crosswalk(args.hfrs_supplement_xlsx, args.finngen_df13_endpoint_xlsx, args.output_tsv).items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
