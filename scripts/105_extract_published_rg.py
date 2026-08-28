#!/usr/bin/env python3
"""Download checksum-pinned supplements and extract comparable published rg rows."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
import re
import ssl
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET


PRIMARY = "PRIMARY_PHASE1"
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL = re.compile(r"([A-Z]+)([0-9]+)$")
ALIASES = {
    "mdd": [("Major depressive disorder", "EXACT_CONCEPT"), ("Major depression", "EXACT_CONCEPT"), ("Depressive symptoms", "RELATED_DEPRESSION_PHENOTYPE")],
    "scz": [("Schizophrenia", "EXACT_CONCEPT")],
    "bipolar": [("Bipolar disorder", "EXACT_CONCEPT")],
    "adhd": [("Attention deficit hyperactivity disorder", "EXACT_CONCEPT"), ("ADHD", "EXACT_CONCEPT")],
    "ibd": [("Inflammatory Bowel Disease", "EXACT_CONCEPT"), ("Inflammatory Bowel Disease (Euro)", "EXACT_CONCEPT")],
    "crohn": [("Crohn's disease", "EXACT_CONCEPT"), ("Crohns disease", "EXACT_CONCEPT")],
    "ra": [("Rheumatoid Arthritis", "EXACT_CONCEPT")],
    "asthma": [("Asthma", "EXACT_CONCEPT")],
    "bmi": [("Body mass index", "EXACT_CONCEPT"), ("Body mass index (BMI)", "EXACT_CONCEPT")],
    "ldl": [("LDL cholesterol", "EXACT_CONCEPT")],
    "hdl": [("HDL cholesterol", "EXACT_CONCEPT")],
    "triglycerides": [("Triglycerides", "EXACT_CONCEPT"), ("Serum total triglycerides", "RELATED_LIPID_PHENOTYPE")],
    "cad": [("Coronary artery disease", "EXACT_CONCEPT"), ("Coronary heart disease", "EXACT_CONCEPT")],
    "stroke": [("Stroke", "EXACT_CONCEPT"), ("Stroke (any)", "EXACT_CONCEPT")],
    "atrial_fibrillation": [("Atrial fibrillation", "EXACT_CONCEPT"), ("Diagnoses - main ICD10: I48 Atrial fibrillation and flutter", "RELATED_ENDPOINT")],
    "sbp": [("Systolic blood pressure", "EXACT_CONCEPT")],
    "longevity": [("Longevity", "EXACT_CONCEPT")],
    "parental_lifespan": [("Parents age at death", "RELATED_SOURCE_DEFINITION"), ("Longevity", "RELATED_SOURCE_DEFINITION")],
    "healthspan": [("Healthspan", "EXACT_CONCEPT")],
    "frailty": [("Frailty", "EXACT_CONCEPT"), ("Frailty index", "EXACT_CONCEPT")],
    "telomere_length": [("Leukocyte telomere length", "EXACT_CONCEPT"), ("Telomere length", "EXACT_CONCEPT")],
    "grip_strength": [("Hand grip strength", "EXACT_CONCEPT"), ("Grip strength", "EXACT_CONCEPT")],
    "breast_cancer": [("Breast cancer", "EXACT_CONCEPT"), ("Breast cancer (overall)", "EXACT_CONCEPT")],
    "prostate_cancer": [("Prostate cancer", "EXACT_CONCEPT"), ("Prostate cancer (overall)", "EXACT_CONCEPT")],
    "colorectal_cancer": [("Colorectal cancer", "EXACT_CONCEPT"), ("Bowel cancer", "RELATED_ENDPOINT")],
    "lung_cancer": [("Lung cancer (all)", "EXACT_CONCEPT"), ("Lung cancer", "EXACT_CONCEPT"), ("Lung adenocarcinoma", "RELATED_HISTOLOGY"), ("Squamous cell lung cancer", "RELATED_HISTOLOGY"), ("Lung cancer (squamous cell)", "RELATED_HISTOLOGY")],
    "ovarian_cancer": [("Ovarian cancer", "EXACT_CONCEPT")],
    "parkinson": [("Parkinsons disease", "EXACT_CONCEPT"), ("Parkinson's disease", "EXACT_CONCEPT")],
}
LAYOUTS = {
    "jansen_2019_supp_table_29": {
        "sheet": "S29_Rg_previous_studies", "header": 4,
        "sleep": {"insomnia": (4, 5, 7)}, "trait": 0, "pmid": 3, "sample_size": 1,
        "published_sleep_n": 1331010,
    },
    "jones_2019_chronotype_supp_data_11": {
        "sheet": "Supp Data 11", "header": 2,
        "sleep": {"chronotype": (4, 5, 7)}, "trait": 0, "pmid": 1, "sample_size": None,
        "published_sleep_n": 697828,
    },
    "dashti_2019_duration_supp_data_18": {
        "sheet": "Supplementary Data 18", "header": 3,
        "sleep": {"sleepdur": (3, 4, 5), "shortsleep": (9, 10, 11), "longsleep": (12, 13, 14)},
        "trait": 1, "pmid": 2, "sample_size": None, "published_sleep_n": 446118,
    },
    "wang_2019_sleepiness_supp_data_6": {
        "sheet": "Sheet1", "header": 3,
        "sleep": {"sleepiness": (6, 7, 8)}, "trait": 2, "pmid": 3, "sample_size": None,
        "published_sleep_n": 452071,
    },
    "jones_2019_actigraphy_supp_data_7": {
        "sheet": "Supp Data 7", "header": 3,
        "sleep": {
            "sleep_efficiency": (5, 6, 8, "Sleep efficiency"),
            "accel_sleep_duration": (5, 6, 8, "Sleep duration"),
            "sleep_timing": (5, 6, 8, "Sleep midpoint"),
        },
        "trait": 1, "pmid": 2, "sample_size": None, "sleep_trait_column": 0,
        "published_sleep_n": 85449,
    },
    "campos_2020_snoring_supp_data_2": {
        "sheet": "rG", "header": 0,
        "sleep": {"snoring": (5, 6, 7)}, "trait": 1, "pmid": 2, "sample_size": None,
        "published_sleep_n": 408317,
    },
    "dashti_2021_napping_supp_data_10": {
        "sheet": "Sheet1", "header": 3,
        "sleep": {"napping": (4, 5, 6)}, "trait": 1, "pmid": 2, "sample_size": None,
        "published_sleep_n": 452633,
    },
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def tls_context() -> ssl.SSLContext:
    for candidate in (os.environ.get("SSL_CERT_FILE"), "/etc/ssl/cert.pem"):
        if candidate and Path(candidate).is_file():
            return ssl.create_default_context(cafile=candidate)
    return ssl.create_default_context()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire(source: dict[str, str], path: Path, refresh: bool) -> None:
    expected_size = int(source["bytes"])
    expected_hash = source["sha256"]
    if path.is_file() and path.stat().st_size == expected_size and sha256(path) == expected_hash:
        return
    if not refresh:
        fail(f"missing or invalid cache file (rerun with --refresh): {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    request = Request(source["download_url"], headers={"User-Agent": "sleep-gwas-atlas/1.0"})
    temporary = path.with_suffix(path.suffix + ".tmp")
    with urlopen(request, timeout=120, context=tls_context()) as response, temporary.open("wb") as handle:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            handle.write(block)
    if temporary.stat().st_size != expected_size or sha256(temporary) != expected_hash:
        temporary.unlink(missing_ok=True)
        fail(f"download identity mismatch: {source['evidence_source_id']}")
    temporary.replace(path)


def column_index(reference: str) -> int:
    match = CELL.match(reference)
    if not match:
        fail(f"invalid XLSX cell reference: {reference}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - 64
    return value - 1


def xlsx_sheet(path: Path, sheet_name: str) -> list[list[object]]:
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{MAIN_NS}}}si"):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t")))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {node.attrib["Id"]: node.attrib["Target"] for node in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")}
        target = None
        for node in workbook.find(f"{{{MAIN_NS}}}sheets"):
            if node.attrib["name"] == sheet_name:
                target = targets[node.attrib[f"{{{REL_NS}}}id"]]
                break
        if target is None:
            fail(f"worksheet not found in {path.name}: {sheet_name}")
        member = target.lstrip("/") if target.startswith("/") else "xl/" + target.lstrip("./")
        root = ET.fromstring(archive.read(member))
        rows = []
        for row_node in root.findall(f".//{{{MAIN_NS}}}sheetData/{{{MAIN_NS}}}row"):
            values: dict[int, object] = {}
            for cell in row_node.findall(f"{{{MAIN_NS}}}c"):
                index = column_index(cell.attrib["r"])
                kind = cell.attrib.get("t")
                value_node = cell.find(f"{{{MAIN_NS}}}v")
                if kind == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t"))
                elif value_node is None:
                    value = ""
                elif kind == "s":
                    value = shared[int(value_node.text)]
                elif kind in {"str", "e"}:
                    value = value_node.text or ""
                elif kind == "b":
                    value = value_node.text == "1"
                else:
                    raw = value_node.text or ""
                    try:
                        number = float(raw)
                        value = int(number) if number.is_integer() else number
                    except ValueError:
                        value = raw
                values[index] = value
            width = max(values, default=-1) + 1
            rows.append([values.get(index, "") for index in range(width)])
        return rows


def value(row: list[object], index: int | None) -> object:
    if index is None or index >= len(row) or row[index] == "":
        return "NA"
    return row[index]


def finite_number(item: object) -> float | None:
    try:
        result = float(item)
    except (TypeError, ValueError):
        return None
    return result if result == result and abs(result) != float("inf") else None


def normalized(item: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(item).casefold()).strip()


def evidence_rows(
    source: dict[str, str], path: Path, significant: set[tuple[str, str]],
    panel_pmid: dict[str, str],
) -> list[dict[str, object]]:
    layout = LAYOUTS[source["evidence_source_id"]]
    sheet = xlsx_sheet(path, layout["sheet"])
    data = sheet[layout["header"] + 1:]
    rows = []
    for sleep_trait, coordinates in layout["sleep"].items():
        rg_col, se_col, p_col, *sleep_filter = coordinates
        for external_trait, aliases in ALIASES.items():
            if (sleep_trait, external_trait) not in significant:
                continue
            alias_map = {normalized(label): (rank, match) for rank, (label, match) in enumerate(aliases)}
            for row_number, item in enumerate(data, layout["header"] + 2):
                if sleep_filter and normalized(value(item, layout["sleep_trait_column"])) != normalized(sleep_filter[0]):
                    continue
                trait_label = value(item, layout["trait"])
                match = alias_map.get(normalized(trait_label))
                if match is None:
                    continue
                rg = finite_number(value(item, rg_col))
                se = finite_number(value(item, se_col))
                p = finite_number(value(item, p_col))
                if rg is None or se is None or se <= 0 or p is None:
                    continue
                trait_pmid_raw = value(item, layout["pmid"])
                trait_pmid = str(int(trait_pmid_raw)) if finite_number(trait_pmid_raw) not in (None, 0) else "NA"
                current_pmid = panel_pmid.get(external_trait, "NA")
                external_dataset_relation = (
                    "SAME_PUBLICATION_IDENTIFIER" if trait_pmid != "NA" and trait_pmid == current_pmid
                    else "DIFFERENT_OR_UNRESOLVED_EXTERNAL_DATASET"
                )
                sleep_relation = (
                    "DIFFERENT_OR_EXPANDED_SLEEP_GWAS" if sleep_trait in {"insomnia", "chronotype"}
                    else "SAME_OR_OVERLAPPING_SLEEP_GWAS"
                )
                if "SAME_OR_OVERLAPPING" in sleep_relation or external_dataset_relation == "SAME_PUBLICATION_IDENTIFIER":
                    comparison_independence = "NOT_FULLY_INDEPENDENT"
                elif "EXPANDED" in sleep_relation:
                    comparison_independence = "PARTIALLY_INDEPENDENT_OR_OVERLAPPING"
                else:
                    comparison_independence = "DIFFERENT_DATASET_COMPARISON"
                rows.append({
                    "sleep_trait": sleep_trait,
                    "external_trait": external_trait,
                    "evidence_source_id": source["evidence_source_id"],
                    "paper_pmid": source["pmid"],
                    "paper_doi": source["doi"],
                    "article_title": source["article_title"],
                    "source_table": source["source_table"],
                    "source_row_number": row_number,
                    "published_trait_label": trait_label,
                    "external_phenotype_match": match[1],
                    "published_trait_pmid": trait_pmid,
                    "published_trait_sample_size": value(item, layout["sample_size"]),
                    "published_sleep_sample_size": layout["published_sleep_n"],
                    "published_rg": rg,
                    "published_rg_se": se,
                    "published_rg_p": p,
                    "sleep_dataset_relation": sleep_relation,
                    "external_dataset_relation": external_dataset_relation,
                    "comparison_independence": comparison_independence,
                    "source_url": source["download_url"],
                    "source_sha256": source["sha256"],
                    "extraction_status": "REVIEWED_TABLE_LEVEL_DIRECT_RG",
                })
    unique = {}
    for row in rows:
        key = (
            row["sleep_trait"], row["external_trait"], row["evidence_source_id"],
            row["source_row_number"], row["published_rg"],
        )
        unique[key] = row
    return list(unique.values())


def tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    sources = read_tsv(root / "config/published_rg_source_workbooks.tsv")
    if len(sources) != 7 or {row["review_status"] for row in sources} != {"REVIEWED_TABLE_LEVEL"}:
        fail("expected seven reviewed published-rg source workbooks")
    master = read_tsv(root / "results/analysis/phase1_master_analysis.tsv")
    significant = {
        (row["sleep_trait"], row["external_trait"])
        for row in master
        if row["primary_or_sensitivity"] == PRIMARY and float(row["fdr"]) <= .05
    }
    if len(significant) != 153:
        fail(f"expected 153 primary discoveries, observed {len(significant)}")
    panel = read_tsv(root / "config/analysis_panel.tsv")
    panel_pmid = {row["trait_id"]: row["pmid"] for row in panel}
    cache = root / "results/_literature_cache/published_rg_workbooks"
    all_rows = []
    for source in sources:
        if source["evidence_source_id"] not in LAYOUTS:
            fail(f"missing extraction layout: {source['evidence_source_id']}")
        path = cache / source["filename"]
        acquire(source, path, args.refresh)
        all_rows.extend(evidence_rows(source, path, significant, panel_pmid))
    all_rows.sort(key=lambda row: (
        row["sleep_trait"], row["external_trait"], row["evidence_source_id"],
        int(row["source_row_number"]),
    ))
    fields = [
        "sleep_trait", "external_trait", "evidence_source_id", "paper_pmid", "paper_doi",
        "article_title", "source_table", "source_row_number", "published_trait_label",
        "external_phenotype_match", "published_trait_pmid", "published_trait_sample_size",
        "published_sleep_sample_size", "published_rg", "published_rg_se", "published_rg_p",
        "sleep_dataset_relation", "external_dataset_relation", "comparison_independence",
        "source_url", "source_sha256", "extraction_status",
    ]
    payload = tsv_text(all_rows, fields)
    output = root / "results/analysis/published_rg_source_evidence.tsv"
    if args.validate_only:
        if not output.is_file() or output.read_text(encoding="utf-8") != payload:
            fail("published-rg extraction drifted")
        print(f"PUBLISHED_RG_EVIDENCE_VALID rows={len(all_rows)} pairs={len({(r['sleep_trait'], r['external_trait']) for r in all_rows})}")
        return 0
    atomic_text(output, payload)
    print(f"PUBLISHED_RG_EVIDENCE_BUILT rows={len(all_rows)} pairs={len({(r['sleep_trait'], r['external_trait']) for r in all_rows})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
