#!/usr/bin/env python3
"""Match Fan et al. 2026 sleep–disease LDSC rows to inherited Brain6 hits.

Uses only Python's standard library. Source workbooks are read-only inputs and
are not copied into the repository. The output is a source-match audit, not a
replication family or a replacement for the frozen novelty crosswalk.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET

NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p": "http://schemas.openxmlformats.org/package/2006/relationships",
}
BASE = "https://static-content.springer.com/esm/art%3A10.1038%2Fs43856-026-01656-w/MediaObjects/"
SLEEP = {
    "insomnia": ("Insomnia - Jansen et al., 2019", "jansen_2019_insomnia_ukb"),
    "sleepdur": ("Sleep duration - Dashti et al., 2019", "dashti_2019_sleep_duration"),
    "shortsleep": ("Short sleep duration - Dashti et al., 2019", "dashti_2019_short_sleep"),
    "longsleep": ("Long sleep duration - Dashti et al., 2019", "dashti_2019_long_sleep"),
    "chronotype": ("", "jones_2019_morning_person_ukb"),
    "sleepiness": ("Narcolepsy - Wang et al., 2019", "wang_2019_daytime_sleepiness"),
    "napping": ("Daytime nap - Dashti et al., 2021", "dashti_2021_daytime_napping"),
    "snoring": ("Snoring - Campos et al., 2020", "campos_2020_snoring"),
    "accel_sleep_duration": (
        "Accelerometer based sleep duration - Jones et al., 2019",
        "jones_2019_accelerometer_sleep",
    ),
}
SLEEPY_METADATA = {
    "jansen_2019_insomnia_ukb": ("Jansen PR", "Insomnia", 386533),
    "dashti_2019_sleep_duration": ("Dashti HS", "Sleep duration", 446118),
    "dashti_2019_short_sleep": ("Dashti HS", "Short sleep duration", 411934),
    "dashti_2019_long_sleep": ("Dashti HS", "Long sleep duration", 339926),
    "wang_2019_daytime_sleepiness": ("Wang H", "Self-report (daytime sleepiness associations)", 452071),
    "dashti_2021_daytime_napping": ("Dashti HS", "Daytime nap", 452633),
    "campos_2020_snoring": ("Campos AI", "Snoring", 408317),
    "jones_2019_accelerometer_sleep": ("Jones SE", "Sleep duration", 85449),
    "jones_2019_morning_person_ukb": ("Jones SE", "Chronotype", 403195),
}
DISEASE = {
    "bipolar": ("Bipolar disorder - Mullins et al., 2021", "mullins_2021_bipolar_eur", "34002096"),
    "scz": ("Schizophrenia - Trubetskoy et al., 2022", "trubetskoy_2022_scz_eur", "35396580"),
    "parkinson": ("Parkinson's disease - Nalls et al., 2019", "nalls_2019_parkinson_public_proxy", "31701892"),
}
CHRONOTYPE_LABELS = (
    "Chronotype - Jones et al., 2019a",
    "Chronotype - Jones et al., 2019b",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref).group(0)  # type: ignore[union-attr]
    value = 0
    for char in letters:
        value = value * 26 + ord(char) - 64
    return value - 1


def read_xlsx_first_sheet(path: Path) -> list[list[str]]:
    """Read a plain, first-sheet XLSX table without third-party packages."""
    with zipfile.ZipFile(path) as book:
        strings: list[str] = []
        if "xl/sharedStrings.xml" in book.namelist():
            root = ET.fromstring(book.read("xl/sharedStrings.xml"))
            strings = ["".join(t.text or "" for t in si.findall(".//m:t", NS))
                       for si in root.findall("m:si", NS)]
        sheet_names = [n for n in book.namelist()
                       if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n)]
        if not sheet_names:
            raise ValueError(f"No worksheet XML found in {path}")
        root = ET.fromstring(book.read(sorted(sheet_names)[0]))
        rows: list[list[str]] = []
        for row in root.findall(".//m:sheetData/m:row", NS):
            cells: dict[int, str] = {}
            for cell in row.findall("m:c", NS):
                ref = cell.attrib.get("r", "")
                if not ref:
                    continue
                kind = cell.attrib.get("t")
                if kind == "inlineStr":
                    value = "".join(t.text or "" for t in cell.findall(".//m:t", NS))
                else:
                    node = cell.find("m:v", NS)
                    value = node.text if node is not None and node.text else ""
                    if kind == "s" and value:
                        value = strings[int(value)]
                cells[_col_index(ref)] = value
            if cells:
                width = max(cells) + 1
                rows.append([cells.get(i, "") for i in range(width)])
        return rows


def records(path: Path, header_row: int = 1) -> list[dict[str, str]]:
    rows = read_xlsx_first_sheet(path)
    if len(rows) <= header_row:
        raise ValueError(f"Workbook is missing expected table header: {path}")
    headers = [x.strip() for x in rows[header_row]]
    return [dict(zip(headers, [v.strip() for v in row] + [""] * max(0, len(headers)-len(row))))
            for row in rows[header_row + 1:] if any(v.strip() for v in row)]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _normalize_label(value: str) -> str:
    """Normalize inconsistent whitespace in publisher workbook labels."""
    return re.sub(r"\s+", " ", value or "").strip()


def build_context(global_tsv: Path, data2: Path, data4: Path, data11: Path) -> list[dict[str, str]]:
    sleep_metadata = records(data2)
    disease_metadata = records(data4)
    rg_rows = records(data11)
    _require("Supplementary Data 2" in read_xlsx_first_sheet(data2)[0][0], "Data2 workbook title mismatch")
    _require("Supplementary Data 4" in read_xlsx_first_sheet(data4)[0][0], "Data4 workbook title mismatch")
    _require("Supplementary Data 11" in read_xlsx_first_sheet(data11)[0][0], "Data11 workbook title mismatch")

    data2_by_key = {(_normalize_label(r.get("Author", "")),
                     _normalize_label(r.get("Phenotype", "")),
                     _normalize_label(r.get("Sample size", ""))): r
                    for r in sleep_metadata}
    data4_by_label = {_normalize_label(r.get("Disorder/trait", "")): r for r in disease_metadata}
    rg_by_pair = {(_normalize_label(r.get("Sleep trait", "")),
                   _normalize_label(r.get("Disease/trait", ""))): r
                  for r in rg_rows}

    # Confirm each locked Brain6 source ID is represented by the expected study
    # and sample size in Fan's cited sleep-GWAS inventory.
    for source_id, (author, phenotype, n) in SLEEPY_METADATA.items():
        _require((author, phenotype, str(n)) in data2_by_key,
                 f"Fan Data2 does not verify expected sleep source {source_id}")
    for label, _source_id, pmid in DISEASE.values():
        _require(_normalize_label(label) in data4_by_label and
                 data4_by_label[_normalize_label(label)].get("PMID") == pmid,
                 f"Fan Data4 does not verify expected disease source {label}")

    with global_tsv.open(encoding="utf-8", newline="") as stream:
        global_rows = list(csv.DictReader(stream, delimiter="\t"))
    hits = [r for r in global_rows if r["significance_under_original_396_family"].lower() == "true"]
    _require(len(hits) == 35, f"Expected 35 inherited Brain6 hits, found {len(hits)}")
    output: list[dict[str, str]] = []
    for hit in hits:
        source_id = hit["source_GWAS_ID_sleep"]
        if hit["sleep_trait"] in SLEEP and source_id in SLEEPY_METADATA:
            expected_n = SLEEPY_METADATA[source_id][2]
            _require(int(hit["sleep_sample_size"]) == expected_n,
                     f"Brain6 sample size does not match Fan Data2 source for {source_id}")
        sleep_label = SLEEP.get(hit["sleep_trait"], ("", ""))[0]
        sleep_source = SLEEP.get(hit["sleep_trait"], ("", ""))[1]
        disease = DISEASE.get(hit["brain_disorder"])
        disease_label, disease_source, _ = disease if disease else ("", "", "")
        candidate_labels: Iterable[str] = ()
        match_status = "NO_EXACT_STUDY_PAIR"
        if hit["sleep_trait"] == "chronotype" and disease:
            candidate_labels = CHRONOTYPE_LABELS
            match_status = "SOURCE_SUFFIX_UNRESOLVED"
        elif sleep_label and disease_label and (sleep_label, disease_label) in rg_by_pair:
            candidate_labels = (sleep_label,)
            match_status = "EXACT_SOURCE_PAIR"
        else:
            candidate_labels = (sleep_label,) if sleep_label else ()

        found = [rg_by_pair[(_normalize_label(s), _normalize_label(disease_label))]
                 for s in candidate_labels
                 if disease_label and (_normalize_label(s), _normalize_label(disease_label)) in rg_by_pair]
        if match_status == "EXACT_SOURCE_PAIR" and len(found) != 1:
            match_status = "NO_EXACT_STUDY_PAIR"
            found = []
        for fan in found:
            fan_rg = float(fan["rg"])
            brain_rg = float(hit["rg"])
            output.append({
                "sleep_trait": hit["sleep_trait"],
                "brain_disorder": hit["brain_disorder"],
                "brain6_rg": hit["rg"],
                "brain6_original_q": hit["original_BH_FDR_q"],
                "sleep_source_id": sleep_source,
                "fan_sleep_source": fan["Sleep trait"],
                "disease_source_id": disease_source,
                "fan_disease_source": fan["Disease/trait"],
                "source_match_status": match_status,
                "fan_rg": fan["rg"],
                "fan_se": fan["se"],
                "fan_p": fan["P_value"],
                "fan_p_fdr": fan["P_FDR"],
                "fan_fdr_005": str(float(fan["P_FDR"]) < 0.05).lower(),
                "direction_concordant": str(fan_rg * brain_rg > 0).lower(),
                "interpretation": "prior direct LDSC context; source data reused, not independent replication",
            })
    output.sort(key=lambda r: (r["sleep_trait"], r["brain_disorder"], r["fan_sleep_source"]))
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--global-map", type=Path, required=True)
    parser.add_argument("--data2", type=Path, required=True, help="Fan Supplementary Data 2 XLSX")
    parser.add_argument("--data4", type=Path, required=True, help="Fan Supplementary Data 4 XLSX")
    parser.add_argument("--data11", type=Path, required=True, help="Fan Supplementary Data 11 XLSX")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    rows = build_context(args.global_map, args.data2, args.data4, args.data11)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with args.output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    exact = [r for r in rows if r["source_match_status"] == "EXACT_SOURCE_PAIR"]
    ambiguous = [r for r in rows if r["source_match_status"] == "SOURCE_SUFFIX_UNRESOLVED"]
    inputs = {}
    for key, path in (("brain6_locked_global_map", args.global_map),
                      ("fan_supplementary_data_2", args.data2),
                      ("fan_supplementary_data_4", args.data4),
                      ("fan_supplementary_data_11", args.data11)):
        item = {"downloaded_basename": path.name, "bytes": path.stat().st_size,
                "sha256": sha256(path)}
        if key.startswith("fan_"):
            num = key.rsplit("_", 1)[1]
            item["source_url"] = f"{BASE}43856_2026_1656_MOESM{int(num)+1}_ESM.xlsx"
        inputs[key] = item
    provenance = {
        "schema_version": 1,
        "status": "PASS_SOURCE_MATCH_AUDIT_NOT_INDEPENDENT_REPLICATION",
        "citation": "Fan Z et al. Communications Medicine. 2026;6:435. doi:10.1038/s43856-026-01656-w",
        "article_url": "https://www.nature.com/articles/s43856-026-01656-w",
        "supplement_license": "CC BY-NC-ND 4.0 (source workbooks retained only as unmodified temporary inputs; not copied into this repository)",
        "generated_utc_date": "2026-09-26",
        "method": "Match locked significant Brain6 pairs to Fan Data11 only where sleep and disease source labels are source-identical; chronotype suffixes remain ambiguous.",
        "inputs": inputs,
        "summary": {
            "brain6_inherited_significant_pairs": len([r for r in csv.DictReader(args.global_map.open(encoding="utf-8"), delimiter="\t") if r["significance_under_original_396_family"].lower()=="true"]),
            "source_exact_pairs": len(exact),
            "exact_pairs_passing_fan_p_fdr_005": sum(r["fan_fdr_005"] == "true" for r in exact),
            "exact_pairs_direction_concordant": sum(r["direction_concordant"] == "true" for r in exact),
            "chronotype_source_suffix_unresolved_rows": len(ambiguous),
            "chronotype_rows_passing_fan_p_fdr_005": sum(r["fan_fdr_005"] == "true" for r in ambiguous),
            "independent_replication": False,
            "novelty_classifications_changed": False,
        },
        "output": {"path": str(args.output), "sha256": sha256(args.output)},
        "limitations": [
            "Source overlap is intrinsic because Fan et al. reused published GWAS; these comparisons are literature context, not independent replication.",
            "The two Jones 2019 chronotype rows share a source paper but the a/b suffix-to-release mapping is not resolved from the supplement metadata; those rows are not called exact source pairs.",
            "Pairs without a matched source-pair row are not evidence that no related published correlation exists.",
        ],
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(provenance["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
