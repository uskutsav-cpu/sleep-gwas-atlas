from __future__ import annotations

import csv
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from audit_fan_2026_rg_context import (
    DISEASE,
    SLEEP,
    SLEEPY_METADATA,
    build_context,
    read_xlsx_first_sheet,
)


def _column_name(index: int) -> str:
    out = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        out = chr(65 + remainder) + out
    return out


def _write_minimal_xlsx(path: Path, rows: list[list[str]]) -> None:
    xml_rows = []
    for row_number, row in enumerate(rows, start=1):
        cells = []
        for column, value in enumerate(row):
            ref = f"{_column_name(column)}{row_number}"
            escaped = (value.replace("&", "&amp;").replace("<", "&lt;")
                       .replace(">", "&gt;").replace('"', "&quot;"))
            cells.append(
                f'<c r="{ref}" t="inlineStr"><is><t>{escaped}</t></is></c>'
            )
        xml_rows.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    sheet = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>'
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as workbook:
        workbook.writestr("xl/worksheets/sheet1.xml", sheet)


def test_xlsx_reader_handles_inline_cells_and_sparse_columns(tmp_path: Path) -> None:
    path = tmp_path / "small.xlsx"
    _write_minimal_xlsx(path, [["title"], ["label", "value"], ["study", "3"]])
    assert read_xlsx_first_sheet(path) == [
        ["title"], ["label", "value"], ["study", "3"]
    ]


def test_build_context_requires_exact_source_pairs_and_keeps_chronotype_ambiguous(
    tmp_path: Path,
) -> None:
    data2 = tmp_path / "data2.xlsx"
    data4 = tmp_path / "data4.xlsx"
    data11 = tmp_path / "data11.xlsx"
    global_tsv = tmp_path / "global.tsv"

    _write_minimal_xlsx(data2, [
        ["Supplementary Data 2: sleep GWAS"],
        ["Author", "Phenotype", "Sample size"],
        *[[author, phenotype, str(n)] for author, phenotype, n in SLEEPY_METADATA.values()],
    ])
    _write_minimal_xlsx(data4, [
        ["Supplementary Data 4: disease GWAS"],
        ["Disorder/trait", "PMID"],
        *[[label, pmid] for label, _source_id, pmid in DISEASE.values()],
    ])

    exact_pairs = [
        ("insomnia", "bipolar"),
        ("sleepdur", "scz"), ("sleepdur", "bipolar"),
        ("longsleep", "scz"), ("longsleep", "bipolar"),
        ("longsleep", "parkinson"),
        ("sleepiness", "scz"), ("sleepiness", "bipolar"),
        ("napping", "scz"), ("napping", "bipolar"),
        ("snoring", "scz"),
        ("accel_sleep_duration", "bipolar"),
    ]
    disease_label = {key: value[0] for key, value in DISEASE.items()}
    article_rows: list[list[str]] = [[
        "Supplementary Data 11: genetic correlations",
    ], ["Sleep trait", "Disease/trait", "rg", "se", "P_value", "P_FDR", "Category"]]
    for sleep, disorder in exact_pairs:
        sleep_label = SLEEP[sleep][0]
        label = disease_label[disorder]
        if disorder == "parkinson":
            label = label.replace(" - ", "  - ")
        article_rows.append([
            sleep_label, label, "-0.10" if disorder == "parkinson" else "0.10", "0.02", "0.001",
            "0.01", "Brain disorders",
        ])
    for sleep_label in (
        "Chronotype - Jones et al., 2019a",
        "Chronotype - Jones et al., 2019b",
    ):
        for disorder in ("scz", "bipolar"):
            article_rows.append([
                sleep_label, disease_label[disorder], "-0.10", "0.02", "0.001",
                "0.01", "Brain disorders",
            ])
    _write_minimal_xlsx(data11, article_rows)

    global_fields = [
        "sleep_trait", "brain_disorder", "rg", "original_BH_FDR_q",
        "significance_under_original_396_family", "source_GWAS_ID_sleep",
        "sleep_sample_size",
    ]
    global_rows: list[dict[str, str]] = []
    source_n = {source: n for source, (_author, _phenotype, n) in SLEEPY_METADATA.items()}
    for sleep, disorder in exact_pairs:
        source_id = SLEEP[sleep][1]
        global_rows.append({
            "sleep_trait": sleep, "brain_disorder": disorder,
            "rg": "-0.08" if disorder == "parkinson" else "0.08",
            "original_BH_FDR_q": "0.001",
            "significance_under_original_396_family": "True",
            "source_GWAS_ID_sleep": source_id, "sleep_sample_size": str(source_n[source_id]),
        })
    for disorder in ("scz", "bipolar"):
        global_rows.append({
            "sleep_trait": "chronotype", "brain_disorder": disorder, "rg": "-0.08",
            "original_BH_FDR_q": "0.001",
            "significance_under_original_396_family": "True",
            "source_GWAS_ID_sleep": "jones_2019_morning_person_ukb",
            "sleep_sample_size": "403195",
        })
    for index in range(21):
        global_rows.append({
            "sleep_trait": f"unmatched_{index}", "brain_disorder": "adhd", "rg": "0.1",
            "original_BH_FDR_q": "0.001",
            "significance_under_original_396_family": "True",
            "source_GWAS_ID_sleep": "not_in_fan_panel", "sleep_sample_size": "1",
        })
    with global_tsv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=global_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(global_rows)

    rows = build_context(global_tsv, data2, data4, data11)
    exact = [row for row in rows if row["source_match_status"] == "EXACT_SOURCE_PAIR"]
    ambiguous = [row for row in rows if row["source_match_status"] == "SOURCE_SUFFIX_UNRESOLVED"]
    assert len(rows) == 16
    assert len(exact) == 12
    assert len(ambiguous) == 4
    assert any(row["sleep_trait"] == "longsleep" and
               row["brain_disorder"] == "parkinson" for row in exact)
    assert all(row["direction_concordant"] == "true" for row in rows)
    assert all(row["fan_fdr_005"] == "true" for row in exact)
