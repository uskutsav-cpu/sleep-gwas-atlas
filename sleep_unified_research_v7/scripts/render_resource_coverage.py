#!/usr/bin/env python3
"""Render the preregistered V7 class counts as a compact vector figure."""
import csv
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "sleep_unified_research_v7" / "tables" / "cross_resource_match_summary.tsv"
OUT = ROOT / "sleep_unified_research_v7" / "figures" / "cross_resource_coverage.svg"
CLASSES = [
    ("A_exact_source_estimand", "A  Exact source/estimand", "#177e89"),
    ("B_related_nonidentical", "B  Related, nonidentical", "#4d8cc9"),
    ("C_incompatible", "C  Incompatible", "#6f7782"),
    ("D_absent_in_snapshot_scope", "D  Absent in snapshot scope", "#c8cdd3"),
    ("E_insufficient_evidence", "E  Insufficient evidence", "#e59b43"),
]
DISPLAY = {
    "human_gwas_atlas": "Human GWAS ATLAS R3",
    "morrison_2024": "Morrison 2024 S10–S15",
    "goodman_2025": "Goodman 2025 SD24",
    "fan_2026": "Fan 2026 SD11",
    "union_accessible_snapshots": "Union of row-level snapshots",
}


def text(x, y, value, size=18, color="#19212b", weight="normal", anchor="start"):
    return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" fill="{color}" font-weight="{weight}" text-anchor="{anchor}">{escape(str(value))}</text>'


def main():
    with TABLE.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    records = [r for r in rows if r["comparator"] != "union_accessible_snapshots"]
    records.append(next(r for r in rows if r["comparator"] == "union_accessible_snapshots"))
    assert len(records) == 5 and all(int(r["atlas_records"]) == 1637 for r in records)

    width, height = 1600, 870
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#ffffff"/>']
    parts += [text(70, 58, "Coverage of the 1,637-record atlas in accessible comparator snapshots", 27, weight="bold"),
              text(70, 91, "Source-aware match classes; absence is limited to each named release or table scope.", 17, color="#48515d")]

    legend_x = 70
    for _, label, color in CLASSES:
        parts.append(f'<rect x="{legend_x}" y="122" width="18" height="18" rx="2" fill="{color}"/>')
        parts.append(text(legend_x + 25, 137, label, 14, color="#28323d"))
        legend_x += 277 if "Incompatible" not in label else 220

    plot_x, plot_width = 390, 1050
    top, bar_h, gap = 198, 52, 39
    parts.append(f'<line x1="{plot_x}" y1="178" x2="{plot_x + plot_width}" y2="178" stroke="#c7cdd4" stroke-width="1"/>')
    for tick in (0, 400, 800, 1200, 1637):
        x = plot_x + plot_width * tick / 1637
        parts.append(f'<line x1="{x:.1f}" y1="182" x2="{x:.1f}" y2="{top + 5* (bar_h+gap) - gap}" stroke="#e8ebee" stroke-width="1"/>')
        parts.append(text(x, 174, f"{tick:,}", 12, color="#69727d", anchor="middle"))
    parts.append(text(plot_x + plot_width/2, 158, "Atlas records (n = 1,637)", 13, color="#48515d", anchor="middle"))

    for i, row in enumerate(records):
        y = top + i * (bar_h + gap)
        name = DISPLAY[row["comparator"]]
        parts.append(text(70, y + 32, name, 17, weight="bold" if i == 4 else "normal"))
        x = plot_x
        for field, _, color in CLASSES:
            n = int(row[field])
            segment = plot_width * n / 1637
            if n:
                parts.append(f'<rect x="{x:.2f}" y="{y}" width="{segment:.2f}" height="{bar_h}" fill="{color}"/>')
                if segment >= 31:
                    fg = "#ffffff" if color in {"#177e89", "#4d8cc9", "#6f7782", "#e59b43"} else "#26313c"
                    parts.append(text(x + segment/2, y + 32, f"{n:,}", 14, fg, "bold", "middle"))
                x += segment
        if row["comparator"] == "union_accessible_snapshots":
            q = int(row["qualified_match_absent_core_extension_across_four_snapshots"])
            label = f"{q:,}*"
        else:
            q = int(row["eligible_core_extension_without_A_B_E_in_snapshot"])
            label = f"{q:,}"
        parts.append(text(1530, y + 32, label, 15, "#26313c", "bold", "middle"))
    parts.append(text(1530, 177, "Qualified*", 13, "#48515d", "bold", "middle"))

    y0 = 687
    notes = [
        "No exact source/estimand match was established. Thirty HGA bibliographic candidates are E because release and estimand compatibility remain unresolved.",
        "Single-snapshot counts are eligible core/extension rows without A/B/E in that comparator. *The union count is 361 match-absent rows across all four snapshots.",
        "Fan's complete SD11 table was joined by source row reference; no Fan rg/SE/P/FDR values are redistributed. SleepChart is carried forward as audited context only.",
        "D denotes release/table-scope absence. CTG-VL had no direct export. No new GWAS, LDSC estimate, or biological evidence was generated.",
    ]
    for i, line in enumerate(notes):
        parts.append(text(70, y0 + i*28, line, 14, "#48515d"))
    parts.append('</svg>')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
