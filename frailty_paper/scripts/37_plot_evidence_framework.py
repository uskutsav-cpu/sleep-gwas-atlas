#!/usr/bin/env python3
"""Build a source-linked working evidence/design framework figure.

This is deliberately a working figure, not a completed PRISMA diagram. The
review panel shows PubMed acquisition and the unscreened queue only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


SLEEP_ORDER = [
    "insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype",
    "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency",
    "accel_sleep_duration", "sleep_timing",
]
SLEEP_GROUPS = {
    "Self-reported symptoms and behavior": ["insomnia", "sleepiness", "napping", "snoring"],
    "Duration and chronotype": ["sleepdur", "shortsleep", "longsleep", "chronotype"],
    "Clinical record phenotype": ["sleep_apnea"],
    "Actigraphy-derived measures": ["sleep_efficiency", "accel_sleep_duration", "sleep_timing"],
}
SLEEP_LABELS = {
    "insomnia": "Insomnia",
    "sleepdur": "Sleep duration (continuous)",
    "shortsleep": "Short sleep (<7 h)",
    "longsleep": "Long sleep (≥9 h)",
    "chronotype": "Chronotype / morningness",
    "sleepiness": "Daytime sleepiness",
    "napping": "Daytime napping",
    "snoring": "Habitual snoring",
    "sleep_apnea": "Sleep apnea",
    "sleep_efficiency": "Sleep efficiency",
    "accel_sleep_duration": "Sleep duration",
    "sleep_timing": "Sleep timing",
}
INPUTS = [
    "frailty_paper/review/prisma_flow.tsv",
    "frailty_paper/review/source_record_counts.tsv",
    "config/analysis_panel.tsv",
    "frailty_paper/config/analysis_plan_v1.yaml",
    "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv",
    "frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv",
    "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv",
    "frailty_paper/analysis/fried_ffs_source_access_audit_2026-09-23.md",
    "frailty_paper/analysis/physical_component_candidate_source_audit_2026-09-23.md",
    "frailty_paper/manifests/targeted_resource_audit.tsv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def reconcile_review_counts(prisma_rows: list[dict[str, str]],
                            source_rows: list[dict[str, str]]) -> dict[str, int]:
    """Use the current review-flow ledger and source-count ledger, checking they agree."""
    prisma = {row["stage"]: int(row["count"]) for row in prisma_rows}
    sources = {row["database"]: row for row in source_rows}
    if "PubMed" not in sources:
        raise SystemExit("PubMed source-record count is missing")
    counts = {
        "pubmed": prisma["Records identified from PubMed searches"],
        "manual": prisma["Records imported from other databases"],
        "duplicates": prisma["Duplicate records removed across all sources"],
        "queue": prisma["Records after deduplication"],
    }
    source_pubmed = int(sources["PubMed"]["records_imported_or_retrieved"])
    source_queue = int(sources["PubMed"]["unique_records_retained_by_priority"])
    if counts["pubmed"] != source_pubmed or counts["queue"] != source_queue:
        raise SystemExit(f"Review flow and source-record counts disagree: {counts}; "
                         f"source PubMed={source_pubmed}, unique={source_queue}")
    if counts["manual"] != 0:
        raise SystemExit("Manual-database records are present; revise the figure's PubMed-only copy")
    if counts["pubmed"] - counts["duplicates"] != counts["queue"]:
        raise SystemExit(f"PubMed source, deduplication, and queue counts do not reconcile: {counts}")
    not_started_stages = [
        "Records screened", "Records excluded at title/abstract", "Reports sought for retrieval",
        "Reports not retrieved", "Reports assessed for eligibility",
        "Reports excluded after full text", "Studies included",
    ]
    if any(prisma[stage] != 0 for stage in not_started_stages):
        raise SystemExit("Review status changed; revise the figure copy and structure before rerendering")
    return counts


def write_immutable(path: Path, content: bytes, *, replace: bool = False) -> None:
    if path.exists():
        if path.read_bytes() != content and not replace:
            raise SystemExit(f"Refusing to overwrite non-identical output: {path}")
        if path.read_bytes() == content:
            return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def wrap(value: str, limit: int) -> list[str]:
    words = value.split()
    lines: list[str] = []
    line = ""
    for word in words:
        if line and len(line) + 1 + len(word) > limit:
            lines.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        lines.append(line)
    return lines


def esc(value: str) -> str:
    return html.escape(value, quote=True)


class Figure:
    def __init__(self) -> None:
        self.parts = [
            '<svg xmlns="http://www.w3.org/2000/svg" width="2000" height="1500" viewBox="0 0 2000 1500">',
            '<rect width="2000" height="1500" fill="#F5F7FA"/>',
            '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#182230}.title{font-size:32px;font-weight:700}.subtitle{font-size:17px;fill:#475467}.badge{font-size:14px;font-weight:700;letter-spacing:.7px;fill:#8A4B08}.panel{fill:#FFFFFF;stroke:#D9E0E8;stroke-width:1.5}.panel-title{font-size:21px;font-weight:700}.panel-sub{font-size:14px;fill:#667085}.card{fill:#F8FAFC;stroke:#CFD8E3;stroke-width:1.2}.blue{fill:#EAF3F7;stroke:#7AA7B8}.teal{fill:#EAF4F1;stroke:#7BAF9F}.amber{fill:#FFF5E8;stroke:#DBA45B}.gray{fill:#F0F2F5;stroke:#B8C0CA}.small{font-size:14px;fill:#475467}.body{font-size:16px}.card-title{font-size:17px;font-weight:700}.count{font-size:26px;font-weight:700;fill:#145C72}.count-amber{font-size:26px;font-weight:700;fill:#9A5A12}.arrow{stroke:#667085;stroke-width:2;fill:none;marker-end:url(#arrow)}</style>',
            '<defs><marker id="arrow" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto"><path d="M0,0 L9,4.5 L0,9 z" fill="#667085"/></marker></defs>',
        ]

    def rect(self, x: int, y: int, w: int, h: int, cls: str = "card", rx: int = 12) -> None:
        self.parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" class="{cls}"/>')

    def text(self, x: int, y: int, value: str, cls: str = "body", anchor: str = "start") -> None:
        self.parts.append(f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{esc(value)}</text>')

    def lines(self, x: int, y: int, value: str, *, limit: int = 70, cls: str = "body", dy: int = 23) -> None:
        rows = wrap(value, limit)
        spans = "".join(f'<tspan x="{x}" dy="{0 if i == 0 else dy}">{esc(line)}</tspan>' for i, line in enumerate(rows))
        self.parts.append(f'<text x="{x}" y="{y}" class="{cls}">{spans}</text>')

    def panel(self, x: int, y: int, w: int, h: int, letter: str, title: str, subtitle: str = "") -> None:
        self.rect(x, y, w, h, "panel", 16)
        self.text(x + 24, y + 39, letter, "count")
        self.text(x + 65, y + 39, title, "panel-title")
        if subtitle:
            self.text(x + 26, y + 65, subtitle, "panel-sub")

    def arrow(self, x1: int, y1: int, x2: int, y2: int) -> None:
        self.parts.append(f'<path d="M{x1},{y1} L{x2},{y2}" class="arrow"/>')

    def bytes(self) -> bytes:
        self.parts.append("</svg>")
        return ("\n".join(self.parts) + "\n").encode("utf-8")


def build_svg(counts: dict[str, int], sleep_rows: dict[str, dict[str, str]], fi_sig: int, latent_sig: int, aging_sig: int) -> bytes:
    f = Figure()
    f.text(44, 54, "Study framework and current evidence status", "title")
    f.text(44, 84, "Sleep and circadian phenotypes × frailty definitions × aging context", "subtitle")
    f.rect(1530, 34, 425, 36, "amber", 18)
    f.text(1742, 58, "WORKING FIGURE · REVIEW SCREENING NOT STARTED", "badge", "middle")

    # A: Identification and deduplication only. The downstream review has not begun.
    f.panel(40, 120, 940, 440, "A", "Evidence acquisition", "PubMed-only snapshot; manual database exports are still absent")
    f.rect(75, 210, 250, 100, "blue")
    f.text(95, 246, f"{counts['pubmed']:,}", "count")
    f.text(95, 277, "records identified in PubMed", "card-title")
    f.text(95, 298, "311 XML files · 3 saved searches", "small")
    f.arrow(336, 260, 385, 260)
    f.rect(400, 210, 225, 100, "gray")
    f.text(420, 246, f"{counts['duplicates']:,}", "count")
    f.text(420, 277, "duplicates linked", "card-title")
    f.text(420, 298, "audit trail retained", "small")
    f.arrow(636, 260, 685, 260)
    f.rect(700, 210, 230, 100, "teal")
    f.text(720, 246, f"{counts['queue']:,}", "count")
    f.text(720, 277, "unique records queued", "card-title")
    f.text(720, 298, "title / abstract screen pending", "small")
    f.rect(75, 342, 405, 122, "amber")
    f.text(96, 374, f"{counts['manual']:,} records", "count-amber")
    f.text(96, 401, "Embase · Scopus · Web of Science · PsycINFO", "card-title")
    f.lines(96, 428, "Licensed searches/exports not supplied; this is not a complete systematic search.", limit=50, cls="small", dy=19)
    f.rect(505, 342, 425, 122, "gray")
    f.text(528, 376, "0 screened · 0 assessed · 0 included", "card-title")
    f.lines(528, 407, "The current deduplicated PubMed queue is ready for screening. Screening and eligibility counts remain pending.", limit=51, cls="small", dy=19)
    f.text(76, 510, "Counts are a reconciliation of retained XML and deduplicated records, not completed screening outcomes.", "small")

    # B: Exact analysis families and the gated downstream chain.
    f.panel(1000, 120, 960, 440, "B", "Study design and analysis sequence", "Distinct families are kept separate; correlations are not causal effects")
    f.rect(1035, 205, 275, 110, "blue")
    f.text(1055, 236, "PRIMARY · FI", "card-title")
    f.lines(1055, 263, "12 sleep × FI pairs; frozen-atlas reuse; 9/12 inherited all-396 BH q≤0.05.", limit=32, cls="small", dy=19)
    f.rect(1340, 205, 275, 110, "teal")
    f.text(1360, 236, "SECONDARY · LATENT FACTORS", "card-title")
    f.lines(1360, 263, "12 × 7 = 84 pairs; 48/84 q≤0.05; sensitivity-only; exact overlap unknown.", limit=33, cls="small", dy=19)
    f.rect(1645, 205, 280, 110, "gray")
    f.text(1665, 236, "AGING CONTEXT", "card-title")
    f.lines(1665, 263, "72 read-only pairs; 26/72 inherited all-396 q≤0.05; not replication.", limit=33, cls="small", dy=19)
    f.arrow(1172, 326, 1172, 356)
    f.arrow(1478, 326, 1478, 356)
    f.arrow(1785, 326, 1785, 356)
    f.rect(1035, 365, 890, 118, "gray")
    f.text(1058, 398, "Independent replication · local sharing · shared loci · fine-mapping · colocalization · QTL/cell context", "card-title")
    f.lines(1058, 430, "Not completed for frailty. Proceed only after source, overlap, input, and multiplicity gates pass; do not infer mechanisms from global rg.", limit=94, cls="small", dy=20)
    f.text(1038, 518, "UK Biobank overlap is expected in 11/12 primary pairs and in latent-family sources; exact intersections remain unknown.", "small")

    # C: Taxonomy derived from the locked panel labels and definitions.
    f.panel(40, 590, 940, 770, "C", "Locked sleep and circadian panel", "12 registered phenotypes; groups organize measures without treating them as interchangeable")
    group_y = 680
    group_heights = [144, 144, 125, 144]
    for (group_name, trait_ids), h in zip(SLEEP_GROUPS.items(), group_heights):
        f.rect(75, group_y, 870, h, "card")
        f.text(100, group_y + 31, group_name, "card-title")
        names = " · ".join(SLEEP_LABELS[t] for t in trait_ids)
        f.lines(100, group_y + 65, names, limit=78, cls="body", dy=24)
        ids = ", ".join(trait_ids)
        f.text(100, group_y + h - 20, f"Panel IDs: {ids}", "small")
        group_y += h + 13
    f.text(76, 1328, "Panel taxonomy follows the locked atlas-v1.0 source registry.", "small")

    # D: Definitions and current analysis eligibility.
    f.panel(1000, 590, 960, 770, "D", "Frailty definitions and evidence status", "The constructs overlap but are not interchangeable")
    cards: Iterable[tuple[str, str, str]] = [
        ("Frailty Index · primary", "Deficit-accumulation measure; current central endpoint. 12 frozen-atlas sleep comparisons are available.", "blue"),
        ("Latent factors · secondary sensitivity", "General frailty plus six specific factors; based on 30 deficits. UK Biobank cohort membership is reported; exact overlap is unknown.", "teal"),
        ("Fried phenotype · construct replication", "Distinct physical phenotype. Exact full summary statistics are not acquired; no sleep-pair estimate is available.", "amber"),
        ("Hospital Frailty Risk Score", "Continuous weighted 109-code EHR score. Publisher hit tables are available; full summary statistics and the exact custom endpoint file remain unverified.", "amber"),
        ("Five physical frailty components", "Weight loss, exhaustion, low activity, slow walking, low grip. Deposited GWAS build/effect/provenance remain unresolved; not harmonized.", "amber"),
    ]
    y = 680
    for label, desc, cls in cards:
        f.rect(1035, y, 890, 116, cls)
        f.text(1060, y + 31, label, "card-title")
        f.lines(1060, y + 61, desc, limit=98, cls="small", dy=20)
        y += 128
    f.text(1037, 1338, "No eligible physical-component, Fried, or HFRS sleep analysis is represented here.", "small")

    # Shared caveat ribbon.
    f.rect(40, 1390, 1920, 66, "gray", 12)
    f.text(65, 1418, "INTERPRETATION", "card-title")
    f.text(225, 1418, "Genome-wide genetic correlation describes shared effects; it does not establish causation, independent replication, local sharing, or mechanism.", "small")
    f.text(225, 1440, "Figure 1 remains provisional until licensed searches, dual screening, eligibility assessment, and the reporting flow are complete.", "small")
    return f.bytes()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--replace-outputs", action="store_true", help="Replace this figure's existing outputs after review")
    args = parser.parse_args()
    repo = args.repo.resolve()
    outdir = (args.outdir or repo / "frailty_paper/analysis").resolve()

    input_paths = {rel: repo / rel for rel in INPUTS}
    missing = [rel for rel, path in input_paths.items() if not path.is_file()]
    if missing:
        raise SystemExit(f"Missing required frozen inputs: {missing}")

    counts = reconcile_review_counts(
        tsv(input_paths["frailty_paper/review/prisma_flow.tsv"]),
        tsv(input_paths["frailty_paper/review/source_record_counts.tsv"]))

    panel = tsv(input_paths["config/analysis_panel.tsv"])
    sleep_rows = {row["trait_id"]: row for row in panel if row["trait_id"] in SLEEP_ORDER}
    if set(sleep_rows) != set(SLEEP_ORDER) or len(sleep_rows) != 12:
        raise SystemExit("The locked 12-trait sleep panel does not match the figure taxonomy")
    if any(row["panel_version"] != "atlas-v1.0" or row["domain"] != "sleep" for row in sleep_rows.values()):
        raise SystemExit("The figure input does not match the frozen sleep panel version/domain")
    grouped_traits = [trait for trait_ids in SLEEP_GROUPS.values() for trait in trait_ids]
    if sorted(grouped_traits) != sorted(SLEEP_ORDER) or len(grouped_traits) != len(set(grouped_traits)):
        raise SystemExit("The sleep taxonomy must include each locked trait exactly once")
    fi_rows = tsv(input_paths["frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv"])
    latent_rows = tsv(input_paths["frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"])
    aging_rows = tsv(input_paths["frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv"])
    if len(fi_rows) != 12:
        raise SystemExit("Unexpected FI frozen result family")
    if len(latent_rows) != 84 or any(row["family_denominator"] != "84" or row["interpretation_status"] != "SENSITIVITY_ONLY" for row in latent_rows):
        raise SystemExit("Unexpected latent-factor result family or status")
    if len(aging_rows) != 72:
        raise SystemExit("Unexpected aging-context extraction size")
    fi_sig = sum(float(row["global_rg_fdr_all_396"]) <= 0.05 for row in fi_rows)
    latent_sig = sum(float(row["q_value"]) <= 0.05 for row in latent_rows)
    aging_sig = sum(float(row["global_rg_fdr_all_396"]) <= 0.05 for row in aging_rows)
    if (fi_sig, latent_sig, aging_sig) != (9, 48, 26):
        raise SystemExit(f"Frozen significance counts changed: FI={fi_sig}, latent={latent_sig}, aging={aging_sig}")

    input_rows: list[dict[str, str]] = []
    def add(panel_name: str, item: str, label: str, value: str, status: str, source: str, note: str = "") -> None:
        input_rows.append({"panel": panel_name, "item": item, "label": label, "value": value, "status": status, "source": source, "note": note})

    for label, count_key, status, note in [
        ("PubMed records identified", "pubmed", "ACQUIRED", "311 XML files; exact query counts and checksums are retained in review manifests"),
        ("Records from other databases", "manual", "NOT_ACQUIRED", "Licensed database exports are absent"),
        ("Duplicate occurrences linked", "duplicates", "DEDUPLICATED_WITH_AUDIT_TRAIL", "PMID, DOI, normalized-title priority"),
        ("Deduplicated screening queue", "queue", "AWAITING_SCREENING", "Screening has not started"),
    ]:
        add("A_review_flow", count_key, label, str(counts[count_key]), status, "frailty_paper/review/prisma_flow.tsv", note)
    add("A_review_flow", "screened", "Records screened", "0", "NOT_STARTED", "frailty_paper/review/prisma_flow.tsv")
    add("A_review_flow", "full_texts", "Reports assessed for eligibility", "0", "NOT_STARTED", "frailty_paper/review/prisma_flow.tsv", "No reports have been sought or retrieved")
    add("A_review_flow", "included", "Studies included", "0", "NOT_STARTED", "frailty_paper/review/prisma_flow.tsv")
    add("B_design", "primary_fi", "Sleep × Frailty Index", "12 pairs; 9/12 inherited BH q≤0.05", "READ_ONLY_FROZEN_ATLAS_REUSE", "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv", "All-396 q-values retained; not independent replication")
    add("B_design", "latent", "Sleep × latent factors", "84 pairs; 48/84 BH q≤0.05", "SECONDARY_SENSITIVITY_ONLY", "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv", "Exact participant overlap unknown")
    add("B_design", "aging", "Aging context", "72 pairs; 26/72 inherited BH q≤0.05", "READ_ONLY_FROZEN_ATLAS_REUSE", "frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv", "Not independent replication")
    add("B_design", "downstream", "Replication, local, locus, molecular", "Not completed", "GATED_OR_NOT_JUSTIFIED", "frailty_paper/config/analysis_plan_v1.yaml", "Only proceed after source, overlap, input and family gates")
    for group, traits in SLEEP_GROUPS.items():
        for trait in traits:
            row = sleep_rows[trait]
            add("C_sleep_taxonomy", trait, SLEEP_LABELS[trait], group, "LOCKED_PANEL", "config/analysis_panel.tsv", row["phenotype_definition"])
    for item, label, value, status, source, note in [
        ("frailty_index", "Frailty Index", "Deficit accumulation; current central endpoint", "PRIMARY; FI ESTIMATES AVAILABLE", "frailty_paper/config/analysis_plan_v1.yaml", "12 frozen-atlas sleep pairs"),
        ("latent_factors", "Latent frailty factors", "General factor + six specific factors from 30 deficits", "SECONDARY; SENSITIVITY ONLY", "frailty_paper/config/analysis_plan_v1.yaml", "Exact participant overlap unknown"),
        ("fried", "Fried Frailty Score", "Physical frailty phenotype", "SOURCE NOT ACQUIRED", "frailty_paper/analysis/fried_ffs_source_access_audit_2026-09-23.md", "Full statistics unavailable from verified route"),
        ("hfrs", "Hospital Frailty Risk Score", "Continuous weighted 109-code EHR score", "EXACT GWAS FILE UNVERIFIED", "frailty_paper/manifests/targeted_resource_audit.tsv", "No sleep-pair estimate available"),
        ("components", "Five physical components", "Weight loss, exhaustion, low activity, slow walking, low grip", "SOURCE QC BLOCKED", "frailty_paper/analysis/physical_component_candidate_source_audit_2026-09-23.md", "Build/effect/provenance unresolved; not harmonized"),
    ]:
        add("D_frailty_definitions", item, label, value, status, source, note)

    headers = ["panel", "item", "label", "value", "status", "source", "note"]
    source_tsv = "\t".join(headers) + "\n" + "".join(
        "\t".join(row[h].replace("\t", " ").replace("\n", " ") for h in headers).rstrip("\t") + "\n"
        for row in input_rows
    )
    tsv_path = outdir / "figure1_evidence_framework_source_data.tsv"
    write_immutable(tsv_path, source_tsv.encode("utf-8"), replace=args.replace_outputs)

    svg = build_svg(counts, sleep_rows, fi_sig, latent_sig, aging_sig)
    svg_path = outdir / "figure1_evidence_framework.svg"
    write_immutable(svg_path, svg, replace=args.replace_outputs)
    pdf_path = outdir / "figure1_evidence_framework.pdf"
    png_path = outdir / "figure1_evidence_framework.png"
    if args.replace_outputs or not pdf_path.exists():
        subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(pdf_path), str(svg_path)], check=True)
    if args.replace_outputs or not png_path.exists():
        subprocess.run(["rsvg-convert", "-w", "4000", "-f", "png", "-b", "#F5F7FA", "-o", str(png_path), str(svg_path)], check=True)

    script_path = Path(__file__).resolve()
    provenance = {
        "schema_version": "frailty_paper_figure1_evidence_framework.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": [sys.executable, script_path.relative_to(repo).as_posix(), *sys.argv[1:]],
        "working_directory": str(repo),
        "script": script_path.relative_to(repo).as_posix(),
        "script_sha256": sha256(script_path),
        "python": platform.python_version(),
        "rsvg_convert_version": subprocess.check_output(["rsvg-convert", "--version"], text=True).strip(),
        "inputs": {rel: sha256(path) for rel, path in input_paths.items()},
        "source_data_tsv": tsv_path.relative_to(repo).as_posix(),
        "source_data_tsv_sha256": sha256(tsv_path),
        "review_counts": counts,
        "analysis_summary": {"fi_pairs": len(fi_rows), "fi_q_le_0_05_all_396": fi_sig, "latent_pairs": len(latent_rows), "latent_q_le_0_05_family_84": latent_sig, "aging_context_pairs": len(aging_rows), "aging_q_le_0_05_all_396": aging_sig},
        "outputs": {path.relative_to(repo).as_posix(): sha256(path) for path in (svg_path, pdf_path, png_path)},
        "rendering": {"engine": "librsvg rsvg-convert", "png_width_px": 4000, "pdf": "vector SVG rendering"},
        "rsvg_convert_executable": subprocess.check_output(["which", "rsvg-convert"], text=True).strip(),
        "render_commands": [
            ["rsvg-convert", "-f", "pdf", "-o", str(pdf_path), str(svg_path)],
            ["rsvg-convert", "-w", "4000", "-f", "png", "-b", "#F5F7FA", "-o", str(png_path), str(svg_path)],
        ],
        "limitations": [
            "Working framework figure, not a completed PRISMA flow or final manuscript figure.",
            f"PubMed only: {counts['queue']:,} unique records remain unscreened; licensed database exports are absent.",
            "FI and aging panels reuse frozen atlas estimates; latent-factor analyses are sensitivity-only with exact overlap unknown.",
            "No independent replication, physical-component, HFRS, local-locus, fine-mapping, colocalization, or molecular follow-up result is shown.",
            "Genetic correlation does not establish causal direction or mechanism.",
        ],
    }
    prov_path = outdir / "figure1_evidence_framework.provenance.json"
    write_immutable(prov_path, (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8"), replace=args.replace_outputs)
    print(json.dumps({"source_rows": len(input_rows), "counts": counts, "analysis_summary": provenance["analysis_summary"], "outputs": provenance["outputs"], "provenance": str(prov_path)}, indent=2))


if __name__ == "__main__":
    main()
