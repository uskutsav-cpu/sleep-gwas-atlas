#!/usr/bin/env python3
"""Generate an evidence-bound Phase 0/1 Markdown report.

The report only interprets values that are present in its input tables. It
marks smoke-test inputs prominently and never substitutes a canned biological
finding for a missing LDSC result.
"""
import argparse
import os
from pathlib import Path

import pandas as pd


def markdown_value(value):
    if pd.isna(value):
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ")


def dataframe_to_markdown(frame):
    if frame.empty:
        return "_No rows available._"
    headers = [markdown_value(column) for column in frame.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(markdown_value(value) for value in row.values) + " |")
    return "\n".join(lines)


def table_subset(frame, columns):
    return frame[[column for column in columns if column in frame.columns]]


def smoke_test_input(paths):
    return any("_smoketest" in Path(path).parts for path in paths if path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--h2", default="results/tables/h2_summary.tsv")
    parser.add_argument("--rg", default="results/tables/rg_matrix.tsv")
    parser.add_argument("--out", default="results/tables/phase1_summary_report.md")
    args = parser.parse_args()

    paths = [args.h2, args.rg]
    synthetic = smoke_test_input(paths)
    lines = ["# Phase 0/1 report: Sleep/Circadian Genetic Atlas", ""]
    if synthetic:
        lines.extend([
            "> **SYNTHETIC SMOKE-TEST INPUT - NOT A SCIENTIFIC RESULT.**",
            "> This report exists only to validate the software path; do not present or interpret its numbers.",
            "",
        ])
    else:
        lines.extend([
            "**Input status:** Derived from the listed local LDSC output tables. "
            "This report does not independently validate source-GWAS provenance or reproduce LDSC.",
            "",
        ])

    lines.extend(["## 1. Trait registry", ""])
    if os.path.exists(args.config):
        config = pd.read_csv(args.config, sep="\t", dtype=str)
        lines.append(f"- Registered traits: **{len(config)}**")
        lines.append("- Traits by domain:")
        for domain, count in config["domain"].value_counts().sort_index().items():
            lines.append(f"  - `{domain}`: {count}")
        if "source_status" in config:
            lines.append("- Source-verification status by declared value:")
            for status, count in config["source_status"].fillna("MISSING").value_counts().sort_index().items():
                lines.append(f"  - `{status}`: {count}")
    else:
        lines.append(f"_Configuration file not found: `{args.config}`._")
    lines.append("")

    lines.extend(["## 2. SNP heritability and QC gate", ""])
    if os.path.exists(args.h2):
        heritability = pd.read_csv(args.h2, sep="\t")
        lines.append(f"- Traits with parsed LDSC h2: **{len(heritability)}**")
        if "verdict" in heritability:
            counts = heritability["verdict"].fillna("UNPARSED").value_counts()
            lines.append(f"- PASS: **{counts.get('PASS', 0)}**; DROP: **{counts.get('DROP', 0)}**")
        lines.extend(["", "Highest h2 Z-scores (descriptive only):", ""])
        top = heritability.sort_values("z", ascending=False, na_position="last").head(10)
        lines.append(dataframe_to_markdown(table_subset(top, [
            "trait", "scale", "h2", "se", "z", "intercept", "verdict",
            "qc_reason", "ldsc_n_eff_h2", "ldsc_n_eff_h2_gt_12000",
            "mixer_univariate_status",
        ])))
        lines.extend([
            "",
            "Gate definition: h2 Z >= 4 and LDSC intercept <= 1.20. A low h2 Z and an inflated intercept are reported as distinct failure modes.",
        ])
    else:
        lines.append(f"_Heritability table not found: `{args.h2}`._")
    lines.append("")

    lines.extend(["## 3. Genome-wide genetic correlation", ""])
    if os.path.exists(args.rg):
        correlation = pd.read_csv(args.rg, sep="\t")
        lines.append(f"- Parsed sleep-disease pairs: **{len(correlation)}**")
        if "fdr" in correlation:
            lines.append(f"- FDR < 0.05: **{int((correlation['fdr'] < 0.05).sum())}**")
        if "rg" in correlation:
            lines.extend(["", "Strongest positive estimates (descriptive only):", ""])
            lines.append(dataframe_to_markdown(table_subset(correlation.sort_values("rg", ascending=False).head(5), ["sleep_trait", "disease_trait", "rg", "se", "p", "fdr"])))
            lines.extend(["", "Strongest negative estimates (descriptive only):", ""])
            lines.append(dataframe_to_markdown(table_subset(correlation.sort_values("rg", ascending=True).head(5), ["sleep_trait", "disease_trait", "rg", "se", "p", "fdr"])))
    else:
        lines.append(f"_Genetic-correlation table not found: `{args.rg}`._")
    lines.append("")

    lines.extend([
        "## 4. Downstream handoff",
        "",
        "1. Use only FDR-significant, well-powered pairs with source provenance retained for Phase 2 local genetic correlation.",
        "2. Confirm binary-trait prevalence citations before reporting liability-scale h2. Treat LDSC N_eff x h2 only as a screening statistic; actual MiXeR eligibility requires univariate MiXeR diagnostics.",
        "3. Preserve each QC ledger, munging log, LDSC log, and table hash with the figure so Phase 1 results remain auditable.",
        "",
    ])

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print(f"Generated Phase 0/1 report: {args.out}")


if __name__ == "__main__":
    main()
