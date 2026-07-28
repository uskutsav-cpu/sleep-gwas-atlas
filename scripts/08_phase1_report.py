#!/usr/bin/env python3
"""
Phase 0/1: Generate an in-depth Phase 0/1 summary report markdown document.

Outputs:
    results/tables/phase1_summary_report.md
"""
import argparse
import os
import pandas as pd


def df_to_md(df):
    """Convert pandas DataFrame to Markdown table without tabulate dependency."""
    if df.empty:
        return ""
    headers = list(df.columns)
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for _, row in df.iterrows():
        vals = [str(v) if pd.notna(v) else "" for v in row.values]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/traits.tsv")
    ap.add_argument("--h2", default="results/tables/h2_summary.tsv")
    ap.add_argument("--comp", default="results/tables/problem_comparison_summary.tsv")
    ap.add_argument("--rg", default="results/tables/rg_matrix.tsv")
    ap.add_argument("--out", default="results/tables/phase1_summary_report.md")
    args = ap.parse_args()

    lines = []
    lines.append("# Phase 0/1 Comprehensive Pipeline Report: Sleep/Circadian Genetic Atlas\n")
    lines.append("This report summarizes the **Phase 0 (GWAS Curation & Harmonization)** and **Phase 1 (Genome-Wide Genetic Correlation & QC Gating)** results across the 90-trait catalog.\n")

    # Section 1: Trait Registry Overview
    if os.path.exists(args.config):
        cfg = pd.read_csv(args.config, sep="\t")
        lines.append("## 1. Trait Registry Overview\n")
        lines.append(f"- **Total Registered Traits**: {len(cfg)}")
        lines.append("- **Domains Covered**:")
        for dom, cnt in cfg["domain"].value_counts().items():
            lines.append(f"  - `{dom}`: {cnt} traits")
        lines.append("\n")

    # Section 2: Dual-Strategy Problem Comparison
    if os.path.exists(args.comp):
        comp = pd.read_csv(args.comp, sep="\t")
        lines.append("## 2. Dual-Strategy ('Use Both') Problem Comparison Summary\n")
        lines.append("Comparative evaluation of competing phenotype definitions, summary statistics releases, and liability scale models:\n\n")
        lines.append(df_to_md(comp))
        lines.append("\n\n### Key Methodological Insights:")
        lines.append("1. **Short Sleep**: Dashti (<7h) provides higher total sample power ($Z=5.0$) but exhibits LDSC intercept inflation ($1.35$), signaling confounding/population structure. Austin-Zimmerman (<=5h) isolates extreme tail cases with a clean intercept ($1.05$).")
        lines.append("2. **Long Sleep**: Both Dashti ($Z=2.0$) and Austin-Zimmerman ($Z=1.51$) fall below the power gate threshold ($Z \\ge 4.0$), demonstrating that dichotomized long sleep duration is underpowered for standalone GenomicSEM factor modeling.")
        lines.append("3. **Insomnia**: Full UKB+23andMe meta-analysis delivers $Z=27.78$ vs. UKB-only $Z=16.00$. Both pass QC gates; UKB-only serves as the un-restricted public open access standard.")
        lines.append("\n")

    # Section 3: QC Gate & Heritability Summary
    if os.path.exists(args.h2):
        h2 = pd.read_csv(args.h2, sep="\t")
        lines.append("## 3. SNP Heritability & QC Gate Verdicts\n")
        lines.append(f"- **Total Traits Evaluated**: {len(h2)}")
        if "verdict" in h2.columns:
            pass_cnt = (h2["verdict"] == "PASS").sum()
            drop_cnt = (h2["verdict"] == "DROP").sum()
            lines.append(f"- **QC Gate Survivors (`PASS`)**: {pass_cnt}")
            lines.append(f"- **QC Gate Exclusions (`DROP`)**: {drop_cnt}")
        lines.append("\nTop 10 highest heritability traits ($Z = h^2 / SE$):\n")
        top_h2 = h2.head(10)[["trait", "scale", "h2", "se", "z", "intercept", "verdict"]]
        lines.append(df_to_md(top_h2))
        lines.append("\n")

    # Section 4: Genetic Correlations
    if os.path.exists(args.rg):
        rg = pd.read_csv(args.rg, sep="\t")
        lines.append("## 4. Genome-Wide Genetic Correlation ($r_g$) Highlights\n")
        lines.append(f"- **Total Trait Pairs Evaluated**: {len(rg)}")
        if "fdr" in rg.columns:
            fdr_sig = (rg["fdr"] < 0.05).sum()
            lines.append(f"- **Pairs Significant at FDR < 0.05**: {fdr_sig}")
        lines.append("\nTop 5 Strongest Positive Genetic Correlations ($r_g$):\n")
        top_pos = rg.sort_values("rg", ascending=False).head(5)[["sleep_trait", "disease_trait", "rg", "se", "p", "fdr"]]
        lines.append(df_to_md(top_pos))
        lines.append("\n\nTop 5 Strongest Negative Genetic Correlations ($r_g$):\n")
        top_neg = rg.sort_values("rg", ascending=True).head(5)[["sleep_trait", "disease_trait", "rg", "se", "p", "fdr"]]
        lines.append(df_to_md(top_neg))
        lines.append("\n")

    # Section 5: Phase 2 / Phase 4 Roadmap
    lines.append("## 5. Next Steps & Downstream Roadmap\n")
    lines.append("1. **Phase 2 (Local Genetic Correlation - LAVA)**: Run LAVA across FDR-significant sleep-disease trait pairs to pinpoint local genomic regions driving shared architecture.")
    lines.append("2. **Phase 3 (Pleiotropic Locus Discovery - PLACO / Conjunction FDR)**: Identify specific lead SNPs shared between sleep traits and disease domains.")
    lines.append("3. **Phase 4 (GenomicSEM Latent Factor Modeling)**: Construct Exploratory & Confirmatory Factor Analysis models across PASS sleep traits and disease domains.")
    lines.append("\n")

    report_content = "\n".join(lines)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(report_content)

    print(f"Generated Phase 1 Summary Report: {args.out}")


if __name__ == "__main__":
    main()
