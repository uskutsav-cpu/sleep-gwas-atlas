#!/usr/bin/env python3
"""Plot the complete tested replication family and the 23 replicated effects."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sleep-gwas-extension-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", os.environ["MPLCONFIGDIR"])
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path("discovery_extension")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def label(row: dict[str, str]) -> str:
    phenotype = row["external_phenotype_name"]
    for prefix in ("K21 ", "J44 "):
        if phenotype.startswith(prefix):
            phenotype = phenotype[len(prefix):]
    return f"{row['sleep_trait'].replace('_', ' ')} — {phenotype}"


def main() -> None:
    results = read_tsv(ROOT / "results/replication/replication_results.tsv")
    tested = [
        row for row in results
        if row["replication_class"] in {"REPLICATED", "DIRECTIONALLY_CONCORDANT", "FAILED_REPLICATION"}
    ]
    replicated = [row for row in tested if row["replication_class"] == "REPLICATED"]
    if len(tested) != 41 or len(replicated) != 23:
        raise SystemExit(f"ERROR: expected 41 tested and 23 replicated rows; got {len(tested)} and {len(replicated)}")
    if any(row["direction_concordant"] != "True" for row in tested):
        raise SystemExit("ERROR: tested replication family contains a direction-discordant estimate")

    plt.rcParams.update({
        "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 10,
        "figure.dpi": 160, "savefig.dpi": 300,
    })
    fig = plt.figure(figsize=(13.2, 12.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=(0.9, 1.5))
    left_grid = grid[0, 0].subgridspec(2, 1, height_ratios=(1.0, 1.0))
    ax_scatter = fig.add_subplot(left_grid[0, 0])
    ax_flow = fig.add_subplot(left_grid[1, 0])
    ax_forest = fig.add_subplot(grid[0, 1])

    colors = {
        "REPLICATED": "#1f6f8b",
        "DIRECTIONALLY_CONCORDANT": "#d08c29",
        "FAILED_REPLICATION": "#b23a48",
    }
    for status in ("DIRECTIONALLY_CONCORDANT", "REPLICATED", "FAILED_REPLICATION"):
        rows = [row for row in tested if row["replication_class"] == status]
        if not rows:
            continue
        ax_scatter.scatter(
            [float(row["discovery_rg"]) for row in rows],
            [float(row["replication_rg"]) for row in rows],
            s=42, alpha=0.88, color=colors[status], edgecolor="white", linewidth=0.55,
            label=status.replace("_", " ").title(), zorder=3,
        )
    bounds = [-0.30, 0.62]
    ax_scatter.plot(bounds, bounds, linestyle="--", linewidth=1, color="#666666", zorder=1)
    ax_scatter.axhline(0, color="#aaaaaa", linewidth=0.7)
    ax_scatter.axvline(0, color="#aaaaaa", linewidth=0.7)
    ax_scatter.set_xlim(bounds)
    ax_scatter.set_ylim(bounds)
    ax_scatter.set_aspect("equal", adjustable="box")
    ax_scatter.set_xlabel("Discovery genetic correlation")
    ax_scatter.set_ylabel("Replication genetic correlation")
    ax_scatter.set_title("A. All 41 h²-eligible replication tests", loc="left", fontweight="bold")
    ax_scatter.legend(frameon=False, loc="lower right")
    ax_scatter.text(
        0.02, 0.98, "41/41 directionally concordant\n23/41 Bonferroni replicated",
        transform=ax_scatter.transAxes, va="top", ha="left",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#cccccc"},
    )

    ax_flow.axis("off")
    ax_flow.set_title("C. Locked replication-family disposition", loc="left", fontweight="bold")
    flow = [
        (0.88, "217\npre-result candidates", "#e8f1f2"),
        (0.65, "58\nsource-available", "#d6e9ed"),
        (0.42, "41\nh²/intercept eligible", "#bddde3"),
        (0.19, "23\nBonferroni replicated", "#8fc5d1"),
    ]
    for index, (position, text, color) in enumerate(flow):
        ax_flow.text(
            0.5, position, text, ha="center", va="center", fontsize=10,
            bbox={"boxstyle": "round,pad=0.6", "facecolor": color, "edgecolor": "#4f6d75"},
            transform=ax_flow.transAxes,
        )
        if index < len(flow) - 1:
            ax_flow.annotate(
                "", xy=(0.5, flow[index + 1][0] + 0.075), xytext=(0.5, position - 0.075),
                xycoords=ax_flow.transAxes, textcoords=ax_flow.transAxes,
                arrowprops={"arrowstyle": "-|>", "color": "#4f6d75", "linewidth": 1.2},
            )
    ax_flow.text(0.76, 0.66, "159 unavailable", color="#666666", transform=ax_flow.transAxes)
    ax_flow.text(0.76, 0.43, "17 underpowered", color="#666666", transform=ax_flow.transAxes)
    ax_flow.text(0.76, 0.20, "18 additional\ndirectional", color="#666666", transform=ax_flow.transAxes)
    ax_flow.text(
        0.5, 0.03, "Replicated set: 18 exact + 5 comparable phenotype matches;\n7/23 with heterogeneity P<0.05",
        ha="center", va="bottom", color="#333333", transform=ax_flow.transAxes,
    )

    replicated.sort(key=lambda row: float(row["discovery_rg"]))
    y = list(range(len(replicated)))
    discovery = [float(row["discovery_rg"]) for row in replicated]
    discovery_se = [float(row["discovery_se"]) for row in replicated]
    replication = [float(row["replication_rg"]) for row in replicated]
    replication_se = [float(row["replication_se"]) for row in replicated]
    ax_forest.errorbar(
        discovery, [value + 0.16 for value in y],
        xerr=[1.96 * value for value in discovery_se], fmt="o", markersize=4.5,
        color="#264653", ecolor="#264653", elinewidth=1, capsize=2, label="Discovery rg (95% CI)",
    )
    ax_forest.errorbar(
        replication, [value - 0.16 for value in y],
        xerr=[1.96 * value for value in replication_se], fmt="s", markersize=4.2,
        color="#e76f51", ecolor="#e76f51", elinewidth=1, capsize=2, label="Replication rg (95% CI)",
    )
    for index, row in enumerate(replicated):
        if float(row["heterogeneity_p"]) < 0.05:
            ax_forest.text(0.635, index, "*", va="center", ha="center", color="#8c1c13", fontsize=12)
    ax_forest.axvline(0, color="#555555", linewidth=0.8)
    ax_forest.set_xlim(-0.32, 0.66)
    ax_forest.set_yticks(y)
    ax_forest.set_yticklabels([label(row) for row in replicated], fontsize=7.3)
    ax_forest.set_xlabel("Genetic correlation (95% CI)")
    ax_forest.set_title("B. Independently replicated Tier-B pairs", loc="left", fontweight="bold")
    ax_forest.grid(axis="x", color="#e5e5e5", linewidth=0.65)
    ax_forest.legend(
        handles=[
            Line2D([0], [0], marker="o", color="#264653", label="Discovery rg (95% CI)", markersize=5, linestyle="-"),
            Line2D([0], [0], marker="s", color="#e76f51", label="Replication rg (95% CI)", markersize=5, linestyle="-"),
            Line2D([0], [0], marker="*", color="#8c1c13", label="Heterogeneity P<0.05", markersize=8, linestyle="None"),
        ],
        frameon=False, loc="lower right",
    )

    fig.suptitle(
        "Independent replication of novelty-enriched sleep–phenotype genetic correlations",
        fontsize=13, fontweight="bold",
    )
    out_png = ROOT / "figures/extension_replication_summary.png"
    out_pdf = ROOT / "figures/extension_replication_summary.pdf"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    caption = ROOT / "figures/extension_replication_summary.caption.txt"
    caption.write_text(
        "Independent replication of the locked novelty-enriched extension. Panel A shows all 41 pairs whose "
        "replication source passed prespecified h2/intercept QC; all were directionally concordant, and 23 met "
        "the frozen Bonferroni threshold 0.05/217. Panel B compares discovery and replication LDSC rg estimates "
        "with 95% confidence intervals for those 23 pairs. Asterisks mark nominal discovery-versus-replication "
        "effect heterogeneity P<0.05. Genetic correlation is not evidence of causality or a shared causal variant.\n",
        encoding="utf-8",
    )
    print(f"REPLICATION_FIGURE_OK tested={len(tested)} replicated={len(replicated)} png={out_png}")


if __name__ == "__main__":
    main()
