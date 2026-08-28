#!/usr/bin/env python3
"""Create the extension discovery-screen overview without implying novelty proof."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter, defaultdict
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/tmp/sleep-gwas-extension-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", os.environ["MPLCONFIGDIR"])
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


NOVELTY_COLORS = {
    "HIGH_NOVELTY_PRIORITY": "#7B2CBF",
    "UNDEREXPLORED": "#2A9D8F",
    "PREVIOUSLY_SCREENED": "#E9C46A",
    "HEAVILY_STUDIED": "#E76F51",
    "RELATED_EVIDENCE_ONLY": "#457B9D",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rg", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_rg_primary.tsv"),
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/figures/extension_discovery_screen.png"),
    )
    parser.add_argument("--caption-out", type=Path)
    parser.add_argument("--max-heatmap-traits", type=int, default=30)
    args = parser.parse_args()
    if not args.rg.is_file():
        raise SystemExit(f"ERROR: primary extension rg table is missing: {args.rg}")
    rows = read_tsv(args.rg)
    required = {
        "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
        "novelty_priority", "rg", "p", "extension_fdr",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("ERROR: extension rg table is empty or lacks required columns")
    pairs = {(row["sleep_trait"], row["extension_trait_id"]) for row in rows}
    if len(pairs) != len(rows):
        raise SystemExit("ERROR: duplicate pair in extension rg table")

    sleeps = list(dict.fromkeys(row["sleep_trait"] for row in rows))
    trait_meta: dict[str, dict[str, str]] = {}
    min_fdr: dict[str, float] = {}
    for row in rows:
        trait_id = row["extension_trait_id"]
        trait_meta[trait_id] = row
        min_fdr[trait_id] = min(min_fdr.get(trait_id, 1.0), float(row["extension_fdr"]))
    selected = sorted(
        trait_meta,
        key=lambda trait_id: (
            min_fdr[trait_id], trait_meta[trait_id]["phenotype_domain"],
            trait_meta[trait_id]["phenotype_name"], trait_id,
        ),
    )[: args.max_heatmap_traits]

    matrix = np.full((len(sleeps), len(selected)), np.nan)
    sleep_index = {trait: index for index, trait in enumerate(sleeps)}
    trait_index = {trait: index for index, trait in enumerate(selected)}
    for row in rows:
        trait_id = row["extension_trait_id"]
        if trait_id in trait_index:
            matrix[sleep_index[row["sleep_trait"]], trait_index[trait_id]] = float(row["rg"])

    significant = [row for row in rows if float(row["extension_fdr"]) < 0.05]
    domain_total, domain_sig = Counter(), Counter()
    novelty_total, novelty_sig = Counter(), Counter()
    for row in rows:
        domain_total[row["phenotype_domain"]] += 1
        novelty_total[row["novelty_priority"]] += 1
    for row in significant:
        domain_sig[row["phenotype_domain"]] += 1
        novelty_sig[row["novelty_priority"]] += 1

    plt.style.use("seaborn-v0_8-whitegrid")
    figure = plt.figure(figsize=(18, 13), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(1.45, 1.0), width_ratios=(1.45, 1.0))

    heat = figure.add_subplot(grid[0, :])
    limit = max(0.2, float(np.nanmax(np.abs(matrix))))
    image = heat.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    heat.set_yticks(range(len(sleeps)), labels=[trait.replace("_", " ") for trait in sleeps])
    labels = [trait_meta[trait]["phenotype_name"][:32] for trait in selected]
    heat.set_xticks(range(len(selected)), labels=labels, rotation=62, ha="right", fontsize=8)
    heat.set_title(
        f"A  Genetic-correlation overview: top {len(selected)} extension traits by minimum extension FDR",
        loc="left", fontweight="bold",
    )
    heat.set_ylabel("Locked core sleep trait")
    colorbar = figure.colorbar(image, ax=heat, fraction=0.018, pad=0.01)
    colorbar.set_label("LDSC $r_g$")

    domains = sorted(domain_total, key=lambda domain: (-domain_sig[domain], domain))
    domain_ax = figure.add_subplot(grid[1, 0])
    positions = np.arange(len(domains))
    rates = [domain_sig[domain] / domain_total[domain] for domain in domains]
    domain_ax.bar(positions, rates, color="#31572C")
    domain_ax.set_xticks(positions, labels=[domain.replace("_", " ") for domain in domains], rotation=55, ha="right", fontsize=8)
    domain_ax.set_ylabel("Fraction of tested pairs at extension FDR < 0.05")
    domain_ax.set_title("B  Signal yield by prespecified domain", loc="left", fontweight="bold")
    domain_ax.set_ylim(0, max(0.05, max(rates, default=0) * 1.15))

    scatter = figure.add_subplot(grid[1, 1])
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["novelty_priority"]].append(row)
    for category in sorted(grouped):
        group = grouped[category]
        x = [abs(float(row["rg"])) for row in group]
        y = [-math.log10(max(float(row["extension_fdr"]), 1e-300)) for row in group]
        scatter.scatter(x, y, s=22, alpha=0.65, color=NOVELTY_COLORS.get(category, "#777777"), label=category.replace("_", " ").title())
    scatter.axhline(-math.log10(0.05), color="#222222", linestyle="--", linewidth=1)
    scatter.axvline(0.15, color="#222222", linestyle=":", linewidth=1)
    scatter.set_xlabel("Absolute LDSC $r_g$")
    scatter.set_ylabel("$-\\log_{10}$(extension FDR)")
    scatter.set_title("C  Discovery thresholds and pre-screen category", loc="left", fontweight="bold")
    scatter.legend(frameon=False, fontsize=8, loc="upper right")

    figure.suptitle(
        "Novelty-Enriched Phenome Discovery Extension\n"
        "Exploratory screen; pair-level novelty and independent replication remain separate gates",
        fontsize=16, fontweight="bold",
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=300, bbox_inches="tight")
    if args.out.suffix.lower() == ".png":
        figure.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)

    caption_out = args.caption_out or args.out.with_suffix(".caption.txt")
    threshold_hits = sum(
        float(row["extension_fdr"]) < 0.05 and abs(float(row["rg"])) >= 0.15
        for row in rows
    )
    caption_out.write_text(
        "Exploratory LDSC genetic-correlation screen between the 12 immutable core sleep traits "
        f"and {len(trait_meta)} extension traits that passed the extension rerun-h2 gate "
        f"({len(rows)} primary pairs; {len(significant)} pairs at extension-only FDR < 0.05; "
        f"{threshold_hits} also with |rg| >= 0.15). Panel A displays up to "
        f"{args.max_heatmap_traits} traits ranked only for visualization by their minimum FDR. "
        "Pre-screen novelty categories are not pair-level novelty claims. Strong novelty requires "
        "the separate literature audit and independent replication gates.\n",
        encoding="utf-8",
    )
    print(
        f"EXTENSION_FIGURE_OK pairs={len(rows)} fdr_hits={len(significant)} "
        f"threshold_hits={threshold_hits} out={args.out}"
    )


if __name__ == "__main__":
    main()
