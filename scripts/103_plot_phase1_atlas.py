#!/usr/bin/env python3
"""Generate reproducible exploratory figures for the Phase-1 deep analysis."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


PRIMARY = "PRIMARY_PHASE1"
SENSITIVITY = "QC_FAILED_SENSITIVITY"
FDR = 0.05
DOMAIN_COLORS = {
    "aging": "#7b3294",
    "cancer": "#c51b7d",
    "cardio": "#d95f02",
    "immune": "#1b9e77",
    "metabolic": "#e6ab02",
    "neuro": "#377eb8",
    "psychiatric": "#7570b3",
}
MODE_COLORS = {
    "SUBJECTIVE_SELF_REPORT": "#4c78a8",
    "OBJECTIVE_ACTIGRAPHY": "#f58518",
    "DISEASE_LIKE_CLINICAL": "#54a24b",
    "CIRCADIAN_TIMING": "#b279a2",
}
RG_CMAP = LinearSegmentedColormap.from_list("rg", ["#2166ac", "#f7f7f7", "#b2182b"])


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def load(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty analysis table: {path}")
    return pd.read_csv(path, sep="\t")


def atomic_figure(fig: plt.Figure, path: Path) -> None:
    buffer = io.BytesIO()
    fig.savefig(
        buffer,
        format="png",
        dpi=320,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "sleep-gwas-atlas scripts/103_plot_phase1_atlas.py"},
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(buffer.getvalue())
    temporary.replace(path)
    plt.close(fig)


def configure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.dpi": 120,
    })


def orders(master: pd.DataFrame) -> tuple[list[str], list[str], dict[str, str], dict[str, str]]:
    sleep = master["sleep_trait"].drop_duplicates().tolist()
    external = master["external_trait"].drop_duplicates().tolist()
    sleep_labels = master.drop_duplicates("sleep_trait").set_index("sleep_trait")["sleep_trait_label"].to_dict()
    external_labels = master.drop_duplicates("external_trait").set_index("external_trait")["external_trait_label"].to_dict()
    if len(sleep) != 12 or len(external) != 33:
        fail(f"expected 12x33 matrix, observed {len(sleep)}x{len(external)}")
    return sleep, external, sleep_labels, external_labels


def heatmap(
    values: np.ndarray,
    sleep: list[str],
    external: list[str],
    sleep_labels: dict[str, str],
    external_labels: dict[str, str],
    title: str,
    path: Path,
    stars: np.ndarray | None = None,
    masked: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(14.5, 6.7))
    limit = max(0.65, float(np.nanmax(np.abs(values))))
    shown = np.ma.masked_invalid(values) if masked else values
    cmap = RG_CMAP.copy()
    cmap.set_bad("#ececec")
    image = ax.imshow(shown, aspect="auto", cmap=cmap, vmin=-limit, vmax=limit)
    ax.set_xticks(range(len(external)), [external_labels[x] + (" †" if x in {"t2d", "melanoma"} else "") for x in external], rotation=62, ha="right")
    ax.set_yticks(range(len(sleep)), [sleep_labels[x] for x in sleep])
    for sensitivity_trait in ("t2d", "melanoma"):
        column = external.index(sensitivity_trait)
        ax.add_patch(Rectangle(
            (column - .5, -.5), 1, len(sleep), fill=False, edgecolor="#252525",
            linewidth=1.0, linestyle=(0, (3, 2)), clip_on=False,
        ))
    if stars is not None:
        for y, x in zip(*np.where(stars)):
            ax.text(x, y, "•", ha="center", va="center", fontsize=5.5, color="#111111")
    ax.set_title(title, loc="left", weight="bold")
    ax.set_xlabel("† QC-failed external trait shown only as sensitivity; excluded from primary summaries")
    bar = fig.colorbar(image, ax=ax, fraction=.022, pad=.012)
    bar.set_label("Genetic correlation (rg)")
    atomic_figure(fig, path)


def figure_heatmaps(master: pd.DataFrame, out: Path) -> list[Path]:
    sleep, external, sleep_labels, external_labels = orders(master)
    matrix = master.pivot(index="sleep_trait", columns="external_trait", values="rg").loc[sleep, external]
    sig = master.assign(sig=(master["primary_or_sensitivity"].eq(PRIMARY) & master["fdr"].le(FDR))).pivot(
        index="sleep_trait", columns="external_trait", values="sig"
    ).loc[sleep, external]
    first = out / "phase1_01_full_rg_heatmap.png"
    heatmap(matrix.to_numpy(float), sleep, external, sleep_labels, external_labels,
            "Full Phase-1 sleep–external-trait genetic-correlation matrix", first, stars=sig.to_numpy(bool))
    values = matrix.where(sig).to_numpy(float)
    second = out / "phase1_02_primary_fdr_heatmap.png"
    heatmap(values, sleep, external, sleep_labels, external_labels,
            "Locked-family primary discoveries (BH FDR ≤ 0.05)", second, stars=np.isfinite(values), masked=True)
    return [first, second]


def figure_sleep_hubs(hubs: pd.DataFrame, out: Path) -> Path:
    frame = hubs.sort_values(["fdr_significant_connections", "sum_abs_rg_significant"])
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    colors = frame["measurement_category"].map(MODE_COLORS)
    ax.barh(frame["sleep_trait_label"], frame["fdr_significant_connections"], color=colors)
    for y, value in enumerate(frame["fdr_significant_connections"]):
        ax.text(value + .25, y, str(value), va="center")
    ax.set_xlim(0, max(frame["fdr_significant_connections"]) + 3)
    ax.set_xlabel("Primary external traits connected at locked BH FDR ≤ 0.05 (of 31)")
    ax.set_title("Sleep-trait hub breadth", loc="left", weight="bold")
    handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=color, markeredgecolor="none", label=key.replace("_", " ").title()) for key, color in MODE_COLORS.items()]
    ax.legend(handles=handles, fontsize=7, loc="lower right")
    path = out / "phase1_03_sleep_trait_hubs.png"
    atomic_figure(fig, path)
    return path


def figure_external_hubs(hubs: pd.DataFrame, out: Path) -> Path:
    frame = hubs.sort_values(["significant_sleep_connections", "maximum_abs_rg"]).tail(20)
    fig, ax = plt.subplots(figsize=(8.8, 6.5))
    colors = frame["external_domain"].map(DOMAIN_COLORS)
    ax.barh(frame["external_trait_label"], frame["significant_sleep_connections"], color=colors)
    ax.set_xlabel("Sleep traits connected at locked BH FDR ≤ 0.05 (of 12)")
    ax.set_title("Broadest external-trait sleep hubs", loc="left", weight="bold")
    handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=color, markeredgecolor="none", label=key.title()) for key, color in DOMAIN_COLORS.items()]
    ax.legend(handles=handles, fontsize=7, ncol=2, loc="lower right")
    path = out / "phase1_04_external_trait_hubs.png"
    atomic_figure(fig, path)
    return path


def figure_domain_matrix(master: pd.DataFrame, out: Path) -> Path:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)].copy()
    primary["sig"] = primary["fdr"].le(FDR)
    matrix = primary.pivot_table(index="sleep_trait_label", columns="external_domain", values="sig", aggfunc="sum")
    sleep_order = primary["sleep_trait_label"].drop_duplicates().tolist()
    domains = sorted(primary["external_domain"].unique())
    matrix = matrix.loc[sleep_order, domains]
    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    image = ax.imshow(matrix.to_numpy(float), cmap="YlGnBu", aspect="auto", vmin=0, vmax=matrix.to_numpy().max())
    ax.set_xticks(range(len(domains)), [x.title() for x in domains], rotation=35, ha="right")
    ax.set_yticks(range(len(sleep_order)), sleep_order)
    for y in range(len(sleep_order)):
        for x in range(len(domains)):
            ax.text(x, y, int(matrix.iloc[y, x]), ha="center", va="center", fontsize=7,
                    color="white" if matrix.iloc[y, x] >= matrix.to_numpy().max() * .55 else "#111111")
    ax.set_title("Locked-FDR connections by sleep trait and clinical domain", loc="left", weight="bold")
    fig.colorbar(image, ax=ax, fraction=.035, pad=.02, label="Significant pair count")
    path = out / "phase1_05_domain_connectivity_matrix.png"
    atomic_figure(fig, path)
    return path


def similarity_matrix(similarity: pd.DataFrame, order: list[str]) -> np.ndarray:
    matrix = np.eye(len(order))
    index = {name: i for i, name in enumerate(order)}
    for row in similarity.itertuples(index=False):
        a, b = index[row.sleep_trait_1], index[row.sleep_trait_2]
        matrix[a, b] = matrix[b, a] = row.pearson_rg_profile
    return matrix


def figure_similarity(master: pd.DataFrame, similarity: pd.DataFrame, dendro_order: pd.DataFrame, out: Path) -> Path:
    order = dendro_order.sort_values("dendrogram_order")["sleep_trait"].tolist()
    labels = master.drop_duplicates("sleep_trait").set_index("sleep_trait")["sleep_trait_label"].to_dict()
    matrix = similarity_matrix(similarity, order)
    fig, ax = plt.subplots(figsize=(7, 6.2))
    image = ax.imshow(matrix, cmap=RG_CMAP, vmin=-1, vmax=1)
    ax.set_xticks(range(len(order)), [labels[x] for x in order], rotation=60, ha="right")
    ax.set_yticks(range(len(order)), [labels[x] for x in order])
    ax.set_title("Similarity of 31-trait genetic-correlation profiles", loc="left", weight="bold")
    fig.colorbar(image, ax=ax, fraction=.035, pad=.02, label="Pearson correlation between rg profiles")
    path = out / "phase1_06_sleep_profile_similarity.png"
    atomic_figure(fig, path)
    return path


def figure_dendrogram(master: pd.DataFrame, linkage: pd.DataFrame, dendro_order: pd.DataFrame, out: Path) -> Path:
    sleep = master["sleep_trait"].drop_duplicates().tolist()
    labels = master.drop_duplicates("sleep_trait").set_index("sleep_trait")["sleep_trait_label"].to_dict()
    leaf_order = dendro_order.sort_values("dendrogram_order")["sleep_trait"].tolist()
    x = {sleep.index(name): float(i) for i, name in enumerate(leaf_order)}
    height = {i: 0.0 for i in range(len(sleep))}
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    for step, row in enumerate(linkage.itertuples(index=False)):
        left, right = int(row.left_node), int(row.right_node)
        parent = len(sleep) + step
        lx, rx = x[left], x[right]
        lh, rh = height[left], height[right]
        top = float(row.distance)
        ax.plot([lx, lx, rx, rx], [lh, top, top, rh], color="#333333", lw=1.2)
        x[parent] = (lx + rx) / 2
        height[parent] = top
    ax.set_xticks(range(len(leaf_order)), [labels[name] for name in leaf_order], rotation=55, ha="right")
    ax.set_ylabel("Average-linkage distance (1 − Pearson r)")
    ax.set_title("Descriptive clustering of sleep rg profiles", loc="left", weight="bold")
    ax.text(.99, .98, "Not a latent-factor model", transform=ax.transAxes, ha="right", va="top", color="#666666")
    path = out / "phase1_07_sleep_profile_dendrogram.png"
    atomic_figure(fig, path)
    return path


def figure_pca(master: pd.DataFrame, pca: pd.DataFrame, out: Path) -> Path:
    modes = master.drop_duplicates("sleep_trait").set_index("sleep_trait")["sleep_measurement_category"].to_dict()
    labels = master.drop_duplicates("sleep_trait").set_index("sleep_trait")["sleep_trait_label"].to_dict()
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for row in pca.itertuples(index=False):
        color = MODE_COLORS[modes[row.sleep_trait]]
        ax.scatter(row.PC1, row.PC2, s=52, color=color, edgecolor="white", linewidth=.7, zorder=3)
        ax.annotate(labels[row.sleep_trait], (row.PC1, row.PC2), xytext=(4, 4), textcoords="offset points", fontsize=7)
    variance1 = pca.iloc[0]["PC1_variance_explained"] * 100
    variance2 = pca.iloc[0]["PC2_variance_explained"] * 100
    ax.axhline(0, color="#dddddd", lw=.8); ax.axvline(0, color="#dddddd", lw=.8)
    ax.set_xlabel(f"PC1 ({variance1:.1f}% descriptive profile variance)")
    ax.set_ylabel(f"PC2 ({variance2:.1f}% descriptive profile variance)")
    ax.set_title("Two-dimensional view of sleep rg profiles", loc="left", weight="bold")
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none", label=key.replace("_", " ").title()) for key, color in MODE_COLORS.items()]
    ax.legend(handles=handles, fontsize=7, loc="best")
    path = out / "phase1_08_sleep_profile_pca.png"
    atomic_figure(fig, path)
    return path


def network_positions(master: pd.DataFrame) -> tuple[dict[str, tuple[float, float]], list[str], list[str]]:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)]
    sleep = primary["sleep_trait"].drop_duplicates().tolist()
    external_frame = primary.drop_duplicates("external_trait").sort_values(["external_domain", "external_trait_label"])
    external = external_frame["external_trait"].tolist()
    positions: dict[str, tuple[float, float]] = {}
    for i, trait in enumerate(sleep):
        positions[f"sleep::{trait}"] = (0.0, 1 - i / max(1, len(sleep) - 1))
    for i, trait in enumerate(external):
        positions[f"external::{trait}"] = (1.0, 1 - i / max(1, len(external) - 1))
    return positions, sleep, external


def figure_network(master: pd.DataFrame, out: Path, negative_only: bool) -> Path:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY) & master["fdr"].le(FDR)].copy()
    if negative_only:
        primary = primary.loc[primary["rg"].lt(0)]
    positions, sleep, external = network_positions(master)
    sleep_labels = master.drop_duplicates("sleep_trait").set_index("sleep_trait")["sleep_trait_label"].to_dict()
    external_meta = master.drop_duplicates("external_trait").set_index("external_trait")[["external_trait_label", "external_domain"]].to_dict("index")
    fig, ax = plt.subplots(figsize=(11, 9))
    for row in primary.itertuples(index=False):
        start = positions[f"sleep::{row.sleep_trait}"]
        end = positions[f"external::{row.external_trait}"]
        color = "#b2182b" if row.rg > 0 else "#2166ac"
        ax.plot([start[0], end[0]], [start[1], end[1]], color=color,
                alpha=.12 + .48 * min(row.abs_rg / .65, 1), lw=.35 + 2.0 * row.abs_rg / .65, zorder=1)
    for trait in sleep:
        x, y = positions[f"sleep::{trait}"]
        ax.scatter(x, y, s=33, color=MODE_COLORS[master.loc[master["sleep_trait"].eq(trait), "sleep_measurement_category"].iloc[0]], zorder=3)
        ax.text(x - .018, y, sleep_labels[trait], ha="right", va="center", fontsize=7)
    connected = set(primary["external_trait"])
    for trait in external:
        x, y = positions[f"external::{trait}"]
        meta = external_meta[trait]
        alpha = 1 if trait in connected else .2
        ax.scatter(x, y, s=27, color=DOMAIN_COLORS[meta["external_domain"]], alpha=alpha, zorder=3)
        ax.text(x + .018, y, meta["external_trait_label"], ha="left", va="center", fontsize=6.5, alpha=alpha)
    ax.set_xlim(-.35, 1.38); ax.set_ylim(-.03, 1.03); ax.axis("off")
    ax.set_title(("Significant negative" if negative_only else "Signed primary significant") + " sleep–trait network", loc="left", weight="bold")
    ax.text(.5, -.02, f"All {len(primary)} locked-FDR edges shown; width scales with |rg|", transform=ax.transAxes, ha="center", color="#555555")
    path = out / ("phase1_10_negative_rg_network.png" if negative_only else "phase1_09_bipartite_network.png")
    atomic_figure(fig, path)
    return path


def figure_measurement(measurement: pd.DataFrame, out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for domain, frame in measurement.groupby("external_domain"):
        ax.scatter(frame["subjective_mean_rg"], frame["objective_mean_rg"], s=42,
                   color=DOMAIN_COLORS[domain], label=domain.title(), alpha=.88, edgecolor="white", linewidth=.5)
    limit = max(abs(measurement["subjective_mean_rg"]).max(), abs(measurement["objective_mean_rg"]).max()) + .03
    ax.plot([-limit, limit], [-limit, limit], color="#777777", ls="--", lw=.8)
    ax.axhline(0, color="#dddddd", lw=.8); ax.axvline(0, color="#dddddd", lw=.8)
    notable = measurement.loc[
        measurement["objective_materially_larger"]
        | measurement["subjective_materially_larger"]
    ].sort_values("objective_minus_subjective_mean_abs_rg")
    for row in notable.itertuples(index=False):
        ax.annotate(row.external_trait_label, (row.subjective_mean_rg, row.objective_mean_rg), xytext=(3, 3), textcoords="offset points", fontsize=6.5)
    ax.set_xlim(-limit, limit); ax.set_ylim(-limit, limit)
    ax.set_xlabel("Mean rg across 7 self-reported sleep traits")
    ax.set_ylabel("Mean rg across 3 actigraphy sleep traits")
    ax.set_title("Self-report versus actigraphy rg profiles", loc="left", weight="bold")
    ax.legend(ncol=2, fontsize=7, loc="best")
    path = out / "phase1_11_objective_subjective_comparison.png"
    atomic_figure(fig, path)
    return path


def figure_known_novel(audit: pd.DataFrame, out: Path) -> Path:
    required = {"sleep_trait", "external_trait", "rg", "abs_rg", "fdr", "novelty_classification"}
    if not required.issubset(audit.columns) or len(audit) != 153:
        fail("literature audit must contain the 153 primary discoveries and required plotting columns")
    direct = {"DIRECT_RG_PREVIOUSLY_REPORTED", "DIRECT_RG_REPLICATION_DIFFERENT_DATASET"}
    audit = audit.copy()
    audit["plot_class"] = np.where(audit["novelty_classification"].isin(direct), "Direct prior rg", "No direct prior rg classified")
    colors = {"Direct prior rg": "#4c78a8", "No direct prior rg classified": "#e45756"}
    fig, ax = plt.subplots(figsize=(7.4, 5.8))
    for group, frame in audit.groupby("plot_class"):
        ax.scatter(frame["abs_rg"], -np.log10(frame["fdr"].clip(lower=np.finfo(float).tiny)),
                   s=30, alpha=.72, color=colors[group], label=f"{group} (n={len(frame)})")
    ax.set_xlabel("Absolute genetic correlation |rg|")
    ax.set_ylabel("−log10 locked-family BH FDR")
    ax.set_title("Effect magnitude and evidence by direct-rg literature status", loc="left", weight="bold")
    ax.legend()
    path = out / "phase1_12_known_vs_novel_effect_size.png"
    atomic_figure(fig, path)
    return path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--quantitative-only", action="store_true", help="Build figures 1–11 before the literature audit is complete")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    analysis = root / "results/analysis"
    output = root / "results/figures/phase1_deep_analysis"
    configure()
    master = load(analysis / "phase1_master_analysis.tsv")
    hubs = load(analysis / "sleep_trait_hubs.tsv")
    external = load(analysis / "external_trait_hubs.tsv")
    similarity = load(analysis / "sleep_profile_similarity.tsv")
    linkage = load(analysis / "sleep_profile_linkage.tsv")
    dendro_order = load(analysis / "sleep_profile_dendrogram_order.tsv")
    pca = load(analysis / "sleep_profile_pca.tsv")
    measurement = load(analysis / "measurement_mode_comparison.tsv")
    generated = figure_heatmaps(master, output)
    generated += [
        figure_sleep_hubs(hubs, output),
        figure_external_hubs(external, output),
        figure_domain_matrix(master, output),
        figure_similarity(master, similarity, dendro_order, output),
        figure_dendrogram(master, linkage, dendro_order, output),
        figure_pca(master, pca, output),
        figure_network(master, output, negative_only=False),
        figure_network(master, output, negative_only=True),
        figure_measurement(measurement, output),
    ]
    if not args.quantitative_only:
        generated.append(figure_known_novel(load(analysis / "literature_novelty_audit.tsv"), output))
    manifest = {
        "schema_version": "sleep-atlas-phase1-figure-manifest.1",
        "generator": "scripts/103_plot_phase1_atlas.py",
        "literature_figure_included": not args.quantitative_only,
        "figures": [
            {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": file_sha256(path)}
            for path in generated
        ],
    }
    manifest_path = output / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PHASE1_FIGURES_BUILT count={len(generated)} output={output.relative_to(root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
