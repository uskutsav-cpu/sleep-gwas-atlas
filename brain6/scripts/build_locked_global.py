#!/usr/bin/env python3
"""Build the Brain6 global map from the immutable atlas tables.

This is a descriptive extraction only. It never recomputes atlas P values or
multiple-testing corrections and refuses to replace any existing output.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import argparse
import platform
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from figure_style import best_contrast_text


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6"
BRAIN = ["adhd", "mdd", "scz", "bipolar", "alz", "parkinson"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_new(path: Path, data: str | bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    if isinstance(data, bytes):
        path.write_bytes(data)
    else:
        path.write_text(data, encoding="utf-8")


def derive_sleep_traits(panel: pd.DataFrame) -> list[str]:
    """Return ordered sleep trait IDs from the locked atlas panel."""
    if not {"domain", "trait_id"}.issubset(panel.columns):
        raise ValueError("analysis panel must include domain and trait_id")
    traits = panel.loc[panel.domain.eq("sleep"), "trait_id"].astype(str).tolist()
    if not traits or len(traits) != len(set(traits)) or any(not trait for trait in traits):
        raise ValueError("locked sleep panel must contain unique, nonempty trait IDs")
    return traits


def save_figure(fig, stem: str, allow_refresh: bool = False) -> None:
    for ext in ("png", "pdf"):
        path = OUT / "results/global/figures" / f"{stem}.{ext}"
        if path.exists() and not allow_refresh:
            raise FileExistsError(f"Refusing to overwrite existing figure: {path}")
        fig.savefig(path, dpi=320 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def record_figure_provenance(refresh: bool) -> Path:
    figure_dir = OUT / "results/global/figures"
    outputs = sorted([*figure_dir.glob("fig1*.png"), *figure_dir.glob("fig1*.pdf")])
    if len(outputs) != 10 or any(path.stat().st_size == 0 for path in outputs):
        raise ValueError(f"Expected ten nonempty global-map figures; found {len(outputs)}")
    source_paths = [
        ROOT / "config/analysis_panel.tsv",
        ROOT / "results/atlas/traits.tsv",
        ROOT / "results/atlas/trait_pairs.tsv",
        ROOT / "results/analysis/phase1_master_analysis.tsv",
        OUT / "results/global/brain6_72_locked.tsv",
        OUT / "results/global/brain6_72_figure_source.tsv",
        OUT / "results/global/disorder_profile_similarity.tsv",
    ]
    font_path = Path(font_manager.findfont("DejaVu Sans", fallback_to_default=True))
    provenance = {
        "schema_version": 1,
        "status": "PASS",
        "script_sha256": sha256(Path(__file__)),
        "style_helper_sha256": sha256(Path(__file__).with_name("figure_style.py")),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name)
                     for name in ("numpy", "pandas", "scipy", "matplotlib")},
        "font": {"family": "DejaVu Sans", "resolved_path": str(font_path),
                 "sha256": sha256(font_path)},
        "sources": {str(path.relative_to(ROOT)): sha256(path) for path in source_paths},
        "outputs": {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                    for path in outputs},
        "render": {"png_dpi": 320, "pdf_vector": True, "refresh_mode": bool(refresh)},
    }
    destination = figure_dir / "figure_build_provenance.json"
    payload = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    if refresh:
        temporary = destination.with_name(f".{destination.name}.partial")
        temporary.write_text(payload, encoding="utf-8")
        temporary.replace(destination)
    else:
        write_new(destination, payload)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-figures", action="store_true",
                        help="re-render figures only after verifying the existing table matches the locked atlas")
    args = parser.parse_args()
    atlas = pd.read_csv(ROOT / "results/atlas/trait_pairs.tsv", sep="\t")
    traits = pd.read_csv(ROOT / "results/atlas/traits.tsv", sep="\t")
    panel = pd.read_csv(ROOT / "config/analysis_panel.tsv", sep="\t")
    detail = pd.read_csv(ROOT / "results/analysis/phase1_master_analysis.tsv", sep="\t")

    sleep_traits = derive_sleep_traits(panel)
    n_sleep = len(sleep_traits)
    n_brain = len(BRAIN)
    all_brain = traits.loc[traits.trait_id.isin(BRAIN)]
    if len(all_brain) != len(BRAIN):
        raise ValueError("One or more Brain6 disorder IDs are missing from results/atlas/traits.tsv")
    if len(atlas) != 396 or atlas.pair_id.duplicated().any():
        raise ValueError("The source atlas is not the unique 396-pair locked family")

    data = atlas.loc[atlas.non_sleep_trait.isin(BRAIN)].copy()
    expected = pd.MultiIndex.from_product([sleep_traits, BRAIN], names=["sleep_trait", "non_sleep_trait"])
    got = pd.MultiIndex.from_frame(data[["sleep_trait", "non_sleep_trait"]])
    expected_rows = n_sleep * n_brain
    if len(data) != expected_rows or got.has_duplicates or set(got) != set(expected):
        raise ValueError(f"Expected exact {n_sleep}x{n_brain} map; found {len(data)} rows")

    keep = [
        "sleep_trait", "external_trait", "rg", "se", "z", "p", "fdr",
        "locked_primary_significant", "sleep_h2", "sleep_h2_se", "external_h2",
        "external_h2_se", "sleep_h2_intercept", "external_h2_intercept",
        "gcov_intercept", "gcov_intercept_se", "SNP_overlap",
        "interpretation_status", "QC_notes", "phenotype_source", "external_source_id",
        "external_dataset_version", "ancestry", "sample_size", "sleep_sample_size",
        "sleep_measurement_category", "sleep_collection_mode", "source_substitution",
        "source_substitution_caveat", "input_log",
    ]
    data = detail.loc[detail.external_trait.isin(BRAIN), keep].copy()
    if len(data) != expected_rows or data.duplicated(["sleep_trait", "external_trait"]).any():
        raise ValueError("Detailed global-map source table is incomplete or duplicated")
    if not np.allclose(
        data.fdr.to_numpy(),
        atlas.set_index(["sleep_trait", "non_sleep_trait"]).loc[
            pd.MultiIndex.from_frame(data[["sleep_trait", "external_trait"]]).set_names(
                ["sleep_trait", "non_sleep_trait"]
            ), "global_rg_fdr_all_396"
        ].to_numpy(),
        rtol=1e-8, atol=1e-300, equal_nan=True,
    ):
        raise ValueError("The detailed table FDR differs from the locked 396-family FDR")
    tier_lookup = atlas.set_index(["sleep_trait", "non_sleep_trait"])["analysis_tier"]
    data["analysis_tier"] = [tier_lookup.loc[(s, d)] for s, d in zip(data.sleep_trait, data.external_trait)]

    trait_lookup = traits.set_index("trait_id")
    panel_lookup = panel.set_index("trait_id")
    data["brain_disorder"] = data.external_trait
    data["brain_disorder_label"] = data.external_trait.map(trait_lookup.label)
    data["original_BH_FDR_q"] = data.fdr
    data["significance_under_original_396_family"] = data.locked_primary_significant.astype(bool)
    data["h2_disorder"] = data.external_h2
    data["h2_disorder_se"] = data.external_h2_se
    data["LDSC_cross_trait_intercept"] = data.gcov_intercept
    data["LDSC_cross_trait_intercept_se"] = data.gcov_intercept_se
    data["source_GWAS_ID_sleep"] = data.sleep_trait.map(trait_lookup.source_id)
    data["source_GWAS_ID_disorder"] = data.external_trait.map(trait_lookup.source_id)
    data["cohort_information"] = data.apply(
        lambda r: f"sleep={panel_lookup.loc[r.sleep_trait, 'source_note']}; "
        f"disorder={panel_lookup.loc[r.external_trait, 'source_note']}", axis=1
    )
    data["notes"] = data.QC_notes.fillna("") + data.source_substitution_caveat.fillna("").map(
        lambda x: f"; source substitution: {x}" if x else ""
    )
    outcols = [
        "sleep_trait", "external_trait", "brain_disorder", "brain_disorder_label",
        "rg", "se", "z", "p", "original_BH_FDR_q",
        "significance_under_original_396_family", "sleep_h2", "sleep_h2_se",
        "h2_disorder", "h2_disorder_se", "sleep_h2_intercept",
        "external_h2_intercept", "LDSC_cross_trait_intercept",
        "LDSC_cross_trait_intercept_se", "SNP_overlap", "analysis_tier",
        "interpretation_status", "ancestry", "sample_size", "sleep_sample_size",
        "sleep_measurement_category", "sleep_collection_mode", "source_GWAS_ID_sleep",
        "source_GWAS_ID_disorder", "external_dataset_version", "phenotype_source",
        "cohort_information", "source_substitution", "notes", "input_log",
    ]
    result = data[outcols].sort_values(
        ["sleep_trait", "external_trait"],
        key=lambda s: s.map({v: i for i, v in enumerate(sleep_traits if s.name == "sleep_trait" else BRAIN)}),
    ).reset_index(drop=True)
    if result.original_BH_FDR_q.lt(0).any() or result.original_BH_FDR_q.gt(1).any():
        raise ValueError("Invalid FDR value in locked source table")

    output = OUT / "results/global/brain6_72_locked.tsv"
    source_map = OUT / "results/global/brain6_72_figure_source.tsv"
    serialized = result.to_csv(sep="\t", index=False, na_rep="NA", lineterminator="\n")
    if args.refresh_figures:
        if not output.exists() or output.read_text(encoding="utf-8") != serialized:
            raise ValueError("Existing global table differs from the locked source; refusing figure refresh")
    else:
        write_new(output, serialized)
        write_new(source_map, serialized)
    (OUT / "results/global/figures").mkdir(parents=True, exist_ok=True)

    q = result.pivot(index="sleep_trait", columns="brain_disorder", values="original_BH_FDR_q").loc[sleep_traits, BRAIN]
    rg = result.pivot(index="sleep_trait", columns="brain_disorder", values="rg").loc[sleep_traits, BRAIN]
    sig = result.pivot(index="sleep_trait", columns="brain_disorder", values="significance_under_original_396_family").loc[sleep_traits, BRAIN]
    labels = [trait_lookup.loc[t, "label"] for t in sleep_traits]
    dlabels = [trait_lookup.loc[t, "label"] for t in BRAIN]

    fig, ax = plt.subplots(figsize=(10.4, 7.2))
    im = ax.imshow(rg, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(range(n_brain), dlabels, rotation=30, ha="right")
    ax.set_yticks(range(n_sleep), labels)
    ax.set_title("Global genetic correlations: sleep × brain disorders")
    fig.colorbar(im, ax=ax, label="Genetic correlation (rg)")
    save_figure(fig, "fig1a_global_rg_heatmap", args.refresh_figures)

    fig, ax = plt.subplots(figsize=(10.4, 7.2))
    im = ax.imshow(np.abs(rg), cmap="magma", vmin=0, vmax=0.5, aspect="auto")
    for i in range(n_sleep):
        for j in range(n_brain):
            magnitude = abs(float(rg.iloc[i, j]))
            color = best_contrast_text(im.cmap(im.norm(magnitude)))
            ax.text(j, i, f"{'+' if rg.iloc[i,j] >= 0 else '−'}{magnitude:.2f}",
                    ha="center", va="center", color=color, fontsize=8)
    ax.set_xticks(range(n_brain), dlabels, rotation=30, ha="right")
    ax.set_yticks(range(n_sleep), labels)
    ax.set_title("Effect-size magnitude (color) and direction/rg (labels)")
    fig.colorbar(im, ax=ax, label="Absolute genetic correlation |rg|")
    save_figure(fig, "fig1b_effect_size_heatmap", args.refresh_figures)

    fig, ax = plt.subplots(figsize=(10.4, 7.2))
    im = ax.imshow(rg, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    for i in range(n_sleep):
        for j in range(n_brain):
            ax.scatter(j, i, s=70, facecolors="none", edgecolors="black", linewidths=1.2) if sig.iloc[i,j] else None
    ax.set_xticks(range(n_brain), dlabels, rotation=30, ha="right")
    ax.set_yticks(range(n_sleep), labels)
    ax.set_title("Original 396-family significance overlay (ring = q < 0.05)")
    fig.colorbar(im, ax=ax, label="Genetic correlation (rg)")
    save_figure(fig, "fig1c_original_fdr_significance", args.refresh_figures)

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    for j, d in enumerate(BRAIN):
        ax.plot(range(n_sleep), rg[d].to_numpy(), marker="o", linewidth=1.4, label=trait_lookup.loc[d, "label"])
    ax.axhline(0, color="#333333", linewidth=.8)
    ax.set_xticks(range(n_sleep), labels, rotation=35, ha="right")
    ax.set_ylabel("Genetic correlation (rg)")
    ax.set_title("Sleep-trait profiles across six brain disorders")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.tight_layout()
    save_figure(fig, "fig1d_sleep_trait_profiles", args.refresh_figures)

    corr = rg.corr(method="pearson")
    corr.index.name = "brain_disorder"
    corr.columns.name = "brain_disorder"
    sim_path = OUT / "results/global/disorder_profile_similarity.tsv"
    similarity_text = corr.to_csv(sep="\t", na_rep="NA", lineterminator="\n")
    if args.refresh_figures:
        if not sim_path.exists() or sim_path.read_text(encoding="utf-8") != similarity_text:
            raise ValueError("Existing disorder-profile table differs; refusing figure refresh")
    else:
        write_new(sim_path, similarity_text)
    z = linkage(squareform(1 - corr.to_numpy(), checks=False), method="average")
    order = [BRAIN[i] for i in leaves_list(z)]
    ordered_corr = corr.loc[order, order]
    fig, ax2 = plt.subplots(figsize=(8.8, 7.2))
    im = ax2.imshow(ordered_corr, vmin=-1, vmax=1, cmap="RdBu_r")
    ax2.set_xticks(range(n_brain), [trait_lookup.loc[t, "label"] for t in order], rotation=35, ha="right")
    ax2.set_yticks(range(n_brain), [trait_lookup.loc[t, "label"] for t in order])
    ax2.set_title("Disorder sleep-genetic profile similarity\n(Pearson r; average-linkage order)")
    for i in range(n_brain):
        for j in range(n_brain):
            value = ordered_corr.iloc[i, j]
            ax2.text(j, i, f"{value:.2f}", ha="center", va="center",
                     color="white" if abs(value) > .58 else "black", fontsize=8)
    fig.colorbar(im, ax=ax2, label=f"Pearson correlation across {n_sleep} sleep traits")
    save_figure(fig, "fig1e_disorder_profile_comparison", args.refresh_figures)

    config = OUT / "config/brain6_locked_family.yaml"
    yaml = (
        "schema_version: 1\n"
        "family_name: original_locked_atlas_396\n"
        "source_traits_manifest: results/atlas/traits.tsv\n"
        "source_pair_table: results/atlas/trait_pairs.tsv\n"
        f"source_rg_matrix_sha256: {sha256(ROOT / 'results/tables/rg_matrix.tsv')}\n"
        f"source_trait_manifest_sha256: {sha256(ROOT / 'results/atlas/traits.tsv')}\n"
        f"source_pair_table_sha256: {sha256(ROOT / 'results/atlas/trait_pairs.tsv')}\n"
        f"sleep_traits: [{', '.join(sleep_traits)}]\n"
        f"brain_disorders: [{', '.join(BRAIN)}]\n"
        f"expected_brain_rows: {expected_rows}\n"
        "fdr_source: global_rg_fdr_all_396\n"
        "fdr_threshold: 0.05\n"
        "recompute_brain6_fdr: false\n"
        "selection_timing: post_atlas_prioritized_deep_follow_up\n"
        "note: Preserve original 396-test significance; profile statistics are descriptive.\n"
    )
    if args.refresh_figures:
        if not config.exists() or config.read_text(encoding="utf-8") != yaml:
            raise ValueError("Existing family config differs; refusing figure refresh")
    else:
        write_new(config, yaml)
        write_new(config.with_suffix(".yaml.sha256"), f"{sha256(config)}  {config.name}\n")
    print(json.dumps({
        "rows": len(result), "original_396_significant": int(result.significance_under_original_396_family.sum()),
        "alz_significant": int(result.loc[result.brain_disorder.eq("alz"), "significance_under_original_396_family"].sum()),
        "output": str(output), "output_sha256": sha256(output), "config_sha256": sha256(config),
        "profile_similarity": str(sim_path),
        "figure_provenance": str(record_figure_provenance(args.refresh_figures)),
    }, indent=2))


if __name__ == "__main__":
    main()
