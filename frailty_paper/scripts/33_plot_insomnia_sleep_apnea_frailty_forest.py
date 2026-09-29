#!/usr/bin/env python3
"""Plot available insomnia and sleep-apnea rg estimates across FI and factors."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path



SLEEP_TRAITS = ["insomnia", "sleep_apnea"]
ENDPOINTS = ["frailty", "frailty_general", "frailty_factor_1", "frailty_factor_2", "frailty_factor_3", "frailty_factor_4", "frailty_factor_5", "frailty_factor_6"]
LABELS = ["Frailty Index", "General factor", "Factor 1 · Limited social support", "Factor 2 · Unhealthy lifestyle", "Factor 3 · Multimorbidity", "Factor 4 · Metabolic problems", "Factor 5 · Poorer cognition", "Factor 6 · Disability"]
SLEEP_LABELS = {"insomnia": "Insomnia", "sleep_apnea": "Sleep apnea"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def collect(repo: Path) -> dict[str, list[dict[str, float | str]]]:
    fi_path = repo / "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv"
    latent_path = repo / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
    fi = read_table(fi_path)
    latent = read_table(latent_path)
    fi_rows = {r["sleep_trait"]: r for r in fi if r["non_sleep_trait"] == "frailty"}
    latent_rows = {(r["sleep_trait"], r["disease_trait"]): r for r in latent}
    if len(fi_rows) != 12 or not set(SLEEP_TRAITS).issubset(fi_rows):
        raise ValueError("Malformed FI table or target sleep traits missing")
    expected_latent = {(s, f) for s in SLEEP_TRAITS for f in ENDPOINTS[1:]}
    if len(latent_rows) != 84 or not expected_latent.issubset(latent_rows):
        raise ValueError("Malformed latent table or target sleep-factor estimates missing")
    output = {}
    for sleep in SLEEP_TRAITS:
        if sleep not in fi_rows:
            raise ValueError(f"Missing FI estimate for {sleep}")
        entries = []
        r = fi_rows[sleep]
        entries.append({"endpoint": "frailty", "rg": float(r["global_rg"]), "se": float(r["global_rg_se"]), "q": float(r["global_rg_fdr_all_396"]), "q_family": "all 396 frozen atlas pairs"})
        for endpoint in ENDPOINTS[1:]:
            r = latent_rows.get((sleep, endpoint))
            if r is None:
                raise ValueError(f"Missing latent estimate for {sleep}/{endpoint}")
            if r["family_denominator"] != "84" or r["interpretation_status"] != "SENSITIVITY_ONLY":
                raise ValueError(f"Unexpected latent result family/status for {sleep}/{endpoint}")
            entries.append({"endpoint": endpoint, "rg": float(r["rg"]), "se": float(r["se"]), "q": float(r["q_value"]), "q_family": "84 latent-factor pairs"})
        output[sleep] = entries
    return output


def resolve_output_dir(requested: Path | None, repo: Path) -> Path:
    """Resolve relative --outdir values from the caller's working directory."""
    return (requested or repo / "frailty_paper/analysis").resolve()


def main() -> None:
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise SystemExit("Figure rendering requires Matplotlib and NumPy; install the pinned dependencies from requirements-pipeline.txt") from exc
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    outdir = resolve_output_dir(args.outdir, repo)
    outdir.mkdir(parents=True, exist_ok=True)
    values = collect(repo)
    y = np.arange(len(ENDPOINTS))
    fig, axes = plt.subplots(ncols=2, figsize=(14, 8), sharey=True, gridspec_kw={"wspace": 0.04})
    fig.patch.set_facecolor("white")
    for ax, sleep in zip(axes, SLEEP_TRAITS):
        ax.set_facecolor("white")
        ax.axvline(0, color="#687386", linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
        for i, row in enumerate(values[sleep]):
            sig = row["q"] <= 0.05
            color = "#145C72" if sig else "#C17A2C"
            ax.errorbar(float(row["rg"]), y[i], xerr=1.96 * float(row["se"]), fmt="o", color=color,
                        markerfacecolor=color if sig else "white", markeredgewidth=1.5, markersize=6,
                        capsize=2.5, linewidth=1.4, zorder=3)
        ax.set_title(SLEEP_LABELS[sleep], loc="left", fontsize=12, weight="bold", color="#182230", pad=12)
        ax.grid(axis="x", color="#E6EAF0", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_color("#98A2B3")
        ax.tick_params(axis="x", colors="#475467", labelsize=9)
        ax.tick_params(axis="y", length=0, pad=8, labelsize=9)
        ax.set_xlim(-0.4, 0.85)
        ax.set_xlabel("Genetic correlation (rg ± 1.96 × SE)", fontsize=9, labelpad=8)
    axes[0].set_yticks(y, LABELS)
    axes[0].invert_yaxis()
    fig.suptitle("Insomnia and sleep apnea across available frailty dimensions", x=0.13, y=0.98, ha="left", weight="bold", fontsize=15, color="#182230")
    fig.text(0.13, 0.935, "FI q uses the frozen 396-pair atlas family; latent-factor q uses its fixed 84-pair family", ha="left", va="bottom", fontsize=9, color="#475467")
    fig.legend(handles=[
        plt.Line2D([], [], marker="o", linestyle="none", color="#145C72", markerfacecolor="#145C72", label="Within-family q ≤ 0.05"),
        plt.Line2D([], [], marker="o", linestyle="none", color="#C17A2C", markerfacecolor="white", markeredgewidth=1.5, label="Within-family q > 0.05"),
    ], loc="lower center", bbox_to_anchor=(0.52, 0.06), ncol=2, frameon=False, fontsize=8.5)
    fig.text(0.13, 0.015, "Descriptive cross-endpoint comparison only. UK Biobank overlap is expected or possible; exact participant intersections are unknown. "
             "Fried, HFRS and physical-component estimates are unavailable and are not represented. Estimates are not causal.",
             ha="left", va="bottom", fontsize=8, color="#475467")
    fig.subplots_adjust(left=0.31, right=0.98, top=0.86, bottom=0.20)

    png = outdir / "insomnia_sleep_apnea_frailty_dimensions_forest.png"
    pdf = outdir / "insomnia_sleep_apnea_frailty_dimensions_forest.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white", metadata={
        "Title": "Insomnia and sleep apnea across available frailty dimensions",
        "Subject": "Descriptive FI and latent-factor rg estimates with distinct locked q-value families",
        "Creator": f"33_plot_insomnia_sleep_apnea_frailty_forest.py; matplotlib {matplotlib.__version__}",
    })
    plt.close(fig)

    fi_path = repo / "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv"
    latent_path = repo / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": "python " + Path(__file__).resolve().relative_to(repo).as_posix() + (" " + " ".join(sys.argv[1:]) if sys.argv[1:] else ""),
        "inputs": {p.relative_to(repo).as_posix(): sha256(p) for p in [fi_path, latent_path]},
        "script": Path(__file__).resolve().relative_to(repo).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "sleep_traits": SLEEP_TRAITS,
        "endpoints": ENDPOINTS,
        "row_count": 16,
        "multiplicity": {"frailty": 396, "latent_factors": 84},
        "outputs": {p.relative_to(repo).as_posix(): sha256(p) for p in [png, pdf]},
        "limitations": ["The endpoints use different locked q-value families.", "Exact participant overlap is unknown.", "Fried, HFRS, and physical components are not available here.", "Forest comparisons are descriptive and are not tests of differences or causal effects."],
    }
    manifest_path = outdir / "insomnia_sleep_apnea_frailty_dimensions_forest.provenance.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"rows": 16, "png": str(png), "pdf": str(pdf), "manifest": str(manifest_path)}))


if __name__ == "__main__":
    main()
