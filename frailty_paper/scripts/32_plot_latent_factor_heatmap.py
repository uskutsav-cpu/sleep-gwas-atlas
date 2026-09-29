#!/usr/bin/env python3
"""Plot the completed secondary 12 × 7 sleep–latent frailty LDSC family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SLEEP_ORDER = ["insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype", "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency", "accel_sleep_duration", "sleep_timing"]
SLEEP_LABELS = ["Insomnia", "Sleep duration", "Short sleep", "Long sleep", "Chronotype", "Sleepiness", "Napping", "Snoring", "Sleep apnea", "Actigraphy efficiency", "Actigraphy duration", "Actigraphy timing"]
FACTOR_ORDER = ["frailty_general", "frailty_factor_1", "frailty_factor_2", "frailty_factor_3", "frailty_factor_4", "frailty_factor_5", "frailty_factor_6"]
FACTOR_LABELS = ["General", "Factor 1\nLimited\nsocial support", "Factor 2\nUnhealthy\nlifestyle", "Factor 3\nMultimorbidity", "Factor 4\nMetabolic\nproblems", "Factor 5\nPoorer\ncognition", "Factor 6\nDisability"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def output_label(path: Path, repo: Path) -> str:
    """Return a portable repo-relative output path, or an absolute external path."""
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def load_family(path: Path) -> tuple[np.ndarray, np.ndarray, list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 84:
        raise ValueError(f"Expected exactly 84 rows in the fixed family, found {len(rows)}")
    pairs = {(r["sleep_trait"], r["disease_trait"]) for r in rows}
    expected = {(s, f) for s in SLEEP_ORDER for f in FACTOR_ORDER}
    if pairs != expected or len(pairs) != 84:
        raise ValueError(f"Input does not contain the exact 12 × 7 family; missing={sorted(expected-pairs)}, unexpected={sorted(pairs-expected)}")
    if any(r["analysis_family"] != "secondary_sleep_x_latent_frailty" or r["family_denominator"] != "84" for r in rows):
        raise ValueError("Unexpected analysis family or multiplicity denominator")
    if any(r["interpretation_status"] != "SENSITIVITY_ONLY" for r in rows):
        raise ValueError("Unexpected interpretation status; review family lock")
    indexed = {(r["sleep_trait"], r["disease_trait"]): r for r in rows}
    rg = np.array([[float(indexed[s, f]["rg"]) for f in FACTOR_ORDER] for s in SLEEP_ORDER])
    q = np.array([[float(indexed[s, f]["q_value"]) for f in FACTOR_ORDER] for s in SLEEP_ORDER])
    return rg, q, rows


def main() -> None:
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
    except ImportError as exc:
        raise SystemExit("Figure rendering requires Matplotlib; install the pinned dependencies from requirements-pipeline.txt") from exc
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    input_path = repo / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
    outdir = (args.outdir or repo / "frailty_paper/analysis").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    rg, q, rows = load_family(input_path)

    fig, ax = plt.subplots(figsize=(11.4, 8.0))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    image = ax.imshow(rg, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.7, vcenter=0, vmax=0.7), aspect="auto")
    ax.set_xticks(np.arange(len(FACTOR_LABELS)), FACTOR_LABELS)
    ax.set_yticks(np.arange(len(SLEEP_LABELS)), SLEEP_LABELS)
    ax.tick_params(axis="both", length=0, labelsize=9)
    ax.set_xticks(np.arange(-0.5, len(FACTOR_LABELS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(SLEEP_LABELS), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.spines[:].set_visible(False)
    for i in range(rg.shape[0]):
        for j in range(rg.shape[1]):
            star = "*" if q[i, j] <= 0.05 else ""
            color = "white" if abs(rg[i, j]) >= 0.42 else "#182230"
            ax.text(j, i, f"{rg[i, j]:.2f}{star}", ha="center", va="center", color=color, fontsize=8.5, weight="normal")
    ax.set_title("Sleep phenotypes × latent frailty dimensions", loc="left", fontsize=15, weight="bold", pad=20, color="#182230")
    ax.text(0, 1.015, "LDSC genetic correlation (rg); * BH q ≤ 0.05 across the fixed 84-pair family", transform=ax.transAxes, ha="left", va="bottom", fontsize=9, color="#475467")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    colorbar.set_label("Genetic correlation (rg)", fontsize=9, labelpad=9)
    colorbar.ax.tick_params(labelsize=8, length=0)
    fig.text(0.125, 0.035,
             "Secondary sensitivity analysis. Exact participant overlap with UK Biobank sleep GWAS is unknown; "
             "estimates are not independent replication or causal effects.",
             ha="left", va="bottom", fontsize=8, color="#475467")
    fig.subplots_adjust(left=0.23, right=0.91, top=0.87, bottom=0.20)

    png = outdir / "sleep_latent_frailty_rg_heatmap.png"
    pdf = outdir / "sleep_latent_frailty_rg_heatmap.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white", metadata={
        "Title": "Sleep phenotypes and latent frailty factors: LDSC genetic correlations",
        "Subject": "Secondary 84-pair family; q-values BH-corrected over all 84; sensitivity-only",
        "Creator": f"32_plot_latent_factor_heatmap.py; matplotlib {matplotlib.__version__}",
    })
    plt.close(fig)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": "python " + Path(__file__).resolve().relative_to(repo).as_posix() + (" " + " ".join(sys.argv[1:]) if sys.argv[1:] else ""),
        "input": input_path.relative_to(repo).as_posix(),
        "input_sha256": sha256(input_path),
        "script": Path(__file__).resolve().relative_to(repo).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "rows": len(rows),
        "analysis_family": "secondary_sleep_x_latent_frailty",
        "family_denominator": 84,
        "q_column": "q_value (BH over all 84 pairs)",
        "interpretation_status": "SENSITIVITY_ONLY",
        "outputs": {output_label(p, repo): sha256(p) for p in (png, pdf)},
        "limitations": ["Exact participant overlap is unknown for all UK Biobank-containing factor comparisons.", "Genetic correlation is not evidence of causation."],
    }
    manifest_path = outdir / "sleep_latent_frailty_rg_heatmap.provenance.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "png": str(png), "pdf": str(pdf), "manifest": str(manifest_path), "input_sha256": manifest["input_sha256"]}))


if __name__ == "__main__":
    main()
