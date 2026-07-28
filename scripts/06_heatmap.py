#!/usr/bin/env python3
"""Figure 2: sleep/circadian x disease genetic correlation heatmap.

This is the direct analogue of Fig. 1a in Grotzinger/Werme et al. 2026 —
rg below the diagonal, asterisks for significance after multiple-testing
correction. Ours is rectangular (sleep on rows, disease on columns) rather
than triangular, because our two axes are different sets of traits.

    python3 06_heatmap.py --rg results/tables/rg_matrix.tsv \
        --config config/traits.tsv --out results/figures/fig2_rg_heatmap.png

Every figure carries a provenance footer: the source rg table, its mtime, and
a hash of its contents. A figure that gets separated from the table that made
it can then still be traced back — or exposed as untraceable. This exists
because a previous render of this script survived deletion of its inputs and
sat in results/figures/ with invented rg values on it, visually
indistinguishable from a real Figure 2.

Synthetic input additionally requires --provenance-label, which stamps a
diagonal watermark across the plot area. The footer is for auditing; the
watermark is because footers get cropped and screenshots get pasted into
slide decks.
"""
import argparse
import hashlib
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DOMAIN_ORDER = ["psychiatric", "neuro", "metabolic", "cardio", "immune",
                "aging", "cancer"]

# Any rg table living under a directory with this name is synthetic by
# construction (see make_test_data.py --logs).
SMOKETEST_DIR = "_smoketest"


def provenance(path):
    """source path + mtime + short content hash, as one footer line."""
    p = Path(path).resolve()
    h = hashlib.sha256(p.read_bytes()).hexdigest()[:12]
    mt = datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    return f"source: {p}  |  mtime: {mt}  |  sha256:{h}"


def stars(fdr):
    if pd.isna(fdr):
        return ""
    return "***" if fdr < 1e-3 else "**" if fdr < 1e-2 else "*" if fdr < 0.05 else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rg", required=True)
    ap.add_argument("--config", default="config/traits.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--annot", choices=["rg", "none"], default="rg")
    ap.add_argument("--provenance-label", default=None,
                    help="text stamped diagonally across the plot, e.g. "
                         "'SYNTHETIC — NOT REAL DATA'. Required for any rg "
                         "table under results/_smoketest/.")
    a = ap.parse_args()

    src = Path(a.rg).resolve()
    if SMOKETEST_DIR in src.parts and not a.provenance_label:
        sys.exit(
            f"ERROR: --rg is under {SMOKETEST_DIR}/, so its numbers are "
            "synthetic, but no --provenance-label was given.\n"
            "       Refusing to render an unlabelled figure from fake input.\n"
            '       Re-run with: --provenance-label "SYNTHETIC — NOT REAL DATA"'
        )

    df = pd.read_csv(a.rg, sep="\t")
    cfg = pd.read_csv(a.config, sep="\t", dtype=str).set_index("trait_id")
    lab = cfg["label"].to_dict()
    dom = cfg["domain"].to_dict()

    mat = df.pivot_table(index="sleep_trait", columns="disease_trait",
                         values="rg", aggfunc="first")
    fdr = df.pivot_table(index="sleep_trait", columns="disease_trait",
                         values="fdr", aggfunc="first").reindex_like(mat)

    # order columns by domain so the block structure is visible
    cols = sorted(mat.columns,
                  key=lambda c: (DOMAIN_ORDER.index(dom.get(c, "cancer"))
                                 if dom.get(c) in DOMAIN_ORDER else 99, c))
    mat, fdr = mat[cols], fdr[cols]

    vmax = float(np.nanmax(np.abs(mat.values))) if mat.size else 1.0
    vmax = max(vmax, 0.1)

    fig, ax = plt.subplots(figsize=(1.05 * len(mat.columns) + 3.2,
                                    0.72 * len(mat.index) + 2.6))
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")

    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels([lab.get(c, c) for c in mat.columns],
                       rotation=40, ha="right", fontsize=9)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels([lab.get(r, r) for r in mat.index], fontsize=9)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if pd.isna(v):
                ax.text(j, i, "·", ha="center", va="center", color="0.6")
                continue
            txt = (f"{v:.2f}\n{stars(fdr.values[i, j])}"
                   if a.annot == "rg" else stars(fdr.values[i, j]))
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.6 * vmax else "black")

    # domain separators
    prev = None
    for j, c in enumerate(mat.columns):
        d = dom.get(c)
        if prev is not None and d != prev:
            ax.axvline(j - 0.5, color="black", lw=1.4)
        prev = d

    ax.set_xticks(np.arange(-0.5, len(mat.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(mat.index), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)

    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("genetic correlation ($r_g$)", fontsize=9)
    ax.set_title("Sleep/circadian × disease genetic correlation (LDSC)\n"
                 "* FDR<0.05  ** FDR<0.01  *** FDR<0.001",
                 fontsize=11, pad=12)

    fig.tight_layout()

    # Watermark, tiled so that cropping to any single row still catches it.
    # Drawn after tight_layout so the axes has its final size, then measured
    # and shrunk to fit: a fixed point size is either clipped mid-word on a
    # narrow matrix or lost on a wide one.
    if a.provenance_label:
        wm = [ax.text(0.5, yf, a.provenance_label, transform=ax.transAxes,
                      rotation=22, ha="center", va="center", fontsize=26,
                      fontweight="bold", color="crimson", alpha=0.28,
                      zorder=5, clip_on=True)
              for yf in (0.16, 0.50, 0.84)]
        fig.canvas.draw()
        avail = ax.get_window_extent().width * 0.94
        got = wm[0].get_window_extent().width
        if got > avail:
            for t in wm:
                t.set_fontsize(max(8.0, 26.0 * avail / got))

    prov = provenance(a.rg)
    footer = f"{a.provenance_label}  |  {prov}" if a.provenance_label else prov
    fig.text(0.005, 0.002, footer, fontsize=5.5, family="monospace",
             color="crimson" if a.provenance_label else "0.4",
             ha="left", va="bottom")

    fig.savefig(a.out, dpi=200, bbox_inches="tight")
    print(f"wrote {a.out}  ({mat.shape[0]} sleep × {mat.shape[1]} disease traits)")
    print(f"  {prov}")
    if a.provenance_label:
        print(f"  watermarked: {a.provenance_label}")


if __name__ == "__main__":
    main()
