#!/usr/bin/env python3
"""Render an explicitly partial PLACO candidate-interval overview.

This figure is descriptive: four pair-QC-passed PLACO outputs are available,
while the protected Track B pair is absent. The displayed intervals are not a
complete-family set of independent or causal loci.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
LOCUS_PATH = ROOT / "brain6/results/loci/placo_candidate_loci_partial.tsv"
PROVENANCE_PATH = ROOT / "brain6/results/loci/placo_candidate_loci_partial.provenance.json"
FAMILY_PATH = ROOT / "brain6/config/placo_family_v3/family_lock.json"
OUTPUT_DIR = ROOT / "brain6/results/placo/figures"
EXPECTED_PAIRS = (
    "insomnia__mdd", "longsleep__scz", "longsleep__bipolar", "longsleep__parkinson",
)
FIGURE_STEM = "figS_partial_placo_candidate_intervals_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows(locus_path: Path = LOCUS_PATH,
              provenance_path: Path = PROVENANCE_PATH,
              family_path: Path = FAMILY_PATH) -> list[dict[str, str]]:
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("status") != "PASS_PARTIAL_FAMILY" or provenance.get("family_complete") is not False:
        raise ValueError("PLACO candidate artifact is not explicitly marked as a partial family")
    if provenance.get("protected_track_b_present") is not False:
        raise ValueError("Protected Track B state is inconsistent with this partial figure")
    locus_record = next((record for record in provenance.get("outputs", [])
                         if Path(record.get("path", "")).name == locus_path.name), None)
    if not locus_record or locus_record.get("sha256") != sha256(locus_path):
        raise ValueError("Candidate-locus table hash does not match its provenance")
    family = json.loads(family_path.read_text(encoding="utf-8"))
    expected = sorted(family["pairs"])
    # Pair identifiers are an explicit visual scope, independent of row order.
    if sorted(EXPECTED_PAIRS) != expected:
        raise ValueError("The locked PLACO pair family differs from the expected four unprotected pairs")

    with locus_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("No partial PLACO candidate loci are available")
    observed = {row["pair_id"] for row in rows}
    if observed != set(EXPECTED_PAIRS):
        raise ValueError(f"Unexpected partial candidate pair set: {sorted(observed)}")
    seen: set[str] = set()
    for row in rows:
        if row["locus_status"] != "PLACO_ONLY_CANDIDATE_LOCUS":
            raise ValueError(f"Unexpected locus interpretation: {row['locus_id']}")
        if row["locus_id"] in seen:
            raise ValueError(f"Duplicate candidate locus id: {row['locus_id']}")
        seen.add(row["locus_id"])
        chrom, start, stop = int(row["CHR"]), int(row["START"]), int(row["STOP"])
        p_value = float(row["lead_P_PLACO"])
        n_variants = int(row["n_candidate_variants"])
        if not 1 <= chrom <= 22 or not 1 <= start <= stop or not 0 < p_value <= 1 or n_variants < 1:
            raise ValueError(f"Invalid candidate-locus coordinates or statistics: {row['locus_id']}")
    return rows


def render(rows: list[dict[str, str]], output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path, pdf_path = output_dir / f"{FIGURE_STEM}.png", output_dir / f"{FIGURE_STEM}.pdf"
    if png_path.exists() or pdf_path.exists():
        raise FileExistsError("Refusing to overwrite an existing candidate-interval figure")

    pairs = list(EXPECTED_PAIRS)
    labels = {
        "insomnia__mdd": "Insomnia × MDD",
        "longsleep__scz": "Long sleep × schizophrenia",
        "longsleep__bipolar": "Long sleep × bipolar disorder",
        "longsleep__parkinson": "Long sleep × Parkinson’s disease",
    }
    figure, axes = plt.subplots(4, 1, figsize=(12.5, 8.2), sharex=True,
                                gridspec_kw={"hspace": 0.34})
    max_bp = max(int(row["STOP"]) for row in rows)
    values = [-math.log10(float(row["lead_P_PLACO"])) for row in rows]
    norm = Normalize(vmin=min(values), vmax=max(values))
    cmap = plt.get_cmap("viridis")
    for axis, pair in zip(axes, pairs, strict=True):
        subset = [row for row in rows if row["pair_id"] == pair]
        chromosomes = sorted({int(row["CHR"]) for row in subset}, reverse=True)
        y_position = {chrom: index for index, chrom in enumerate(chromosomes)}
        for row in subset:
            y = y_position[int(row["CHR"])]
            left, right = int(row["START"])/1e6, int(row["STOP"])/1e6
            x_lead = (int(row["START"]) + int(row["STOP"])) / 2e6
            color = cmap(norm(-math.log10(float(row["lead_P_PLACO"]))))
            axis.plot([left, right], [y, y], color=color, linewidth=3.4, alpha=0.9,
                      solid_capstyle="round", zorder=2)
            count = int(row["n_candidate_variants"])
            axis.scatter([x_lead], [y], s=22 + 10*math.sqrt(count), color=[color], edgecolor="white",
                         linewidth=0.7, zorder=3)
        axis.set_yticks(range(len(chromosomes)), [str(chrom) for chrom in chromosomes], fontsize=8)
        axis.set_ylim(len(chromosomes)-0.45, -0.55)
        axis.set_ylabel("Chromosome", fontsize=8, labelpad=7)
        axis.text(0.0, 1.06, labels[pair], transform=axis.transAxes, ha="left", va="bottom",
                  fontsize=9.5, fontweight="bold")
        axis.grid(axis="x", color="#d8dee8", linewidth=0.55)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(axis="y", length=0, pad=5)
        axis.set_xlim(0, max_bp / 1e6 * 1.02)
    axes[-1].set_xlabel("GRCh37 chromosome position (Mb; within-chromosome scale)", fontsize=9)
    axes[-1].set_xticks(range(0, int(max_bp/1e6)+25, 25))
    figure.suptitle("PLACO candidate intervals in four available Brain6 pairs",
                    x=0.12, y=0.985, ha="left", fontsize=14, fontweight="bold")
    figure.text(0.12, 0.955,
                "Partial family: protected Track B is unavailable. Intervals are PLACO-only candidates, not final independent or causal loci.",
                ha="left", va="top", fontsize=9, color="#414b5a")
    colorbar = figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                               fraction=0.025, pad=0.02)
    colorbar.set_label("−log₁₀ PLACO P for the best lead signal", fontsize=8.5)
    colorbar.ax.tick_params(labelsize=8)
    figure.text(0.12, 0.012,
                "Bars show candidate intervals using frozen 500-kb lead flanks (overlapping regions merged); point area reflects candidate-variant count. The across-pair family is incomplete.",
                ha="left", va="bottom", fontsize=8, color="#414b5a")
    figure.subplots_adjust(left=0.12, right=0.91, top=0.90, bottom=0.07)
    figure.savefig(png_path, dpi=300, facecolor="white", metadata={"Software": "Matplotlib"})
    figure.savefig(pdf_path, facecolor="white", metadata={"Creator": "Brain6 candidate-locus figure"})
    plt.close(figure)
    return png_path, pdf_path


def validate_existing(output_dir: Path) -> dict | None:
    png_path = output_dir / f"{FIGURE_STEM}.png"
    pdf_path = output_dir / f"{FIGURE_STEM}.pdf"
    provenance_path = output_dir / f"{FIGURE_STEM}.provenance.json"
    paths = (png_path, pdf_path, provenance_path)
    if not any(path.exists() for path in paths):
        return None
    if not all(path.is_file() for path in paths):
        raise ValueError("Partial candidate-figure output set found; refusing overwrite or repair")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    expected = {
        str(LOCUS_PATH.relative_to(ROOT)): sha256(LOCUS_PATH),
        str(PROVENANCE_PATH.relative_to(ROOT)): sha256(PROVENANCE_PATH),
        str(FAMILY_PATH.relative_to(ROOT)): sha256(FAMILY_PATH),
        str(Path(__file__).resolve().relative_to(ROOT)): sha256(Path(__file__)),
    }
    if provenance.get("status") != "PASS_PARTIAL_FAMILY" or provenance.get("inputs") != expected:
        raise ValueError("Existing figure provenance is stale or does not match its sources")
    recorded = {Path(row["path"]).name: row for row in provenance.get("outputs", [])}
    for path in (png_path, pdf_path):
        record = recorded.get(path.name)
        if not record or record["bytes"] != path.stat().st_size or record["sha256"] != sha256(path):
            raise ValueError(f"Existing figure output hash mismatch: {path}")
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    rows = load_rows()
    existing = validate_existing(args.output_dir)
    if existing:
        print(f"Validated existing partial candidate figure: {FIGURE_STEM} "
              f"({existing['candidate_locus_count']} intervals)")
        return
    png, pdf = render(rows, args.output_dir)
    inputs = {
        str(LOCUS_PATH.relative_to(ROOT)): sha256(LOCUS_PATH),
        str(PROVENANCE_PATH.relative_to(ROOT)): sha256(PROVENANCE_PATH),
        str(FAMILY_PATH.relative_to(ROOT)): sha256(FAMILY_PATH),
        str(Path(__file__).resolve().relative_to(ROOT)): sha256(Path(__file__)),
    }
    provenance = {
        "status": "PASS_PARTIAL_FAMILY",
        "figure_id": FIGURE_STEM,
        "scope": "Four pair-QC-passed PLACO outputs; protected insomnia-ADHD Track B output is unavailable.",
        "interpretation": "Descriptive PLACO-only candidate intervals; not a complete-family independent-locus result, mechanism, or causal claim.",
        "source_table": str(LOCUS_PATH.relative_to(ROOT)),
        "source_table_sha256": sha256(LOCUS_PATH),
        "source_locus_provenance_sha256": sha256(PROVENANCE_PATH),
        "family_lock_sha256": sha256(FAMILY_PATH),
        "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__)),
        "inputs": inputs,
        "candidate_locus_count": len(rows),
        "pairs": list(EXPECTED_PAIRS),
        "outputs": [{"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                     "sha256": sha256(path)} for path in (png, pdf)],
    }
    provenance_path = args.output_dir / f"{FIGURE_STEM}.provenance.json"
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=args.output_dir,
                                     prefix=f".{provenance_path.name}.", delete=False) as stream:
        temp_path = Path(stream.name)
        stream.write(json.dumps(provenance, indent=2) + "\n")
    os.replace(temp_path, provenance_path)
    print(f"Wrote {len(rows)} partial candidate intervals: {png} and {pdf}")


if __name__ == "__main__":
    main()
