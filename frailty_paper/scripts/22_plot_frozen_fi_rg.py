#!/usr/bin/env python3
"""Plot the existing frozen sleep–FI LDSC extract without refitting estimates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


LABELS = {
    "insomnia": "Insomnia",
    "sleepdur": "Sleep duration (continuous)",
    "shortsleep": "Short sleep (<7 h)",
    "longsleep": "Long sleep (≥9 h)",
    "chronotype": "Chronotype (morningness)",
    "sleepiness": "Daytime sleepiness",
    "napping": "Daytime napping",
    "snoring": "Snoring",
    "sleep_apnea": "Sleep apnea",
    "sleep_efficiency": "Actigraphy sleep efficiency",
    "accel_sleep_duration": "Actigraphy sleep duration",
    "sleep_timing": "Actigraphy sleep timing",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 12:
        raise ValueError(f"Expected the 12-row frozen FI extract, found {len(rows)}")
    unknown = sorted({row["sleep_trait"] for row in rows} - LABELS.keys())
    if unknown:
        raise ValueError(f"No display labels defined for: {', '.join(unknown)}")
    return rows


def read_environment_pins(path: Path) -> dict[str, str]:
    pins: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["component"] in {"python_workflow", "matplotlib"}:
                pins[row["component"]] = row["version_or_commit"]
    return pins


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path)
    parser.add_argument(
        "--allow-unpinned", action="store_true",
        help="Allow a clearly provenance-labeled preview when runtime versions differ from environment/tool_versions.tsv.",
    )
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    args.outdir = args.outdir or repo_root / "frailty_paper/analysis"
    args.outdir.mkdir(parents=True, exist_ok=True)

    manifest_path = repo_root / "frailty_paper/manifests/frozen_atlas_frailty_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("artifact_role") != "READ_ONLY_EXTRACT_OF_FROZEN_CORE_ATLAS_RESULTS_NOT_NEW_ANALYSIS":
        raise SystemExit("Unexpected frozen FI extract manifest role")
    input_path = repo_root / manifest["extracted_table"]
    input_sha256 = sha256(input_path)
    if input_sha256 != manifest.get("extracted_table_sha256"):
        raise SystemExit("Frozen FI extract SHA-256 does not match its audit manifest")
    if manifest.get("frailty_pair_count") != 12 or "396" not in manifest.get("fdr_denominator", ""):
        raise SystemExit("Frozen FI extract manifest has an unexpected row count or q-value family")

    pins_path = repo_root / "environment/tool_versions.tsv"
    pins = read_environment_pins(pins_path)
    actual_versions = {"python_workflow": platform.python_version(), "matplotlib": matplotlib.__version__}
    mismatches = {
        component: (pins.get(component), actual)
        for component, actual in actual_versions.items()
        if pins.get(component) != actual
    }
    if mismatches and not args.allow_unpinned:
        expected = ", ".join(f"{name}={value[0]} (found {value[1]})" for name, value in mismatches.items())
        raise SystemExit(f"Pinned figure environment mismatch: {expected}. Use the locked environment or --allow-unpinned for a preview.")

    rows = read_rows(input_path)
    names = [LABELS[row["sleep_trait"]] for row in rows]
    rg = np.array([float(row["global_rg"]) for row in rows])
    se = np.array([float(row["global_rg_se"]) for row in rows])
    q = np.array([float(row["global_rg_fdr_all_396"]) for row in rows])
    significant = q <= 0.05
    y = np.arange(len(rows))
    ci = 1.96 * se

    fig, (ax, qax) = plt.subplots(
        ncols=2, figsize=(12.5, 7.8), sharey=True,
        gridspec_kw={"width_ratios": [5.8, 1.0], "wspace": 0.02},
    )
    fig.subplots_adjust(left=0.30, right=0.97, top=0.88, bottom=0.16)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axvline(0, color="#687386", linewidth=1.1, linestyle=(0, (4, 3)), zorder=1)
    for i in range(len(rows)):
        color = "#145C72" if significant[i] else "#C17A2C"
        ax.errorbar(
            rg[i], y[i], xerr=ci[i], fmt="o", color=color,
            markerfacecolor=color if significant[i] else "white",
            markeredgewidth=1.8, markersize=7, capsize=3, linewidth=1.7,
            elinewidth=1.5, zorder=3,
        )
        qax.text(0.06, y[i], f"{q[i]:.2g}", ha="left", va="center",
                 fontsize=8.5, color="#344054")

    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xlim(-0.36, 0.76)
    ax.set_xlabel("Global genetic correlation (rg), with approximate 95% CI", labelpad=9)
    fig.suptitle("Sleep phenotypes and the Frailty Index", x=0.30, y=0.97,
                 ha="left", weight="bold", fontsize=16)
    fig.text(0.30, 0.925,
             "Frozen sleep-atlas estimates; q-values use the original 396-pair family",
             ha="left", va="bottom", fontsize=9.5, color="#475467")
    qax.set_xlim(0, 1)
    qax.set_ylim(ax.get_ylim())
    qax.axis("off")
    qax.text(0.06, 1.01, "All-396 q", transform=qax.transAxes, ha="left",
             va="bottom", fontsize=8.5, weight="bold", color="#344054")
    ax.grid(axis="x", color="#E6EAF0", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#98A2B3")
    ax.tick_params(axis="y", length=0, pad=8, labelsize=9)
    ax.tick_params(axis="x", colors="#475467")
    ax.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="none", color="#145C72",
                       markerfacecolor="#145C72", label="All-396 q ≤ 0.05"),
            plt.Line2D([], [], marker="o", linestyle="none", color="#C17A2C",
                       markerfacecolor="white", markeredgewidth=1.8,
                       label="All-396 q > 0.05"),
        ],
        loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2,
        frameon=False, fontsize=8.5,
    )
    fig.text(0.30, 0.015,
             "Descriptive extract only; not a new LDSC analysis or independent reprocessing. "
             "CI = rg ± 1.96 × SE. Direction depends on the source phenotype coding.",
             ha="left", va="bottom", fontsize=8, color="#475467")

    png = args.outdir / "preliminary_frozen_fi_global_rg_forest.png"
    pdf = args.outdir / "preliminary_frozen_fi_global_rg_forest.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight", metadata={
        "Title": "Sleep phenotypes and the Frailty Index: frozen atlas estimates",
        "Subject": "Descriptive extraction; original all-396-pair q-values retained",
        "Creator": f"22_plot_frozen_fi_rg.py; matplotlib {matplotlib.__version__}",
    })
    plt.close(fig)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": "python " + Path(__file__).resolve().relative_to(repo_root).as_posix()
        + (" " + " ".join(sys.argv[1:]) if sys.argv[1:] else ""),
        "input": input_path.relative_to(repo_root).as_posix(),
        "input_sha256": input_sha256,
        "frozen_audit_manifest": manifest_path.relative_to(repo_root).as_posix(),
        "frozen_audit_manifest_sha256": sha256(manifest_path),
        "registered_input_sha256": manifest["extracted_table_sha256"],
        "script": Path(__file__).resolve().relative_to(repo_root).as_posix(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "environment_pin_file": pins_path.relative_to(repo_root).as_posix(),
        "environment_pins": pins,
        "environment_gate_passed": not mismatches,
        "preview_only": bool(mismatches),
        "environment_mismatches": mismatches,
        "rows": len(rows),
        "q_column": "global_rg_fdr_all_396 (unchanged from frozen atlas extract)",
        "outputs": {path.resolve().relative_to(repo_root).as_posix(): sha256(path) for path in (png, pdf)},
        "limitations": [
            "Uses the frozen extract; does not regenerate LDSC or verify absent harmonized/munged source inputs.",
            "The 95% intervals are displayed as rg ± 1.96 × SE.",
            "Sleep-trait direction depends on the registered source phenotype coding.",
        ],
    }
    (args.outdir / "preliminary_frozen_fi_global_rg_forest.provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"png": str(png), "pdf": str(pdf), "rows": len(rows),
                      "registered_input_hash_match": True,
                      "input_sha256": provenance["input_sha256"]}))


if __name__ == "__main__":
    main()
