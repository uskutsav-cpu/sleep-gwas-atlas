#!/usr/bin/env python3
"""Plot the one frozen external global-rg replication, without re-estimating it."""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/replication/figures"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    source = ROOT / "results/track_b/replication/pair_b_ldsc.tsv"
    provenance_path = ROOT / "results/track_b/replication/pair_b_ldsc.provenance.json"
    global_path = ROOT / "brain6/results/global/brain6_72_locked.tsv"
    gwas_path = ROOT / "brain6/manifests/gwas_master.tsv"
    lock_path = ROOT / "results/track_b/replication_source.lock.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if sha256(source) != provenance["output_sha256"]:
        raise SystemExit("Archived replication table hash differs from its frozen provenance")
    if sha256(lock_path) != provenance["input_sha256"]["results/track_b/replication_source.lock.json"]:
        raise SystemExit("Replication source lock hash differs from its frozen provenance")

    replicated = read_tsv(source)
    replicated = [r for r in replicated if r["pair_id"] == "B"]
    discovery = [r for r in read_tsv(global_path)
                 if r["sleep_trait"] == "insomnia" and r["brain_disorder"] == "adhd"]
    if len(replicated) != 1 or len(discovery) != 1:
        raise SystemExit("Expected exactly one locked insomnia–ADHD estimate in each source table")
    rep, disc = replicated[0], discovery[0]
    finngen = [r for r in read_tsv(gwas_path) if r["trait"] == "adhd_replication_finngen_r13"]
    if len(finngen) != 1 or finngen[0]["GWAS_accession"] != rep["replication_source_id"]:
        raise SystemExit("FinnGen endpoint metadata does not match the archived replication source")
    if int(finngen[0]["cases"]) + int(finngen[0]["controls"]) != int(finngen[0]["sample_size"]):
        raise SystemExit("FinnGen case/control counts do not sum to the recorded sample size")
    if rep["direction_vs_discovery"] != "CONCORDANT" or rep["replication_class"] != "DIRECTIONAL_REPLICATION":
        raise SystemExit("Archived replication classification is not the expected directional result")

    estimates = np.array([float(disc["rg"]), float(rep["rg"])])
    ses = np.array([float(disc["se"]), float(rep["rg_se"])])
    labels = ["Locked atlas discovery\nDemontis et al. 2023 ADHD GWAS",
              f"FinnGen R13 external ADHD cohort\n{int(finngen[0]['cases']):,} cases; {int(finngen[0]['controls']):,} controls"]
    colors = ["#355C7D", "#C45B35"]
    y = np.array([1, 0])
    lo, hi = estimates - 1.96 * ses, estimates + 1.96 * ses

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.titlesize": 15, "axes.labelsize": 11})
    fig = plt.figure(figsize=(9.3, 5.0), facecolor="white")
    ax = fig.add_axes([0.34, 0.34, 0.43, 0.43])
    for i in range(2):
        ax.errorbar(estimates[i], y[i],
                    xerr=[[estimates[i] - lo[i]], [hi[i] - estimates[i]]],
                    fmt="o", markersize=9, markeredgecolor="white", markeredgewidth=1.1,
                    color=colors[i], ecolor=colors[i], elinewidth=2.2, capsize=4, zorder=3)
        ax.text(0.635, y[i],
                f"rᵍ = {estimates[i]:.4f}  ({lo[i]:.3f}, {hi[i]:.3f})",
                va="center", ha="left", fontsize=10, color="#252525")

    ax.axvline(0, color="#333333", linewidth=1, zorder=1)
    ax.set_yticks(y, labels)
    ax.set_xlim(-0.04, 0.96)
    ax.set_ylim(-0.75, 1.65)
    ax.set_xticks(np.arange(0, 0.61, 0.1))
    ax.set_xlabel("Global genetic correlation (rᵍ), with pointwise 95% CI")
    fig.text(0.34, 0.94, "External global-rg check for insomnia–ADHD",
             ha="left", va="top", fontsize=15, weight="bold", color="#202124")
    fig.text(0.34, 0.87,
             "The same insomnia GWAS is used in both rows; the external cohort replaces the ADHD GWAS.",
             ha="left", va="top", fontsize=9.3, color="#484848")
    ax.grid(axis="x", color="#D9DEE5", linewidth=0.8, zorder=0)
    ax.tick_params(axis="y", length=0, pad=11)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#8A8F98")
    fig.text(0.34, 0.13,
            "Directional concordance only. Finnish founder-population and phenotype differences apply;\n"
            "individual-level overlap was not verified. This is not a fully independent two-trait replication.",
            ha="left", va="top", fontsize=9.1, color="#484848",
            bbox={"boxstyle": "round,pad=0.55", "facecolor": "#F4F6F8", "edgecolor": "#DDE2E8"})

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        path = OUT / f"fig_replication_insomnia_adhd.{ext}"
        fig.savefig(path, dpi=300 if ext == "png" else None, bbox_inches="tight",
                    facecolor="white", metadata={"Title": "Insomnia-ADHD external global-rg replication",
                                                  "CreationDate": None, "ModDate": None})
    plt.close(fig)
    output_paths = [OUT / "fig_replication_insomnia_adhd.png", OUT / "fig_replication_insomnia_adhd.pdf"]
    source_paths = [source, provenance_path, lock_path, global_path, gwas_path, Path(__file__)]
    figure_provenance = {
        "schema_version": 1,
        "status": "PASS",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "pair": "insomnia__adhd",
        "replication_class": rep["replication_class"],
        "direction_vs_discovery": rep["direction_vs_discovery"],
        "same_sleep_gwas_used_for_both_estimates": True,
        "confidence_interval": "pointwise normal 95% interval: rg +/- 1.96 * SE",
        "estimates": [
            {"evidence_stage": "LOCKED_GLOBAL_DISCOVERY", "rg": estimates[0], "se": ses[0],
             "ci95_low": float(lo[0]), "ci95_high": float(hi[0])},
            {"evidence_stage": "EXTERNAL_ADHD_COHORT_DIRECTIONAL_REPLICATION", "rg": estimates[1], "se": ses[1],
             "ci95_low": float(lo[1]), "ci95_high": float(hi[1]),
             "cases": int(finngen[0]["cases"]), "controls": int(finngen[0]["controls"])},
        ],
        "inputs": [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": sha256(p)}
                   for p in source_paths],
        "outputs": [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": sha256(p)}
                    for p in output_paths],
        "caveats": ["Finnish founder-population and phenotype differences",
                    "individual-level overlap not verified",
                    "same insomnia GWAS used; not fully independent two-trait replication"],
    }
    provenance_out = OUT / "fig_replication_insomnia_adhd.provenance.json"
    provenance_out.write_text(json.dumps(figure_provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "pair": "insomnia__adhd",
                      "replication_output_sha256": provenance["output_sha256"],
                      "figure_outputs": figure_provenance["outputs"],
                      "provenance": str(provenance_out)}, indent=2))


if __name__ == "__main__":
    main()
