#!/usr/bin/env python3
"""Plot sensitivity-only sleep–latent frailty LDSC estimates with 95% CIs."""
from __future__ import annotations
import csv, hashlib, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

SLEEP = ["insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype", "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency", "accel_sleep_duration", "sleep_timing"]
SLEEP_LABELS = ["Insomnia", "Sleep duration", "Short sleep", "Long sleep", "Chronotype", "Sleepiness", "Napping", "Snoring", "Sleep apnea", "Actigraphy efficiency", "Actigraphy duration", "Actigraphy timing"]
FACTORS = ["frailty_general", "frailty_factor_1", "frailty_factor_2", "frailty_factor_3", "frailty_factor_4", "frailty_factor_5", "frailty_factor_6"]
FACTOR_LABELS = ["General factor", "F1 Social support", "F2 Unhealthy lifestyle", "F3 Multimorbidity", "F4 Metabolic", "F5 Cognition", "F6 Disability"]
SOURCE = Path("frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv")
SCRIPT = Path(__file__).resolve()
REPO = SCRIPT.parents[2]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    source = REPO / SOURCE
    with source.open(newline="", encoding="utf-8") as f: rows = list(csv.DictReader(f, delimiter="\t"))
    expected = {(s, d) for s in SLEEP for d in FACTORS}
    indexed = {(r["sleep_trait"], r["disease_trait"]): r for r in rows}
    if len(rows) != 84 or set(indexed) != expected: raise SystemExit("Expected exact locked 12×7 latent-factor family")
    if any(r["analysis_family"] != "secondary_sleep_x_latent_frailty" or r["family_denominator"] != "84" or r["interpretation_status"] != "SENSITIVITY_ONLY" for r in rows): raise SystemExit("Latent-family lock/status mismatch")
    out = REPO / "frailty_paper/analysis"
    out.mkdir(parents=True, exist_ok=True)
    # Keep the complete plotted data and transparent Wald intervals beside the figure.
    data_path = out / "figure3_latent_factor_forest_source_data.tsv"
    fields = ["sleep_trait", "sleep_label", "frailty_factor", "frailty_label", "rg", "se", "ci95_lower", "ci95_upper", "p_value", "q_value", "q_family", "family_denominator", "cohort_overlap_status", "interpretation_status", "claim_limit"]
    with data_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader()
        for s, sl in zip(SLEEP, SLEEP_LABELS):
            for factor, fl in zip(FACTORS, FACTOR_LABELS):
                r = indexed[s, factor]; rg, se = float(r["rg"]), float(r["se"])
                w.writerow({"sleep_trait":s,"sleep_label":sl,"frailty_factor":factor,"frailty_label":fl,"rg":r["rg"],"se":r["se"],"ci95_lower":f"{rg-1.96*se:.6f}","ci95_upper":f"{rg+1.96*se:.6f}","p_value":r["p_value"],"q_value":r["q_value"],"q_family":"BH across fixed 84-pair secondary family","family_denominator":"84","cohort_overlap_status":r["cohort_overlap_status"],"interpretation_status":r["interpretation_status"],"claim_limit":r["claim_limit"]})
    fig, axes = plt.subplots(1, 7, figsize=(18, 9), sharey=True)
    fig.patch.set_facecolor("white")
    colors = ["#145C72", "#AD5D33", "#527A43", "#75558A", "#9A6C14", "#3D718A", "#A14C59"]
    for j, (ax, factor, title) in enumerate(zip(axes, FACTORS, FACTOR_LABELS)):
        ax.axvline(0, color="#7B8794", lw=0.8, zorder=0)
        ax.set_title(title, fontsize=10, weight="bold", pad=10, color="#182230")
        for i, s in enumerate(SLEEP):
            r = indexed[s, factor]; x, se = float(r["rg"]), float(r["se"])
            ax.errorbar(x, i, xerr=1.96*se, fmt="o", ms=4.2, capsize=1.5, lw=1.1, color=colors[j], ecolor=colors[j], zorder=2)
            if float(r["q_value"]) <= .05:
                ax.text(x + (1.96*se) + .018, i, "*", va="center", ha="left", fontsize=9, weight="bold", color="#182230")
        ax.set_xlim(-.45, .85); ax.set_xticks([-.4, 0, .4, .8]); ax.tick_params(axis="x", labelsize=7, length=3, colors="#475467")
        ax.set_ylim(11.6, -.6); ax.grid(axis="y", color="#E6EAF0", lw=.5); ax.grid(axis="x", color="#EEF1F4", lw=.5)
        ax.spines[["top", "right", "left"]].set_visible(False); ax.spines["bottom"].set_color("#C5CDD6")
        if j == 0:
            ax.set_yticks(range(12), SLEEP_LABELS); ax.tick_params(axis="y", labelsize=8, length=0, colors="#344054")
        else: ax.tick_params(axis="y", length=0, labelleft=False)
    fig.suptitle("Sleep–frailty decomposition: latent-factor LDSC estimates", x=.055, y=.985, ha="left", fontsize=15, weight="bold", color="#182230")
    fig.text(.055, .951, "Secondary sensitivity family · points show rg; bars show Wald 95% CIs · * BH q ≤ 0.05 across 84 pairs", fontsize=9, color="#475467")
    fig.text(.055, .036, "Insomnia is one of the 30 deficits used to construct the general factor; that pairing has direct part–whole dependence. Exact participant overlap is unknown. Estimates do not establish independent replication, differences between correlated endpoints, or causation.", fontsize=8, color="#475467", wrap=True)
    fig.text(.055, .016, "Physical components, Fried phenotype and HFRS are not represented: eligible sleep-pair estimates are unavailable.", fontsize=8, color="#475467")
    fig.subplots_adjust(left=.16, right=.99, top=.88, bottom=.14, wspace=.18)
    paths = [out/"figure3_latent_factor_forest.png", out/"figure3_latent_factor_forest.pdf", out/"figure3_latent_factor_forest.svg"]
    fig.savefig(paths[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(paths[1], bbox_inches="tight", facecolor="white", metadata={"Title":"Sleep–frailty latent-factor LDSC sensitivity estimates","Subject":"Fixed 84-pair family; 95% Wald CIs; sensitivity-only"})
    fig.savefig(paths[2], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    manifest = {"created_utc":datetime.now(timezone.utc).isoformat(),"command":"python " + SCRIPT.relative_to(REPO).as_posix(),"script":SCRIPT.relative_to(REPO).as_posix(),"script_sha256":sha(SCRIPT),"python":platform.python_version(),"matplotlib":matplotlib.__version__,"input":SOURCE.as_posix(),"input_sha256":sha(source),"rows":84,"analysis_family":"secondary_sleep_x_latent_frailty","family_denominator":84,"interpretation_status":"SENSITIVITY_ONLY","interval":"rg ± 1.96 × SE (Wald 95% CI)","outputs":{p.relative_to(REPO).as_posix():sha(p) for p in paths+[data_path]},"limitations":["Exact participant overlap is unknown.","Insomnia contributes to construction of the general factor; this pair has part–whole dependence.","No direct tests compare correlated frailty endpoints.","Physical component, Fried and HFRS sleep-pair estimates are unavailable and are not shown.","Genetic correlation is not evidence of causation."]}
    manifest_path = out/"figure3_latent_factor_forest.provenance.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"rows":84,"outputs":[p.as_posix() for p in paths],"source_data":data_path.as_posix(),"provenance":manifest_path.as_posix()}))
if __name__ == "__main__": main()
