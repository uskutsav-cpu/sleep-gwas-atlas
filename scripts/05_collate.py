#!/usr/bin/env python3
"""Turn LDSC .log files into tidy tables, and apply the h2 QC gate.

The gate is the important part. It decides which traits are powered enough to
carry GenomicSEM downstream, and its output is your first real scoping
recommendation to your mentor.

    python3 05_collate.py --mode h2 --logdir results/logs --out results/tables/h2_summary.tsv
    python3 05_collate.py --mode rg --logdir results/logs --out results/tables/rg_matrix.tsv
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

# --- QC gate thresholds -------------------------------------------------------
# Z = h2/se. Convention in the field is that rg is uninterpretable below about
# Z=4; LDSC's own docs suggest h2 Z > 4 before trusting genetic correlations.
Z_MIN = 4.0
# LDSC intercept should sit near 1. Well above 1 suggests confounding
# (stratification) rather than polygenicity.
INTERCEPT_MAX = 1.20

H2_PAT = re.compile(
    r"Total (Observed|Liability) scale h2:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)")
INT_PAT = re.compile(r"^Intercept:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)", re.M)
CHI_PAT = re.compile(r"Mean Chi\^2:\s*(-?[\d.eE+-]+)")
LAM_PAT = re.compile(r"Lambda GC:\s*(-?[\d.eE+-]+)")
RAT_PAT = re.compile(r"^Ratio:\s*(-?[\d.eE+-]+)", re.M)


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def parse_h2(logdir):
    rows = []
    for path in sorted(glob.glob(os.path.join(logdir, "h2_*.log"))):
        trait = os.path.basename(path)[3:-4]
        txt = open(path).read()
        m = H2_PAT.search(txt)
        if not m:
            rows.append({"trait": trait, "note": "NO h2 IN LOG - check for errors"})
            continue
        h2, se = f(m.group(2)), f(m.group(3))
        i = INT_PAT.search(txt)
        z = h2 / se if se else float("nan")
        rows.append({
            "trait": trait,
            "scale": m.group(1).lower(),
            "h2": h2, "se": se, "z": round(z, 2),
            "intercept": f(i.group(1)) if i else float("nan"),
            "intercept_se": f(i.group(2)) if i else float("nan"),
            "lambda_gc": f(LAM_PAT.search(txt).group(1)) if LAM_PAT.search(txt) else float("nan"),
            "mean_chi2": f(CHI_PAT.search(txt).group(1)) if CHI_PAT.search(txt) else float("nan"),
            "ratio": f(RAT_PAT.search(txt).group(1)) if RAT_PAT.search(txt) else float("nan"),
        })
    df = pd.DataFrame(rows)
    if "h2" not in df:
        return df
    df["pass_z"] = df["z"] >= Z_MIN
    df["pass_intercept"] = df["intercept"] <= INTERCEPT_MAX
    df["verdict"] = ["PASS" if (a and b) else "DROP"
                     for a, b in zip(df["pass_z"].fillna(False),
                                     df["pass_intercept"].fillna(False))]
    return df.sort_values("z", ascending=False)


def parse_rg(logdir):
    """LDSC prints a fixed-width table after 'Summary of Genetic Correlation'."""
    frames = []
    for path in sorted(glob.glob(os.path.join(logdir, "rg_*.log"))):
        lines = open(path).read().splitlines()
        try:
            i = next(k for k, l in enumerate(lines)
                     if "Summary of Genetic Correlation Results" in l)
        except StopIteration:
            print(f"  !! no results table in {path}")
            continue
        header = lines[i + 1].split()
        body = []
        for l in lines[i + 2:]:
            if not l.strip():
                break
            body.append(l.split())
        if body:
            frames.append(pd.DataFrame(body, columns=header))
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    for c in ["rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "gcov_int"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # strip paths -> bare trait ids
    for c in ["p1", "p2"]:
        df[c] = (df[c].str.replace(r".*/", "", regex=True)
                      .str.replace(".sumstats.gz", "", regex=False))
    df = df.rename(columns={"p1": "sleep_trait", "p2": "disease_trait"})
    # Benjamini-Hochberg across all pairs tested
    d = df.dropna(subset=["p"]).sort_values("p").copy()
    m = len(d)
    raw = d["p"].values * m / np.arange(1, m + 1)
    # BH is a step-UP procedure: enforce monotonicity by taking the running
    # minimum from the LARGEST p-value backwards, not forwards.
    d["fdr"] = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
    df = df.merge(d[["fdr"]], left_index=True, right_index=True, how="left")
    return df.sort_values("p")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["h2", "rg"], required=True)
    ap.add_argument("--logdir", default="results/logs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    df = parse_h2(a.logdir) if a.mode == "h2" else parse_rg(a.logdir)
    if df.empty:
        print("!! nothing parsed - check that LDSC actually ran")
        raise SystemExit(1)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    print(df.to_string(index=False))
    print(f"\nwrote {a.out}")
    if a.mode == "h2" and "verdict" in df:
        drop = df.loc[df.verdict == "DROP", "trait"].tolist()
        print(f"\nQC GATE (Z>={Z_MIN}, intercept<={INTERCEPT_MAX})")
        print(f"  PASS: {df.loc[df.verdict=='PASS','trait'].tolist()}")
        print(f"  DROP: {drop}")
        print("\n  -> Set status=PASS in config/traits.tsv for the survivors.")
        print("  -> The DROP list IS your meeting talking point: these traits")
        print("     cannot support a GenomicSEM factor, so say so explicitly.")
