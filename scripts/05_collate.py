#!/usr/bin/env python3
"""Turn LDSC .log files into tidy tables, apply the h2 QC gate, and generate
dual-strategy problem comparison tables.

The gate decides which traits are powered enough to carry GenomicSEM downstream,
and the problem comparison table contrasts competing phenotype tail definitions,
insomnia summary stats releases, and liability scale ascertainment.

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
SYNTHETIC_MARKER = "SYNTHETIC SMOKE-TEST OUTPUT - NOT REAL LDSC RESULTS"


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_traits_meta(config_path="config/traits.tsv"):
    if os.path.exists(config_path):
        metadata = pd.read_csv(config_path, sep="\t")
        if metadata["trait_id"].duplicated().any():
            dupes = metadata.loc[metadata["trait_id"].duplicated(), "trait_id"].tolist()
            raise ValueError(f"duplicate trait_id values in config: {dupes}")
        return metadata.set_index("trait_id")
    return None


def refuse_synthetic_in_real_directory(logdir):
    """Never let a copied smoke-test log turn into a plausible result table."""
    log_paths = glob.glob(os.path.join(logdir, "*.log"))
    contains_synthetic = any(SYNTHETIC_MARKER in open(path, encoding="utf-8").read()
                             for path in log_paths)
    if contains_synthetic and "_smoketest" not in os.path.normpath(logdir).split(os.sep):
        raise SystemExit(
            "ERROR: synthetic LDSC logs were found outside results/_smoketest/. "
            "Refusing to collate fabricated values as real analysis output."
        )


def parse_h2(logdir, config_path="config/traits.tsv"):
    refuse_synthetic_in_real_directory(logdir)
    meta_df = load_traits_meta(config_path)
    rows = []
    for path in sorted(glob.glob(os.path.join(logdir, "h2_*.log"))):
        trait = os.path.basename(path)[3:-4]
        with open(path, encoding="utf-8") as handle:
            txt = handle.read()
        m = H2_PAT.search(txt)
        if not m:
            rows.append({"trait": trait, "input_log": path,
                         "note": "NO h2 IN LOG - check for errors"})
            continue
        h2, se = f(m.group(2)), f(m.group(3))
        i = INT_PAT.search(txt)
        z = h2 / se if se else float("nan")

        n_eff_h2 = float("nan")
        if meta_df is not None and trait in meta_df.index:
            row_meta = meta_df.loc[trait]
            n_tot = f(row_meta.get("n_total", float("nan")))
            n_case = f(row_meta.get("ncase", float("nan")))
            n_ctrl = f(row_meta.get("ncontrol", float("nan")))
            if not np.isnan(n_case) and not np.isnan(n_ctrl) and n_case > 0 and n_ctrl > 0:
                neff = 4.0 / (1.0 / n_case + 1.0 / n_ctrl)
            else:
                neff = n_tot
            if not np.isnan(neff) and not np.isnan(h2):
                n_eff_h2 = neff * h2

        rows.append({
            "trait": trait, "input_log": path,
            "scale": m.group(1).lower(),
            "h2": h2, "se": se, "z": round(z, 2),
            "intercept": f(i.group(1)) if i else float("nan"),
            "intercept_se": f(i.group(2)) if i else float("nan"),
            "lambda_gc": f(LAM_PAT.search(txt).group(1)) if LAM_PAT.search(txt) else float("nan"),
            "mean_chi2": f(CHI_PAT.search(txt).group(1)) if CHI_PAT.search(txt) else float("nan"),
            "ratio": f(RAT_PAT.search(txt).group(1)) if RAT_PAT.search(txt) else float("nan"),
            "n_eff_h2": round(n_eff_h2, 1) if not np.isnan(n_eff_h2) else float("nan"),
            "mixer_pass": (n_eff_h2 > 12000) if not np.isnan(n_eff_h2) else False
        })
    df = pd.DataFrame(rows)
    if "h2" not in df:
        return df
    df["pass_z"] = df["z"] >= Z_MIN
    df["pass_intercept"] = df["intercept"] <= INTERCEPT_MAX
    df["verdict"] = ["PASS" if (a and b) else "DROP"
                     for a, b in zip(df["pass_z"].fillna(False),
                                     df["pass_intercept"].fillna(False))]
    df["qc_reason"] = np.select(
        [df["z"].isna(), ~df["pass_z"], df["intercept"].isna(), ~df["pass_intercept"]],
        ["h2_or_se_missing", f"h2_z_below_{Z_MIN:g}", "intercept_missing", f"intercept_above_{INTERCEPT_MAX:g}"],
        default="pass",
    )
    return df.sort_values("z", ascending=False)


def generate_problem_comparison_table(df_h2, out_path):
    """Generate side-by-side comparison table for competing problem variants."""
    pairs = [
        ("Short Sleep", "shortsleep_dashti", "<7h (Dashti 2019)", "shortsleep_az", "<=5h (Austin-Zimmerman 2023)"),
        ("Long Sleep",  "longsleep_dashti",  ">=9h (Dashti 2019)", "longsleep_az",  ">=10h (Austin-Zimmerman 2023)"),
        ("Insomnia",    "insomnia_ukb",      "UKB-only (Public)",  "insomnia_full", "UKB+23andMe (Full)")
    ]

    comp_rows = []
    df_idx = df_h2.set_index("trait") if "trait" in df_h2 else df_h2

    for pheno_label, t1, opt1_name, t2, opt2_name in pairs:
        r1 = df_idx.loc[t1] if t1 in df_idx.index else None
        r2 = df_idx.loc[t2] if t2 in df_idx.index else None

        h2_1 = f"{r1['h2']:.4f} ({r1['se']:.4f})" if r1 is not None and 'h2' in r1 else "NA"
        z1 = f"{r1['z']:.2f}" if r1 is not None and 'z' in r1 else "NA"
        icept1 = f"{r1['intercept']:.4f}" if r1 is not None and 'intercept' in r1 else "NA"
        v1 = r1['verdict'] if r1 is not None and 'verdict' in r1 else "NA"
        mixer1 = "PASS" if r1 is not None and r1.get('mixer_pass', False) else "DROP"

        h2_2 = f"{r2['h2']:.4f} ({r2['se']:.4f})" if r2 is not None and 'h2' in r2 else "NA"
        z2 = f"{r2['z']:.2f}" if r2 is not None and 'z' in r2 else "NA"
        icept2 = f"{r2['intercept']:.4f}" if r2 is not None and 'intercept' in r2 else "NA"
        v2 = r2['verdict'] if r2 is not None and 'verdict' in r2 else "NA"
        mixer2 = "PASS" if r2 is not None and r2.get('mixer_pass', False) else "DROP"

        comp_rows.append({
            "Phenotype": pheno_label,
            "Option 1": opt1_name,
            "Opt1_h2": h2_1, "Opt1_Z": z1, "Opt1_Intercept": icept1, "Opt1_QC": v1, "Opt1_MiXeR": mixer1,
            "Option 2": opt2_name,
            "Opt2_h2": h2_2, "Opt2_Z": z2, "Opt2_Intercept": icept2, "Opt2_QC": v2, "Opt2_MiXeR": mixer2,
        })

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nwrote problem comparison summary table: {out_path}")
    print(comp_df.to_string(index=False))
    return comp_df


def parse_rg(logdir):
    """LDSC prints a fixed-width table after 'Summary of Genetic Correlation'."""
    refuse_synthetic_in_real_directory(logdir)
    frames = []
    for path in sorted(glob.glob(os.path.join(logdir, "rg_*.log"))):
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        try:
            i = next(k for k, l in enumerate(lines)
                     if "Summary of Genetic Correlation Results" in l)
        except StopIteration:
            print(f"  !! no results table in {path}")
            continue
        header = lines[i + 1].split()
        if not {"p1", "p2", "rg", "se", "p"}.issubset(header):
            print(f"  !! unexpected rg header in {path}: {header}")
            continue
        body = []
        for l in lines[i + 2:]:
            if not l.strip():
                break
            fields = l.split()
            if len(fields) != len(header):
                # LDSC writes a prose footer after the table in some versions.
                break
            body.append(fields)
        if body:
            frame = pd.DataFrame(body, columns=header)
            frame["input_log"] = path
            frames.append(frame)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    for c in ["rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "gcov_int"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ["p1", "p2"]:
        df[c] = (df[c].str.replace(r".*/", "", regex=True)
                      .str.replace(".sumstats.gz", "", regex=False))
    df = df.rename(columns={"p1": "sleep_trait", "p2": "disease_trait"})

    d = df.dropna(subset=["p"]).sort_values("p").copy()
    m = len(d)
    if m == 0:
        return pd.DataFrame()
    raw = d["p"].values * m / np.arange(1, m + 1)
    d["fdr"] = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
    df["fdr"] = np.nan
    df.loc[d.index, "fdr"] = d["fdr"]
    return df.sort_values("p")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["h2", "rg"], required=True)
    ap.add_argument("--logdir", default="results/logs")
    ap.add_argument("--config", default="config/traits.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    df = parse_h2(a.logdir, a.config) if a.mode == "h2" else parse_rg(a.logdir)
    if df.empty:
        print("!! nothing parsed - check that LDSC actually ran")
        raise SystemExit(1)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    print(f"Parsed {len(df)} {a.mode} result rows")
    print(f"wrote {a.out}")

    if a.mode == "h2":
        comp_out = os.path.join(os.path.dirname(a.out), "problem_comparison_summary.tsv")
        generate_problem_comparison_table(df, comp_out)
