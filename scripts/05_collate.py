#!/usr/bin/env python3
"""Turn LDSC logs into locked-panel tables and apply the h2 QC gate.

The gate decides which traits are powered enough to carry downstream. Excluded
alternative phenotype definitions are not parsed into production comparisons.

    python3 05_collate.py --mode h2 --logdir results/logs --out results/tables/h2_summary.tsv
    python3 05_collate.py --mode rg --logdir results/logs --out results/tables/rg_matrix.tsv
"""
import argparse
import glob
import itertools
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


def load_traits_meta(config_path="config/analysis_panel.tsv"):
    if not os.path.isfile(config_path):
        raise SystemExit(f"ERROR: locked panel manifest does not exist: {config_path}")
    metadata = pd.read_csv(config_path, sep="\t")
    if metadata["trait_id"].duplicated().any():
        dupes = metadata.loc[metadata["trait_id"].duplicated(), "trait_id"].tolist()
        raise ValueError(f"duplicate trait_id values in config: {dupes}")
    return metadata.set_index("trait_id")


def refuse_synthetic_in_real_directory(logdir):
    """Never let a copied smoke-test log turn into a plausible result table."""
    log_paths = glob.glob(os.path.join(logdir, "*.log"))

    def contains_marker(path):
        with open(path, encoding="utf-8") as handle:
            return SYNTHETIC_MARKER in handle.read()

    contains_synthetic = any(contains_marker(path) for path in log_paths)
    if contains_synthetic and "_smoketest" not in os.path.normpath(logdir).split(os.sep):
        raise SystemExit(
            "ERROR: synthetic LDSC logs were found outside results/_smoketest/. "
            "Refusing to collate fabricated values as real analysis output."
        )


def parse_h2(logdir, config_path="config/analysis_panel.tsv", expected_traits=None):
    refuse_synthetic_in_real_directory(logdir)
    meta_df = load_traits_meta(config_path)
    rows = []
    if expected_traits:
        if len(expected_traits) != len(set(expected_traits)):
            raise SystemExit("ERROR: duplicate trait in the requested h2 family")
        log_paths = [os.path.join(logdir, f"h2_{trait}.log") for trait in expected_traits]
        missing_logs = [path for path in log_paths if not os.path.isfile(path)]
        if missing_logs:
            raise SystemExit(f"ERROR: requested h2 logs are missing: {missing_logs}")
    else:
        log_paths = sorted(glob.glob(os.path.join(logdir, "h2_*.log")))
    for path in log_paths:
        trait = os.path.basename(path)[3:-4]
        if meta_df is not None and trait not in meta_df.index:
            raise SystemExit(f"ERROR: h2 log is outside the locked panel: {path}")
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
            # This reproduces CDG3's LDSC screening arithmetic only. It is not
            # MiXeR eligibility; that decision must come from a successful
            # univariate MiXeR fit and its own diagnostics.
            "ldsc_n_eff_h2": round(n_eff_h2, 1) if not np.isnan(n_eff_h2) else float("nan"),
            "ldsc_n_eff_h2_gt_12000": (n_eff_h2 > 12000) if not np.isnan(n_eff_h2) else False,
            "mixer_univariate_status": "NOT_RUN",
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


def parse_rg(logdir, config_path="config/analysis_panel.tsv", inclusion_path=None):
    """LDSC prints a fixed-width table after 'Summary of Genetic Correlation'."""
    refuse_synthetic_in_real_directory(logdir)
    expected = None
    if inclusion_path:
        inclusion = pd.read_csv(inclusion_path, sep="\t", dtype=str).fillna("")
        required = {"trait_id", "domain", "include_phase1"}
        missing = required.difference(inclusion.columns)
        if missing:
            raise SystemExit(f"ERROR: inclusion table missing columns: {sorted(missing)}")
        included = inclusion[inclusion["include_phase1"].str.lower().eq("true")]
        sleeps = included.loc[included["domain"].eq("sleep"), "trait_id"].tolist()
        diseases = included.loc[~included["domain"].eq("sleep"), "trait_id"].tolist()
        expected = set(itertools.product(sleeps, diseases))
        # Matrix logs contain all selected non-sleep traits for one sleep
        # trait. Controlled one-pair reruns use rg_SLEEP__DISEASE.log; never
        # mix those diagnostics into the readiness-selected matrix family.
        log_paths = [os.path.join(logdir, f"rg_{trait}.log") for trait in sleeps]
        missing_logs = [path for path in log_paths if not os.path.isfile(path)]
        if missing_logs:
            raise SystemExit(f"ERROR: readiness-selected rg logs are missing: {missing_logs}")
    else:
        log_paths = sorted(glob.glob(os.path.join(logdir, "rg_*.log")))
    frames = []
    for path in log_paths:
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
                # LDSC's fixed-width table does not quote p1/p2. Paths with
                # spaces (for example an external volume named "Extreme SSD")
                # therefore split into extra fields. Recover the two input
                # paths by their stable .sumstats.gz suffix, then split only
                # the numeric tail. Keep the old prose-footer behavior when
                # the row does not match this shape.
                match = re.match(
                    r"^(.*?)\.sumstats\.gz\s+(.*?)\.sumstats\.gz\s+(.+)$", l
                )
                if match:
                    numeric = match.group(3).split()
                    if len(numeric) == len(header) - 2:
                        fields = [
                            match.group(1) + ".sumstats.gz",
                            match.group(2) + ".sumstats.gz",
                            *numeric,
                        ]
                    else:
                        break
                else:
                    break
            body.append(fields)
        if body:
            frame = pd.DataFrame(body, columns=header)
            # Some LDSC versions round the fixed-width summary-table p column
            # to 0.0000 for small values even though the preceding scalar
            # ``P:`` lines retain scientific notation. LDSC emits one scalar
            # for every result row in the same order, so preserve the complete
            # sequence rather than turning finite p-values into zero/FDR zero.
            scalar_p = [line.split(":", 1)[1].strip() for line in lines[:i]
                        if line.startswith("P:")]
            if len(scalar_p) == len(frame):
                frame.loc[:, "p"] = scalar_p
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

    meta_df = load_traits_meta(config_path)
    if meta_df is not None:
        panel_ids = set(meta_df.index)
        observed_ids = set(df["sleep_trait"]).union(df["disease_trait"])
        unknown = observed_ids.difference(panel_ids)
        if unknown:
            raise SystemExit(f"ERROR: rg logs contain traits outside the locked panel: {sorted(unknown)}")
        wrong_sleep = [
            trait for trait in df["sleep_trait"].unique()
            if meta_df.loc[trait, "domain"] != "sleep"
        ]
        wrong_disease = [
            trait for trait in df["disease_trait"].unique()
            if meta_df.loc[trait, "domain"] == "sleep"
        ]
        if wrong_sleep or wrong_disease:
            raise SystemExit(
                f"ERROR: rg pair roles disagree with the panel: "
                f"sleep={wrong_sleep}, non_sleep={wrong_disease}"
            )
    pair_columns = ["sleep_trait", "disease_trait"]
    if df.duplicated(pair_columns).any():
        duplicates = df.loc[df.duplicated(pair_columns, keep=False), pair_columns]
        raise SystemExit(
            "ERROR: duplicate rg pairs found across logs: "
            + ", ".join(f"{a}__{b}" for a, b in duplicates.drop_duplicates().itertuples(index=False))
        )

    if expected is not None:
        actual = set(df[pair_columns].itertuples(index=False, name=None))
        if actual != expected:
            missing_pairs = sorted(expected.difference(actual))
            extra_pairs = sorted(actual.difference(expected))
            raise SystemExit(
                f"ERROR: rg logs do not match the readiness-selected family: "
                f"expected={len(expected)}, actual={len(actual)}, "
                f"missing={missing_pairs[:5]}, extra={extra_pairs[:5]}"
            )

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
    ap.add_argument("--config", default="config/analysis_panel.tsv")
    ap.add_argument("--inclusion", help="Phase 1 inclusion table defining the exact rg family")
    ap.add_argument("--traits", nargs="+", help="Exact trait family to collate in h2 mode")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    if a.mode == "h2" and a.inclusion:
        raise SystemExit("ERROR: --inclusion is only valid with --mode rg")
    if a.mode == "rg" and a.traits:
        raise SystemExit("ERROR: --traits is only valid with --mode h2")
    df = (
        parse_h2(a.logdir, a.config, a.traits)
        if a.mode == "h2"
        else parse_rg(a.logdir, a.config, a.inclusion)
    )
    if df.empty:
        print("!! nothing parsed - check that LDSC actually ran")
        raise SystemExit(1)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    print(f"Parsed {len(df)} {a.mode} result rows")
    print(f"wrote {a.out}")
