#!/usr/bin/env python3
"""Build the quantitative deep-analysis layer from the immutable Phase-1 atlas."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata


PRIMARY = "PRIMARY_PHASE1"
SENSITIVITY = "QC_FAILED_SENSITIVITY"
LOCKED_FDR = 0.05
PAIR_FIELDS = ["sleep_trait", "disease_trait"]
SHA256 = re.compile(r"[0-9a-f]{64}$")
SNP_OVERLAP = re.compile(
    r"Computing rg for phenotype \d+/\d+\s+"
    r"Reading summary statistics from data/munged/([a-z0-9_]+)\.sumstats\.gz "
    r".*?(\d+) SNPs with valid alleles\.",
    re.S,
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def frame_text(frame: pd.DataFrame) -> str:
    return frame.to_csv(
        sep="\t", index=False, line_terminator="\n", na_rep="NA", float_format="%.15g"
    )


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def bh(values: np.ndarray) -> np.ndarray:
    if values.ndim != 1 or not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        fail("BH input must be a finite probability vector")
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(values) / np.arange(1, len(values) + 1))[::-1]
    )[::-1].clip(max=1.0)
    result = np.empty(len(values), dtype=float)
    result[order] = adjusted
    return result


def bh_dict(values: dict[str, float]) -> dict[str, float]:
    keys = list(values)
    adjusted = bh(np.array([values[key] for key in keys], dtype=float))
    return dict(zip(keys, adjusted))


def finite_column(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        fail(f"{column} contains a non-finite value")
    return values


def entropy(values: list[str]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    counts = np.array(list(Counter(values).values()), dtype=float)
    proportions = counts / counts.sum()
    raw = float(-(proportions * np.log(proportions)).sum())
    maximum = math.log(len(counts)) if len(counts) > 1 else 0.0
    return raw, raw / maximum if maximum else 0.0


def label_join(rows: pd.DataFrame, value: str, label: str) -> str:
    if rows.empty:
        return "NA"
    return f"{rows.iloc[0][value]}__{rows.iloc[0][label]}"


def parse_overlap_logs(root: Path, sleep_traits: list[str]) -> tuple[dict[tuple[str, str], int], list[Path]]:
    overlap: dict[tuple[str, str], int] = {}
    paths: list[Path] = []
    families = [root / "results/logs", root / "results/logs/rg_qc_failed_sensitivity"]
    for directory in families:
        for sleep in sleep_traits:
            path = directory / f"rg_{sleep}.log"
            if not path.is_file():
                fail(f"missing LDSC pair-family log: {path}")
            paths.append(path)
            for external, count in SNP_OVERLAP.findall(path.read_text(encoding="utf-8")):
                key = (sleep, external)
                if key in overlap:
                    fail(f"duplicate SNP-overlap evidence for {sleep} x {external}")
                overlap[key] = int(count)
    return overlap, paths


def validate_inputs(root: Path) -> tuple[dict[str, object], dict[str, Path], list[Path]]:
    paths = {
        "rg_matrix": root / "results/tables/rg_matrix.tsv",
        "h2_summary": root / "results/tables/h2_summary.tsv",
        "trait_readiness": root / "results/tables/trait_readiness.tsv",
        "panel": root / "config/analysis_panel.tsv",
        "panel_lock": root / "config/analysis_panel.lock.json",
        "public_sources": root / "config/public_gwas_sources.tsv",
        "substitutions": root / "config/public_gwas_substitution_candidates.tsv",
        "measurement_modes": root / "config/sleep_measurement_modes.tsv",
        "core_provenance": root / "results/atlas/core.provenance.json",
        "atlas_traits": root / "results/atlas/traits.tsv",
        "atlas_pairs": root / "results/atlas/trait_pairs.tsv",
        "legacy_phase1_report": root / "results/tables/phase1_summary_report.md",
        "requested_phase1_completion": root / "results/tables/phase1_completion.json",
    }
    for name, path in paths.items():
        if name != "requested_phase1_completion" and not path.is_file():
            fail(f"required Phase-1 evidence is absent: {path}")

    panel = read_tsv(paths["panel"])
    rg = read_tsv(paths["rg_matrix"])
    h2 = read_tsv(paths["h2_summary"])
    readiness = read_tsv(paths["trait_readiness"])
    modes = read_tsv(paths["measurement_modes"])
    sleep = panel.loc[panel["domain"].eq("sleep"), "trait_id"].tolist()
    external = panel.loc[~panel["domain"].eq("sleep"), "trait_id"].tolist()
    expected = {(a, b) for a in sleep for b in external}
    observed = set(rg[PAIR_FIELDS].itertuples(index=False, name=None))
    tiers = rg["analysis_tier"].value_counts().to_dict()
    if (
        len(panel) != 45 or len(sleep) != 12 or len(external) != 33
        or panel["trait_id"].duplicated().any()
        or len(rg) != 396 or rg.duplicated(PAIR_FIELDS).any() or observed != expected
        or tiers != {PRIMARY: 372, SENSITIVITY: 24}
        or set(rg.loc[rg["analysis_tier"].eq(SENSITIVITY), "disease_trait"]) != {"t2d", "melanoma"}
        or len(h2) != 45 or set(h2["trait"]) != set(panel["trait_id"])
        or len(readiness) != 45 or set(readiness["trait_id"]) != set(panel["trait_id"])
        or len(modes) != 12 or set(modes["trait_id"]) != set(sleep)
    ):
        fail("canonical Phase-1 panel/pair/QC family differs from the locked 45/12/33/396/372/24 contract")
    if set(rg.loc[rg["analysis_tier"].eq(SENSITIVITY), "interpretation_status"]) != {
        "EXCLUDED_FROM_PRIMARY_INFERENCE"
    }:
        fail("QC-failed rows are not cleanly segregated from primary inference")

    numeric = ["rg", "se", "z", "p", "fdr", "gcov_int", "gcov_int_se"]
    for column in numeric:
        values = finite_column(rg, column)
        if column in {"p", "fdr"} and ((values < 0) | (values > 1)).any():
            fail(f"{column} is outside [0,1]")
    if (pd.to_numeric(rg["se"]) <= 0).any() or (pd.to_numeric(rg["gcov_int_se"]) <= 0).any():
        fail("pair standard errors must be positive")

    p_values = pd.to_numeric(rg["p"]).to_numpy(dtype=float)
    recomputed_all = bh(p_values)
    observed_all = pd.to_numeric(rg["fdr"]).to_numpy(dtype=float)
    if not np.allclose(recomputed_all, observed_all, rtol=1e-12, atol=1e-15):
        fail("reported locked-396 BH values do not reproduce from P values")
    primary_mask = rg["analysis_tier"].eq(PRIMARY).to_numpy()
    recomputed_primary = bh(p_values[primary_mask])
    observed_primary = pd.to_numeric(
        rg.loc[primary_mask, "fdr_primary_phase1"], errors="coerce"
    ).to_numpy(dtype=float)
    if not np.allclose(recomputed_primary, observed_primary, rtol=1e-12, atol=1e-15):
        fail("reported primary-372 BH sensitivity values do not reproduce")
    locked_primary_significant = int(((recomputed_all <= LOCKED_FDR) & primary_mask).sum())
    primary_only_significant = int((recomputed_primary <= LOCKED_FDR).sum())
    if locked_primary_significant != 153 or primary_only_significant != 155:
        fail("Phase-1 discovery counts do not reproduce as locked 153 / sensitivity 155")

    recomputed_z = pd.to_numeric(rg["rg"]) / pd.to_numeric(rg["se"])
    z_difference = (recomputed_z - pd.to_numeric(rg["z"])).abs()
    if z_difference.max() > 0.05:
        fail("rg/SE Z differs from the rounded reported LDSC Z by more than 0.05")

    overlap, log_paths = parse_overlap_logs(root, sleep)
    if set(overlap) != expected or min(overlap.values()) <= 0:
        fail("SNP-overlap evidence does not cover the exact 396-pair family")

    core = json.loads(paths["core_provenance"].read_text(encoding="utf-8"))
    if (
        core.get("rg_matrix_sha256") != sha256(paths["rg_matrix"])
        or core.get("h2_summary_sha256") != sha256(paths["h2_summary"])
        or core.get("panel_sha256") != sha256(paths["panel"])
        or core.get("outputs", {}).get("results/atlas/traits.tsv") != sha256(paths["atlas_traits"])
        or core.get("outputs", {}).get("results/atlas/trait_pairs.tsv") != sha256(paths["atlas_pairs"])
    ):
        fail("atlas core provenance hash chain is invalid")
    diagnostics = {
        "schema_version": "sleep-atlas-phase1-deep-analysis-inputs.1",
        "verification_status": "PASS_WITH_DOCUMENTED_METADATA_DISCREPANCIES",
        "trait_count": 45,
        "sleep_trait_count": 12,
        "external_trait_count": 33,
        "pair_count": 396,
        "unique_pair_count": len(observed),
        "primary_pair_count": int(primary_mask.sum()),
        "sensitivity_pair_count": int((~primary_mask).sum()),
        "locked_396_primary_fdr_le_0_05": locked_primary_significant,
        "primary_372_sensitivity_fdr_le_0_05": primary_only_significant,
        "all_396_fdr_le_0_05_including_sensitivity": int((recomputed_all <= LOCKED_FDR).sum()),
        "maximum_abs_recomputed_vs_reported_z_difference": float(z_difference.max()),
        "minimum_snp_overlap": min(overlap.values()),
        "maximum_snp_overlap": max(overlap.values()),
        "sensitivity_external_traits": ["t2d", "melanoma"],
        "requested_phase1_completion": {
            "path": "results/tables/phase1_completion.json",
            "status": "ABSENT",
            "diagnosis": "No file by this name exists in the worktree or Git history; results/atlas/core.provenance.json is the checksum-bound canonical completion equivalent.",
        },
        "core_provenance_count_semantics": {
            "reported_field": "primary_fdr_significant_pair_count",
            "reported_value": core.get("primary_fdr_significant_pair_count"),
            "diagnosis": "The field counts the separate primary-only 372-test BH family (155), not the locked all-396 primary discovery family (153). The immutable core is not modified; deep analysis uses rg_matrix.fdr.",
        },
        "legacy_phase1_summary": {
            "status": "STALE_NOT_USED_FOR_ANALYSIS",
            "diagnosis": "The legacy Markdown predates the complete rerun and reports 27 h2 traits / 136 rg pairs; canonical TSVs and core provenance supersede it.",
        },
        "primary_definition": "analysis_tier == PRIMARY_PHASE1 and locked all-396 BH fdr <= 0.05",
    }
    return diagnostics, paths, log_paths


def build_master(root: Path, paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    panel = read_tsv(paths["panel"])
    rg = read_tsv(paths["rg_matrix"])
    h2 = read_tsv(paths["h2_summary"])
    readiness = read_tsv(paths["trait_readiness"])
    modes = read_tsv(paths["measurement_modes"])
    substitutions = read_tsv(paths["substitutions"])
    overlap, _ = parse_overlap_logs(
        root, panel.loc[panel["domain"].eq("sleep"), "trait_id"].tolist()
    )
    panel_by = panel.set_index("trait_id").to_dict("index")
    h2_by = h2.set_index("trait").to_dict("index")
    ready_by = readiness.set_index("trait_id").to_dict("index")
    mode_by = modes.set_index("trait_id").to_dict("index")
    substitution_by = substitutions.set_index("trait_id").to_dict("index")
    rows = []
    for record in rg.to_dict("records"):
        sleep, external = record["sleep_trait"], record["disease_trait"]
        s, e = panel_by[sleep], panel_by[external]
        sh, eh = h2_by[sleep], h2_by[external]
        sr, er = ready_by[sleep], ready_by[external]
        mode = mode_by[sleep]
        source_caveat = substitution_by.get(external, {}).get("phenotype_caveat", "")
        qc_notes = []
        if sh["qc_reason"] != "pass":
            qc_notes.append(f"sleep_h2:{sh['qc_reason']}")
        if eh["qc_reason"] != "pass":
            qc_notes.append(f"external_h2:{eh['qc_reason']}")
        for prefix, item in (("sleep", sr), ("external", er)):
            for field in ("source_issues", "harmonization_issues", "ldsc_issues", "phase1_issues"):
                if item.get(field, "") not in {"", "none", "NA"}:
                    qc_notes.append(f"{prefix}_{field}:{item[field]}")
        if source_caveat:
            qc_notes.append(f"source_substitution:{source_caveat}")
        sample_size = e["n_total"] if e["n_total"] not in {"", "NA"} else "NA"
        rg_value, se_value = float(record["rg"]), float(record["se"])
        rows.append({
            "sleep_trait": sleep,
            "sleep_trait_label": s["label"],
            "external_trait": external,
            "external_trait_label": e["label"],
            "external_domain": e["domain"],
            "rg": rg_value,
            "se": se_value,
            "z": rg_value / se_value,
            "reported_z": float(record["z"]),
            "z_abs_difference": abs(rg_value / se_value - float(record["z"])),
            "p": float(record["p"]),
            "fdr": float(record["fdr"]),
            "fdr_primary_372_sensitivity": (
                float(record["fdr_primary_phase1"])
                if record["fdr_primary_phase1"] not in {"", "NA"} else np.nan
            ),
            "locked_primary_significant": (
                record["analysis_tier"] == PRIMARY and float(record["fdr"]) <= LOCKED_FDR
            ),
            "abs_rg": abs(rg_value),
            "direction": "POSITIVE" if rg_value > 0 else ("NEGATIVE" if rg_value < 0 else "ZERO"),
            "gcov_intercept": float(record["gcov_int"]),
            "gcov_intercept_se": float(record["gcov_int_se"]),
            "gcov_intercept_z": float(record["gcov_int"]) / float(record["gcov_int_se"]),
            "SNP_overlap": overlap[(sleep, external)],
            "sleep_h2": float(sh["h2"]),
            "sleep_h2_se": float(sh["se"]),
            "sleep_h2_z": float(sh["z"]),
            "sleep_h2_scale": sh["scale"],
            "external_h2": float(eh["h2"]),
            "external_h2_se": float(eh["se"]),
            "external_h2_z": float(eh["z"]),
            "external_h2_scale": eh["scale"],
            "sleep_h2_intercept": float(sh["intercept"]),
            "sleep_h2_intercept_se": float(sh["intercept_se"]),
            "sleep_h2_attenuation_ratio": float(sh["ratio"]) if sh["ratio"] not in {"", "NA"} else np.nan,
            "external_h2_intercept": float(eh["intercept"]),
            "external_h2_intercept_se": float(eh["intercept_se"]),
            "external_h2_attenuation_ratio": float(eh["ratio"]) if eh["ratio"] not in {"", "NA"} else np.nan,
            "primary_or_sensitivity": record["analysis_tier"],
            "interpretation_status": record["interpretation_status"],
            "QC_notes": "; ".join(qc_notes) if qc_notes else "none",
            "phenotype_source": e["source_note"],
            "external_source_id": e["source_id"],
            "external_dataset_version": e["dataset_version"],
            "ancestry": e["ancestry"],
            "sample_size": sample_size,
            "sleep_sample_size": s["n_total"] if s["n_total"] not in {"", "NA"} else "NA",
            "sleep_measurement_category": mode["measurement_category"],
            "sleep_collection_mode": mode["collection_mode"],
            "sleep_phenotype_family": mode["phenotype_family"],
            "sleep_higher_value_interpretation": mode["higher_value_interpretation"],
            "source_substitution": external in substitution_by,
            "source_substitution_caveat": source_caveat or "NA",
            "input_log": record["input_log"],
        })
    master = pd.DataFrame(rows)
    sleep_order = {value: index for index, value in enumerate(panel.loc[panel["domain"].eq("sleep"), "trait_id"])}
    external_order = {value: index for index, value in enumerate(panel.loc[~panel["domain"].eq("sleep"), "trait_id"])}
    master["_sleep_order"] = master["sleep_trait"].map(sleep_order)
    master["_external_order"] = master["external_trait"].map(external_order)
    master.sort_values(["_sleep_order", "_external_order"], inplace=True)
    master.drop(columns=["_sleep_order", "_external_order"], inplace=True)
    return master, panel, h2, modes


def global_summary(master: pd.DataFrame) -> tuple[dict[str, object], str]:
    groups = {
        "PRIMARY": master.loc[master["primary_or_sensitivity"].eq(PRIMARY)],
        "QC_FAILED_SENSITIVITY": master.loc[master["primary_or_sensitivity"].eq(SENSITIVITY)],
    }
    payload: dict[str, object] = {
        "schema_version": "sleep-atlas-phase1-global-summary.1",
        "primary_significance_definition": "PRIMARY_PHASE1 and BH FDR <= 0.05 in the locked 396-pair family",
        "sensitivity_rows_excluded_from_headline_counts": True,
        "strata": {},
    }
    for name, frame in groups.items():
        significant = frame.loc[frame["fdr"].le(LOCKED_FDR)]
        by_z = frame.iloc[frame["z"].abs().argmax()]
        by_fdr = frame.sort_values(["fdr", "p"]).iloc[0]
        by_abs = frame.sort_values("abs_rg", ascending=False).iloc[0]
        positive = frame.sort_values("rg", ascending=False).iloc[0]
        negative = frame.sort_values("rg").iloc[0]
        payload["strata"][name] = {
            "pair_count": len(frame),
            "positive_rg_count": int(frame["rg"].gt(0).sum()),
            "negative_rg_count": int(frame["rg"].lt(0).sum()),
            "fdr_significant_count": len(significant),
            "fdr_significant_positive_count": int(significant["rg"].gt(0).sum()),
            "fdr_significant_negative_count": int(significant["rg"].lt(0).sum()),
            "median_rg": float(frame["rg"].median()),
            "median_abs_rg": float(frame["abs_rg"].median()),
            "rg_quantiles": {str(q): float(frame["rg"].quantile(q)) for q in (0, .05, .25, .5, .75, .95, 1)},
            "abs_rg_quantiles": {str(q): float(frame["abs_rg"].quantile(q)) for q in (0, .05, .25, .5, .75, .95, 1)},
            "largest_positive": {"pair": f"{positive.sleep_trait}__{positive.external_trait}", "rg": float(positive.rg)},
            "largest_negative": {"pair": f"{negative.sleep_trait}__{negative.external_trait}", "rg": float(negative.rg)},
            "strongest_by_abs_z": {"pair": f"{by_z.sleep_trait}__{by_z.external_trait}", "z": float(by_z.z)},
            "strongest_by_fdr": {"pair": f"{by_fdr.sleep_trait}__{by_fdr.external_trait}", "fdr": float(by_fdr.fdr)},
            "strongest_by_abs_rg": {"pair": f"{by_abs.sleep_trait}__{by_abs.external_trait}", "abs_rg": float(by_abs.abs_rg)},
        }
    primary = payload["strata"]["PRIMARY"]
    lines = [
        "# Global Phase-1 atlas summary", "",
        "Primary inference uses the locked 396-test BH family. The 24 T2D/melanoma rows are reported separately and never enter headline counts.", "",
        f"- Primary pairs: **{primary['pair_count']}**",
        f"- Primary FDR-significant connections: **{primary['fdr_significant_count']}**",
        f"- Significant positive / negative: **{primary['fdr_significant_positive_count']} / {primary['fdr_significant_negative_count']}**",
        f"- All primary positive / negative estimates: **{primary['positive_rg_count']} / {primary['negative_rg_count']}**",
        f"- Median rg / median |rg|: **{primary['median_rg']:.4f} / {primary['median_abs_rg']:.4f}**",
        f"- Largest positive: **{primary['largest_positive']['pair']}**, rg={primary['largest_positive']['rg']:.4f}",
        f"- Largest negative: **{primary['largest_negative']['pair']}**, rg={primary['largest_negative']['rg']:.4f}",
        f"- Strongest by |Z|: **{primary['strongest_by_abs_z']['pair']}**, Z={primary['strongest_by_abs_z']['z']:.3f}",
        "", "## QC-failed sensitivity stratum", "",
        f"The sensitivity stratum contains {payload['strata']['QC_FAILED_SENSITIVITY']['pair_count']} pairs and is not pooled with primary inference.", "",
    ]
    return payload, "\n".join(lines)


def sleep_hubs(master: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)]
    rows = []
    for sleep, frame in primary.groupby("sleep_trait", sort=False):
        significant = frame.loc[frame["fdr"].le(LOCKED_FDR)]
        positive = significant.loc[significant["rg"].gt(0)].sort_values("rg", ascending=False)
        negative = significant.loc[significant["rg"].lt(0)].sort_values("rg")
        raw_entropy, normalized_entropy = entropy(significant["external_domain"].tolist())
        rows.append({
            "sleep_trait": sleep,
            "sleep_trait_label": frame.iloc[0]["sleep_trait_label"],
            "measurement_category": frame.iloc[0]["sleep_measurement_category"],
            "tested_external_traits": len(frame),
            "fdr_significant_connections": len(significant),
            "percent_tested_connected": 100 * len(significant) / len(frame),
            "significant_positive": int(significant["rg"].gt(0).sum()),
            "significant_negative": int(significant["rg"].lt(0).sum()),
            "sum_abs_rg_significant": significant["abs_rg"].sum(),
            "mean_abs_rg_significant": significant["abs_rg"].mean() if len(significant) else np.nan,
            "median_abs_rg_significant": significant["abs_rg"].median() if len(significant) else np.nan,
            "maximum_abs_rg_significant": significant["abs_rg"].max() if len(significant) else np.nan,
            "strongest_positive_connection": label_join(positive, "external_trait", "external_trait_label"),
            "strongest_positive_rg": positive.iloc[0]["rg"] if len(positive) else np.nan,
            "strongest_negative_connection": label_join(negative, "external_trait", "external_trait_label"),
            "strongest_negative_rg": negative.iloc[0]["rg"] if len(negative) else np.nan,
            "connected_domain_count": significant["external_domain"].nunique(),
            "connected_domains": ";".join(sorted(significant["external_domain"].unique())),
            "domain_entropy": raw_entropy,
            "normalized_domain_entropy": normalized_entropy,
            "mean_external_h2_connected": significant["external_h2"].mean() if len(significant) else np.nan,
            "connections_fdr_le_0_01": int(significant["fdr"].le(.01).sum()),
            "connections_fdr_le_0_001": int(significant["fdr"].le(.001).sum()),
            "proportion_significant_surviving_fdr_0_01": significant["fdr"].le(.01).mean() if len(significant) else np.nan,
            "proportion_significant_surviving_fdr_0_001": significant["fdr"].le(.001).mean() if len(significant) else np.nan,
        })
    result = pd.DataFrame(rows)
    rank_columns = {
        "rank_by_connection_count": ("fdr_significant_connections", False),
        "rank_by_sum_abs_rg": ("sum_abs_rg_significant", False),
        "rank_by_mean_abs_rg": ("mean_abs_rg_significant", False),
        "rank_by_domain_diversity": ("normalized_domain_entropy", False),
        "rank_by_negative_connections": ("significant_negative", False),
    }
    for output, (column, ascending) in rank_columns.items():
        result[output] = result[column].rank(method="min", ascending=ascending).astype("Int64")
    order = {value: index for index, value in enumerate(panel.loc[panel["domain"].eq("sleep"), "trait_id"])}
    result["_order"] = result["sleep_trait"].map(order)
    return result.sort_values("_order").drop(columns="_order")


def external_hubs(master: pd.DataFrame, modes: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)]
    objective = set(modes.loc[modes["measurement_category"].eq("OBJECTIVE_ACTIGRAPHY"), "trait_id"])
    subjective = set(modes.loc[modes["measurement_category"].eq("SUBJECTIVE_SELF_REPORT"), "trait_id"])
    rows = []
    for external, frame in primary.groupby("external_trait", sort=False):
        significant = frame.loc[frame["fdr"].le(LOCKED_FDR)]
        strongest = frame.sort_values("abs_rg", ascending=False).iloc[0]
        rows.append({
            "external_trait": external,
            "external_trait_label": frame.iloc[0]["external_trait_label"],
            "external_domain": frame.iloc[0]["external_domain"],
            "significant_sleep_connections": len(significant),
            "significant_positive": int(significant["rg"].gt(0).sum()),
            "significant_negative": int(significant["rg"].lt(0).sum()),
            "average_rg_all_sleep_traits": frame["rg"].mean(),
            "average_abs_rg_all_sleep_traits": frame["abs_rg"].mean(),
            "maximum_abs_rg": frame["abs_rg"].max(),
            "minimum_rg": frame["rg"].min(),
            "maximum_rg": frame["rg"].max(),
            "most_strongly_connected_sleep_trait": strongest["sleep_trait"],
            "most_strongly_connected_sleep_trait_label": strongest["sleep_trait_label"],
            "strongest_rg": strongest["rg"],
            "objective_sleep_phenotypes_connected": int(significant["sleep_trait"].isin(objective).sum()),
            "subjective_sleep_phenotypes_connected": int(significant["sleep_trait"].isin(subjective).sum()),
            "connected_to_at_least_3_sleep_traits": len(significant) >= 3,
            "connected_to_at_least_6_sleep_traits": len(significant) >= 6,
            "connected_to_at_least_9_sleep_traits": len(significant) >= 9,
        })
    result = pd.DataFrame(rows)
    result["hub_rank"] = result["significant_sleep_connections"].rank(method="min", ascending=False).astype("Int64")
    external_order = {value: index for index, value in enumerate(panel.loc[~panel["domain"].eq("sleep"), "trait_id"])}
    result["_order"] = result["external_trait"].map(external_order)
    return result.sort_values("_order").drop(columns="_order")


def domain_architecture(master: pd.DataFrame, permutations: int, seed: int) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)].copy()
    observed_rows = []
    for domain, frame in primary.groupby("external_domain", sort=True):
        significant = frame.loc[frame["fdr"].le(LOCKED_FDR)]
        strongest = frame.sort_values("abs_rg", ascending=False).iloc[0]
        by_sleep = significant.groupby("sleep_trait").size() if len(significant) else pd.Series(dtype=int)
        max_count = int(by_sleep.max()) if len(by_sleep) else 0
        top_sleep = ";".join(sorted(by_sleep.loc[by_sleep.eq(max_count)].index)) if max_count else "NA"
        observed_rows.append({
            "external_domain": domain,
            "external_trait_count": frame["external_trait"].nunique(),
            "tested_pairs": len(frame),
            "fdr_significant_pairs": len(significant),
            "significant_proportion": len(significant) / len(frame),
            "significant_positive": int(significant["rg"].gt(0).sum()),
            "significant_negative": int(significant["rg"].lt(0).sum()),
            "mean_rg": frame["rg"].mean(),
            "mean_abs_rg": frame["abs_rg"].mean(),
            "median_abs_rg": frame["abs_rg"].median(),
            "strongest_pair": f"{strongest.sleep_trait}__{strongest.external_trait}",
            "strongest_pair_rg": strongest["rg"],
            "sleep_trait_with_highest_domain_connectivity": top_sleep,
            "highest_sleep_connectivity_count": max_count,
        })
    result = pd.DataFrame(observed_rows)

    trait_scores = primary.groupby("external_trait")["locked_primary_significant"].sum().astype(int)
    trait_domains = primary.drop_duplicates("external_trait").set_index("external_trait")["external_domain"]
    traits = trait_scores.index.to_numpy()
    labels = trait_domains.loc[traits].to_numpy()
    scores = trait_scores.loc[traits].to_numpy()
    domains = sorted(set(labels))
    observed = {domain: int(scores[labels == domain].sum()) for domain in domains}
    random = np.random.default_rng(seed)
    null = {domain: np.empty(permutations, dtype=int) for domain in domains}
    for iteration in range(permutations):
        shuffled = random.permutation(labels)
        for domain in domains:
            null[domain][iteration] = int(scores[shuffled == domain].sum())
    enrichment_p = {
        domain: (1 + int((null[domain] >= observed[domain]).sum())) / (permutations + 1)
        for domain in domains
    }
    depletion_p = {
        domain: (1 + int((null[domain] <= observed[domain]).sum())) / (permutations + 1)
        for domain in domains
    }
    enrichment_q = bh_dict(enrichment_p)
    for index, row in result.iterrows():
        domain = row["external_domain"]
        result.loc[index, "permutation_expected_significant_pairs"] = null[domain].mean()
        result.loc[index, "permutation_sd_significant_pairs"] = null[domain].std(ddof=1)
        result.loc[index, "permutation_enrichment_p"] = enrichment_p[domain]
        result.loc[index, "permutation_enrichment_bh_q"] = enrichment_q[domain]
        result.loc[index, "permutation_depletion_p"] = depletion_p[domain]
        result.loc[index, "permutation_count"] = permutations
        result.loc[index, "permutation_unit"] = "external-trait 12-rg profile"
        result.loc[index, "dependence_caveat"] = (
            "Domain labels were permuted across whole external-trait profiles, preserving all 12 correlated sleep edges within a trait; residual dependence among external traits is not modeled."
        )
    return result


def measurement_comparison(master: pd.DataFrame, modes: pd.DataFrame) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)].copy()
    categories = modes.groupby("measurement_category")["trait_id"].apply(set).to_dict()
    subjective = categories["SUBJECTIVE_SELF_REPORT"]
    objective = categories["OBJECTIVE_ACTIGRAPHY"]
    rows = []
    for external, frame in primary.groupby("external_trait", sort=False):
        sub = frame.loc[frame["sleep_trait"].isin(subjective)]
        obj = frame.loc[frame["sleep_trait"].isin(objective)]
        sub_sig = sub.loc[sub["fdr"].le(LOCKED_FDR)]
        obj_sig = obj.loc[obj["fdr"].le(LOCKED_FDR)]
        if len(sub_sig) and len(obj_sig):
            connection_class = "BOTH_SUBJECTIVE_AND_OBJECTIVE"
        elif len(sub_sig):
            connection_class = "SUBJECTIVE_ONLY"
        elif len(obj_sig):
            connection_class = "OBJECTIVE_ONLY"
        else:
            connection_class = "NEITHER"
        opposite_pairs = []
        for left in sub_sig.itertuples():
            for right in obj_sig.itertuples():
                if np.sign(left.rg) != np.sign(right.rg):
                    opposite_pairs.append(f"{left.sleep_trait}:{left.rg:.4g}|{right.sleep_trait}:{right.rg:.4g}")
        difference = obj["abs_rg"].mean() - sub["abs_rg"].mean()
        rows.append({
            "external_trait": external,
            "external_trait_label": frame.iloc[0]["external_trait_label"],
            "external_domain": frame.iloc[0]["external_domain"],
            "subjective_trait_count": len(sub),
            "objective_trait_count": len(obj),
            "subjective_significant_count": len(sub_sig),
            "objective_significant_count": len(obj_sig),
            "connection_class": connection_class,
            "subjective_mean_rg": sub["rg"].mean(),
            "objective_mean_rg": obj["rg"].mean(),
            "subjective_mean_abs_rg": sub["abs_rg"].mean(),
            "objective_mean_abs_rg": obj["abs_rg"].mean(),
            "objective_minus_subjective_mean_abs_rg": difference,
            "group_mean_sign_reversal": np.sign(sub["rg"].mean()) != np.sign(obj["rg"].mean()),
            "significant_pair_direction_discordance": bool(opposite_pairs),
            "discordant_significant_pairs": ";".join(opposite_pairs) if opposite_pairs else "NA",
            "objective_materially_larger": difference >= .10,
            "subjective_materially_larger": difference <= -.10,
            "material_difference_definition": "absolute difference in group mean |rg| >= 0.10; descriptive, not an independence-based test",
        })
    return pd.DataFrame(rows)


def similarity_and_clustering(master: pd.DataFrame, panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)]
    sleep_order = panel.loc[panel["domain"].eq("sleep"), "trait_id"].tolist()
    external_order = panel.loc[
        (~panel["domain"].eq("sleep")) & (~panel["trait_id"].isin(["t2d", "melanoma"])), "trait_id"
    ].tolist()
    matrix = primary.pivot(index="sleep_trait", columns="external_trait", values="rg").loc[sleep_order, external_order]
    significant = primary.assign(sig=primary["fdr"].le(LOCKED_FDR)).pivot(
        index="sleep_trait", columns="external_trait", values="sig"
    ).loc[sleep_order, external_order]
    values = matrix.to_numpy(dtype=float)
    rows = []
    pearson_matrix = np.corrcoef(values)
    for i, left in enumerate(sleep_order):
        for j in range(i + 1, len(sleep_order)):
            right = sleep_order[j]
            left_values, right_values = values[i], values[j]
            pearson = pearson_matrix[i, j]
            spearman = np.corrcoef(rankdata(left_values), rankdata(right_values))[0, 1]
            cosine = float(left_values @ right_values / (np.linalg.norm(left_values) * np.linalg.norm(right_values)))
            left_sig = set(np.array(external_order)[significant.loc[left].to_numpy(dtype=bool)])
            right_sig = set(np.array(external_order)[significant.loc[right].to_numpy(dtype=bool)])
            union = left_sig | right_sig
            intersection = left_sig & right_sig
            same = sum(
                np.sign(matrix.loc[left, trait]) == np.sign(matrix.loc[right, trait])
                for trait in intersection
            )
            opposite = len(intersection) - same
            rows.append({
                "sleep_trait_1": left, "sleep_trait_2": right,
                "pearson_rg_profile": pearson, "spearman_rg_profile": spearman,
                "cosine_rg_profile": cosine,
                "fdr_significant_jaccard": len(intersection) / len(union) if union else 1.0,
                "signed_jaccard": (same - opposite) / len(union) if union else 1.0,
                "shared_significant_external_traits": len(intersection),
                "same_direction_shared_significant": same,
                "opposite_direction_shared_significant": opposite,
                "primary_external_traits_compared": len(external_order),
            })
    similarities = pd.DataFrame(rows)
    distance = np.clip(1 - pearson_matrix, 0, 2)
    np.fill_diagonal(distance, 0)
    sleep_linkage = linkage(squareform(distance, checks=False), method="average")
    sleep_link = pd.DataFrame(sleep_linkage, columns=["left_node", "right_node", "distance", "cluster_size"])
    sleep_link.insert(0, "step", np.arange(len(sleep_link)))
    order = pd.DataFrame({
        "dendrogram_order": np.arange(1, len(sleep_order) + 1),
        "sleep_trait": np.array(sleep_order)[leaves_list(sleep_linkage)],
    })
    centered = values - values.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :2] * singular[:2]
    variance = singular**2 / np.sum(singular**2)
    pca = pd.DataFrame({
        "sleep_trait": sleep_order, "PC1": scores[:, 0], "PC2": scores[:, 1],
        "PC1_variance_explained": variance[0], "PC2_variance_explained": variance[1],
        "interpretation_limit": "Descriptive rg-profile coordinates; not biological latent factors",
    })

    external_values = matrix.T.to_numpy(dtype=float)
    external_corr = np.corrcoef(external_values)
    external_distance = np.clip(1 - external_corr, 0, 2)
    np.fill_diagonal(external_distance, 0)
    external_linkage = linkage(squareform(external_distance, checks=False), method="average")
    best = None
    for clusters in range(2, min(9, len(external_order))):
        labels = fcluster(external_linkage, clusters, criterion="maxclust")
        silhouette = silhouette_precomputed(external_distance, labels)
        candidate = (silhouette, -clusters, labels)
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    assert best is not None
    labels = stable_labels(best[2], external_order)
    domain_by = panel.set_index("trait_id")["domain"].to_dict()
    cluster_domains = {
        cluster: sorted({domain_by[external_order[i]] for i in range(len(external_order)) if labels[i] == cluster})
        for cluster in sorted(set(labels))
    }
    external_clusters = pd.DataFrame({
        "external_trait": external_order,
        "external_domain": [domain_by[trait] for trait in external_order],
        "data_driven_cluster": labels,
        "selected_cluster_count": len(set(labels)),
        "mean_silhouette": best[0],
        "cluster_domains": [";".join(cluster_domains[label]) for label in labels],
        "cross_domain_cluster": [len(cluster_domains[label]) >= 2 for label in labels],
        "clustering_basis": "average-linkage correlation distance across 12 primary sleep-rg estimates; k selected by maximum mean silhouette over 2..8",
    })
    return similarities, sleep_link, order, pca, external_clusters


def silhouette_precomputed(distance: np.ndarray, labels: np.ndarray) -> float:
    values = []
    for i, label in enumerate(labels):
        same = np.where(labels == label)[0]
        same = same[same != i]
        if len(same) == 0:
            values.append(0.0)
            continue
        a = float(distance[i, same].mean())
        b = min(float(distance[i, np.where(labels == other)[0]].mean()) for other in set(labels) if other != label)
        values.append((b - a) / max(a, b) if max(a, b) else 0.0)
    return float(np.mean(values))


def stable_labels(labels: np.ndarray, names: list[str]) -> np.ndarray:
    groups = defaultdict(list)
    for label, name in zip(labels, names):
        groups[int(label)].append(name)
    mapping = {
        old: index + 1
        for index, old in enumerate(sorted(groups, key=lambda key: min(groups[key])))
    }
    return np.array([mapping[int(label)] for label in labels], dtype=int)


def effect_tiers(master: pd.DataFrame) -> pd.DataFrame:
    frame = master.loc[
        master["primary_or_sensitivity"].eq(PRIMARY) & master["fdr"].le(LOCKED_FDR)
    ].copy()
    bins = [-np.inf, .1, .2, .3, .4, .5, np.inf]
    labels = ["LT_0.10", "0.10_TO_LT_0.20", "0.20_TO_LT_0.30", "0.30_TO_LT_0.40", "0.40_TO_LT_0.50", "GE_0.50"]
    frame["effect_size_tier"] = pd.cut(frame["abs_rg"], bins=bins, labels=labels, right=False)
    for threshold in (.1, .2, .3, .4, .5):
        frame[f"abs_rg_ge_{str(threshold).replace('.', '_')}"] = frame["abs_rg"].ge(threshold)
    return frame[[
        "sleep_trait", "external_trait", "external_domain", "rg", "se", "z", "p", "fdr", "abs_rg",
        "effect_size_tier", "abs_rg_ge_0_1", "abs_rg_ge_0_2", "abs_rg_ge_0_3",
        "abs_rg_ge_0_4", "abs_rg_ge_0_5",
    ]].sort_values(["abs_rg", "fdr"], ascending=[False, True])


def result_confidence(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in master.itertuples(index=False):
        reasons = []
        if row.primary_or_sensitivity == SENSITIVITY:
            confidence = "SENSITIVITY_ONLY"
            reasons.append("external trait failed the locked h2 QC gate")
        elif row.sleep_h2 > 1 or row.external_h2 > 1 or min(row.sleep_h2_z, row.external_h2_z) < 6:
            confidence = "QC_CAUTION"
            if row.sleep_h2 > 1 or row.external_h2 > 1:
                reasons.append("nonphysical liability-scale h2 magnitude above one")
            if min(row.sleep_h2_z, row.external_h2_z) < 6:
                reasons.append("at least one h2 Z is below 6 despite passing the locked Z>=4 gate")
        else:
            moderate = []
            if min(row.sleep_h2_z, row.external_h2_z) < 10:
                moderate.append("at least one h2 Z is below 10")
            if max(row.sleep_h2_intercept, row.external_h2_intercept) > 1.10:
                moderate.append("at least one h2 intercept is above 1.10 but within the locked <=1.20 gate")
            if row.source_substitution:
                moderate.append("external phenotype uses a reviewed public source substitution")
            if row.external_source_id.startswith("finngen_"):
                moderate.append("external phenotype is a Finnish register endpoint")
            if row.SNP_overlap < 1_000_000:
                moderate.append("pair SNP overlap is below one million")
            confidence = "MODERATE_CONFIDENCE" if moderate else "HIGH_CONFIDENCE"
            reasons.extend(moderate or ["both traits have h2 Z>=10, intercepts<=1.10, >=1M overlapping SNPs, and no source substitution"])
        rows.append({
            "sleep_trait": row.sleep_trait, "external_trait": row.external_trait,
            "external_domain": row.external_domain, "rg": row.rg, "se": row.se,
            "p": row.p, "fdr": row.fdr, "primary_or_sensitivity": row.primary_or_sensitivity,
            "confidence_class": confidence, "confidence_reasons": "; ".join(reasons),
            "sleep_h2_z": row.sleep_h2_z, "external_h2_z": row.external_h2_z,
            "sleep_h2_intercept": row.sleep_h2_intercept,
            "external_h2_intercept": row.external_h2_intercept,
            "sleep_attenuation_ratio": row.sleep_h2_attenuation_ratio,
            "external_attenuation_ratio": row.external_h2_attenuation_ratio,
            "SNP_overlap": row.SNP_overlap, "gcov_intercept": row.gcov_intercept,
            "gcov_intercept_se": row.gcov_intercept_se,
            "gcov_intercept_z": row.gcov_intercept_z,
            "ancestry": row.ancestry, "measurement_category": row.sleep_measurement_category,
            "source_substitution": row.source_substitution,
            "source_caveat": row.source_substitution_caveat,
            "classification_rule": "Transparent rule hierarchy documented in scripts/102_analyze_phase1_atlas.py; no composite score",
        })
    return pd.DataFrame(rows)


def network_tables(master: pd.DataFrame, panel: pd.DataFrame, modes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)].copy()
    sleep = panel.loc[panel["domain"].eq("sleep"), "trait_id"].tolist()
    external = panel.loc[
        (~panel["domain"].eq("sleep")) & (~panel["trait_id"].isin(["t2d", "melanoma"])), "trait_id"
    ].tolist()
    names = [f"sleep::{name}" for name in sleep] + [f"external::{name}" for name in external]
    index = {name: i for i, name in enumerate(names)}
    significant = primary.loc[primary["fdr"].le(LOCKED_FDR)]
    adjacency = [[] for _ in names]
    weights = np.zeros((len(names), len(names)), dtype=float)
    edges = []
    for row in primary.itertuples(index=False):
        left, right = f"sleep::{row.sleep_trait}", f"external::{row.external_trait}"
        is_significant = row.fdr <= LOCKED_FDR
        edges.append({
            "sleep_trait": row.sleep_trait, "external_trait": row.external_trait,
            "external_domain": row.external_domain, "rg": row.rg, "abs_rg": row.abs_rg,
            "direction": row.direction, "p": row.p, "fdr": row.fdr,
            "locked_fdr_significant": is_significant, "abs_rg_ge_0_2": row.abs_rg >= .2,
            "positive_network": is_significant and row.rg > 0,
            "negative_network": is_significant and row.rg < 0,
        })
        if is_significant:
            a, b = index[left], index[right]
            adjacency[a].append(b); adjacency[b].append(a)
            weights[a, b] = weights[b, a] = abs(row.rg)
    betweenness = brandes(adjacency)
    communities = spectral_communities(weights, names)
    panel_by = panel.set_index("trait_id").to_dict("index")
    mode_by = modes.set_index("trait_id").to_dict("index")
    community_domains = defaultdict(set)
    for node, community in zip(names, communities):
        kind, trait = node.split("::", 1)
        if kind == "external":
            community_domains[int(community)].add(panel_by[trait]["domain"])
    node_rows = []
    signed_strength = defaultdict(float)
    for row in significant.itertuples(index=False):
        signed_strength[f"sleep::{row.sleep_trait}"] += row.rg
        signed_strength[f"external::{row.external_trait}"] += row.rg
    for i, node in enumerate(names):
        kind, trait = node.split("::", 1)
        meta = panel_by[trait]
        community = int(communities[i])
        node_rows.append({
            "node_id": node, "trait_id": trait, "node_type": kind.upper(),
            "label": meta["label"], "domain": meta["domain"],
            "measurement_category": mode_by[trait]["measurement_category"] if kind == "sleep" else "NA",
            "significant_degree": len(adjacency[i]),
            "weighted_degree_abs_rg": weights[i].sum(),
            "signed_strength_rg": signed_strength[node],
            "unweighted_betweenness": betweenness[i],
            "community_id": community,
            "community_external_domains": ";".join(sorted(community_domains[community])) or "NA",
            "cross_domain_community": len(community_domains[community]) >= 2,
            "centrality_limit": "Descriptive network topology; not biological causality",
        })
    return pd.DataFrame(node_rows), pd.DataFrame(edges)


def brandes(adjacency: list[list[int]]) -> np.ndarray:
    centrality = np.zeros(len(adjacency), dtype=float)
    for source in range(len(adjacency)):
        stack = []
        predecessors = [[] for _ in adjacency]
        sigma = np.zeros(len(adjacency)); sigma[source] = 1
        distance = np.full(len(adjacency), -1, dtype=int); distance[source] = 0
        queue = deque([source])
        while queue:
            vertex = queue.popleft(); stack.append(vertex)
            for neighbor in adjacency[vertex]:
                if distance[neighbor] < 0:
                    queue.append(neighbor); distance[neighbor] = distance[vertex] + 1
                if distance[neighbor] == distance[vertex] + 1:
                    sigma[neighbor] += sigma[vertex]; predecessors[neighbor].append(vertex)
        dependency = np.zeros(len(adjacency))
        while stack:
            vertex = stack.pop()
            for predecessor in predecessors[vertex]:
                dependency[predecessor] += (sigma[predecessor] / sigma[vertex]) * (1 + dependency[vertex])
            if vertex != source:
                centrality[vertex] += dependency[vertex]
    return centrality / 2


def spectral_communities(weights: np.ndarray, names: list[str]) -> np.ndarray:
    degree = weights.sum(axis=1)
    isolated = degree == 0
    safe = np.where(isolated, 1.0, degree)
    laplacian = np.eye(len(weights)) - np.diag(1 / np.sqrt(safe)) @ weights @ np.diag(1 / np.sqrt(safe))
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    maximum = min(8, len(weights) - 1)
    candidates = range(2, maximum + 1)
    k = max(candidates, key=lambda value: (eigenvalues[value] - eigenvalues[value - 1], -value))
    embedding = eigenvectors[:, :k]
    norm = np.linalg.norm(embedding, axis=1, keepdims=True)
    embedding = embedding / np.where(norm == 0, 1, norm)
    labels = deterministic_kmeans(embedding, k)
    return stable_labels(labels, names)


def deterministic_kmeans(values: np.ndarray, k: int) -> np.ndarray:
    seeds = [0]
    while len(seeds) < k:
        distances = np.min(
            np.stack([np.sum((values - values[index]) ** 2, axis=1) for index in seeds]), axis=0
        )
        distances[seeds] = -1
        seeds.append(int(np.argmax(distances)))
    centers = values[seeds].copy()
    labels = np.zeros(len(values), dtype=int)
    for _ in range(100):
        distances = np.stack([np.sum((values - center) ** 2, axis=1) for center in centers], axis=1)
        updated = np.argmin(distances, axis=1)
        if np.array_equal(updated, labels) and _ > 0:
            break
        labels = updated
        for cluster in range(k):
            members = values[labels == cluster]
            if len(members):
                centers[cluster] = members.mean(axis=0)
    return labels + 1


def bridges(sleep: pd.DataFrame, external_clusters: pd.DataFrame, network_nodes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in sleep.itertuples(index=False):
        if row.connected_domain_count >= 3:
            rows.append({
                "bridge_type": "SLEEP_DOMAIN_BRIDGE", "trait_id": row.sleep_trait,
                "domain": "sleep", "bridge_scope": row.connected_domains,
                "domain_count": row.connected_domain_count,
                "significant_connection_count": row.fdr_significant_connections,
                "bridge_strength": row.sum_abs_rg_significant,
                "bridge_tier": "BROAD_5PLUS_DOMAINS" if row.connected_domain_count >= 5 else "MULTIDOMAIN_3_TO_4",
                "objective_basis": "Locked-FDR connections span at least three external domains",
            })
    external_nodes = network_nodes.loc[network_nodes["node_type"].eq("EXTERNAL")].set_index("trait_id")
    for row in external_clusters.itertuples(index=False):
        if row.cross_domain_cluster:
            node = external_nodes.loc[row.external_trait]
            rows.append({
                "bridge_type": "EXTERNAL_PROFILE_CROSS_DOMAIN_CLUSTER_MEMBER",
                "trait_id": row.external_trait, "domain": row.external_domain,
                "bridge_scope": row.cluster_domains,
                "domain_count": len(row.cluster_domains.split(";")),
                "significant_connection_count": int(node.significant_degree),
                "bridge_strength": float(node.weighted_degree_abs_rg),
                "bridge_tier": f"PROFILE_CLUSTER_{row.data_driven_cluster}",
                "objective_basis": "Data-driven 12-rg profile cluster contains multiple predefined clinical domains",
            })
    return pd.DataFrame(rows)


def surprising_pairs(master: pd.DataFrame, measurement: pd.DataFrame) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)].copy()
    measurement_by = measurement.set_index("external_trait").to_dict("index")
    family_members = primary.groupby("sleep_phenotype_family")["sleep_trait"].unique().to_dict()
    rows = []
    for row in primary.itertuples(index=False):
        reasons = []
        metrics = []
        if row.abs_rg >= .4:
            reasons.append("VERY_LARGE_ABS_RG_GE_0.40"); metrics.append(f"abs_rg={row.abs_rg:.4g}")
        siblings = primary.loc[
            primary["external_trait"].eq(row.external_trait)
            & primary["sleep_trait"].isin(family_members[row.sleep_phenotype_family])
            & ~primary["sleep_trait"].eq(row.sleep_trait)
        ]
        if len(siblings):
            median = siblings["rg"].median()
            median_abs = siblings["abs_rg"].median()
            if row.abs_rg >= .1 and abs(median) >= .1 and np.sign(row.rg) != np.sign(median):
                reasons.append("OPPOSITE_SIGN_FROM_RELATED_SLEEP_FAMILY")
                metrics.append(f"related_median_rg={median:.4g}")
            if row.abs_rg >= median_abs + .2:
                reasons.append("ABS_RG_AT_LEAST_0.20_ABOVE_RELATED_FAMILY_MEDIAN")
                metrics.append(f"related_median_abs_rg={median_abs:.4g}")
            if row.fdr <= LOCKED_FDR and row.abs_rg >= .2 and siblings["abs_rg"].le(.05).any():
                reasons.append("STRONG_FOR_ONE_SLEEP_BUT_NEAR_ZERO_FOR_RELATED_PHENOTYPE")
            strong_siblings = siblings.loc[siblings["fdr"].le(LOCKED_FDR) & siblings["abs_rg"].ge(.2)]
            if row.fdr > LOCKED_FDR and row.abs_rg <= .05 and len(strong_siblings) >= 2:
                reasons.append("NOMINAL_NULL_SURROUNDED_BY_STRONG_RELATED_SLEEP_EFFECTS")
                metrics.append("strong_related=" + ";".join(strong_siblings["sleep_trait"]))
        measurement_row = measurement_by[row.external_trait]
        if measurement_row["significant_pair_direction_discordance"] and row.fdr <= LOCKED_FDR:
            reasons.append("SUBJECTIVE_OBJECTIVE_SIGN_DISCORDANCE_FOR_EXTERNAL_TRAIT")
        domain = primary.loc[
            primary["sleep_trait"].eq(row.sleep_trait) & primary["external_domain"].eq(row.external_domain)
        ]
        if row.fdr <= LOCKED_FDR and len(domain) >= 2 and domain["fdr"].le(LOCKED_FDR).sum() == 1:
            reasons.append("ONLY_SIGNIFICANT_CONNECTION_IN_SLEEP_BY_DOMAIN_BLOCK")
        if reasons:
            rows.append({
                "sleep_trait": row.sleep_trait, "external_trait": row.external_trait,
                "external_domain": row.external_domain, "rg": row.rg, "se": row.se,
                "p": row.p, "fdr": row.fdr, "abs_rg": row.abs_rg,
                "surprise_criteria": ";".join(sorted(set(reasons))),
                "criterion_metrics": ";".join(metrics) if metrics else "NA",
                "interpretation_limit": "Algorithmic outlier flag requiring biological/literature review; not itself evidence of novelty",
            })
    return pd.DataFrame(rows).sort_values(["abs_rg", "fdr"], ascending=[False, True])


def negative_rg(master: pd.DataFrame) -> pd.DataFrame:
    frame = master.loc[
        master["primary_or_sensitivity"].eq(PRIMARY)
        & master["fdr"].le(LOCKED_FDR) & master["rg"].lt(0)
    ].copy()
    frame["prior_observational_expectation"] = "NOT_YET_LITERATURE_CLASSIFIED"
    frame["direct_negative_rg_previously_reported"] = "NOT_YET_LITERATURE_CLASSIFIED"
    frame["antagonistic_local_followup_candidate"] = frame["abs_rg"].ge(.1)
    frame["claim_limit"] = "Negative global rg does not demonstrate local antagonistic pleiotropy"
    return frame[[
        "sleep_trait", "sleep_trait_label", "external_trait", "external_trait_label",
        "external_domain", "rg", "se", "z", "p", "fdr", "abs_rg",
        "prior_observational_expectation", "direct_negative_rg_previously_reported",
        "antagonistic_local_followup_candidate", "claim_limit",
    ]].sort_values("rg")


def multiple_testing(master: pd.DataFrame) -> pd.DataFrame:
    frame = master.copy()
    frame["locked_396_fdr_significant"] = frame["primary_or_sensitivity"].eq(PRIMARY) & frame["fdr"].le(.05)
    frame["primary_372_fdr_significant_sensitivity"] = (
        frame["primary_or_sensitivity"].eq(PRIMARY) & frame["fdr_primary_372_sensitivity"].le(.05)
    )
    frame["bonferroni_396_significant"] = frame["primary_or_sensitivity"].eq(PRIMARY) & frame["p"].le(.05 / 396)
    frame["locked_fdr_le_0_01"] = frame["primary_or_sensitivity"].eq(PRIMARY) & frame["fdr"].le(.01)
    frame["locked_fdr_le_0_001"] = frame["primary_or_sensitivity"].eq(PRIMARY) & frame["fdr"].le(.001)
    return frame[[
        "sleep_trait", "external_trait", "external_domain", "primary_or_sensitivity",
        "rg", "se", "p", "fdr", "fdr_primary_372_sensitivity",
        "locked_396_fdr_significant", "primary_372_fdr_significant_sensitivity",
        "bonferroni_396_significant", "locked_fdr_le_0_01", "locked_fdr_le_0_001",
    ]]


def leave_one_out(master: pd.DataFrame) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)]
    rows = []
    omissions = [("NONE", "ALL_12_SLEEP_TRAITS", primary)]
    omissions += [("SLEEP_TRAIT", trait, primary.loc[~primary["sleep_trait"].eq(trait)]) for trait in primary["sleep_trait"].unique()]
    omissions += [("EXTERNAL_DOMAIN", domain, primary.loc[~primary["external_domain"].eq(domain)]) for domain in sorted(primary["external_domain"].unique())]
    for kind, omitted, frame in omissions:
        significant = frame.loc[frame["fdr"].le(LOCKED_FDR)]
        sleep_counts = significant.groupby("sleep_trait").size()
        domain_counts = significant.groupby("external_domain").size()
        rows.append({
            "omission_type": kind, "omitted_unit": omitted, "remaining_pairs": len(frame),
            "significant_connections": len(significant),
            "significant_positive": int(significant["rg"].gt(0).sum()),
            "significant_negative": int(significant["rg"].lt(0).sum()),
            "mean_abs_rg_all_remaining": frame["abs_rg"].mean(),
            "mean_abs_rg_significant": significant["abs_rg"].mean() if len(significant) else np.nan,
            "top_remaining_sleep_hub": ";".join(sorted(sleep_counts.loc[sleep_counts.eq(sleep_counts.max())].index)) if len(sleep_counts) else "NA",
            "top_remaining_domain": ";".join(sorted(domain_counts.loc[domain_counts.eq(domain_counts.max())].index)) if len(domain_counts) else "NA",
            "locked_threshold_reused": True,
        })
    baseline = rows[0]["significant_connections"]
    for row in rows:
        row["change_from_all_12_significant_connections"] = row["significant_connections"] - baseline
    return pd.DataFrame(rows)


def source_substitution(master: pd.DataFrame, paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = read_tsv(paths["substitutions"])
    caveat = candidates.set_index("trait_id")["phenotype_caveat"].to_dict()
    replacement = set(candidates["trait_id"])
    pair_rows = master.loc[master["external_trait"].isin(replacement)].copy()
    pair_rows["replacement_group"] = np.where(
        pair_rows["external_source_id"].str.startswith("finngen_"), "FINNGEN_REPLACEMENT", "NON_FINNGEN_REPLACEMENT"
    )
    pair_rows["phenotype_caveat"] = pair_rows["external_trait"].map(caveat)
    pair_rows = pair_rows[[
        "sleep_trait", "external_trait", "external_domain", "replacement_group",
        "primary_or_sensitivity", "rg", "se", "p", "fdr", "abs_rg", "gcov_intercept",
        "gcov_intercept_se", "SNP_overlap", "external_h2", "external_h2_se",
        "external_h2_z", "external_h2_intercept", "phenotype_caveat",
    ]]
    summary_rows = []
    groups = {
        "ALL_REPLACEMENTS": master["external_trait"].isin(replacement),
        "FINNGEN_REPLACEMENTS": master["external_trait"].isin(
            set(candidates.loc[candidates["candidate_source_id"].str.startswith("finngen_"), "trait_id"])
        ),
        "NON_REPLACEMENT_EXTERNAL_TRAITS": ~master["external_trait"].isin(replacement),
    }
    for name, mask in groups.items():
        for tier in ("ALL_ROWS", PRIMARY, SENSITIVITY):
            selected = master.loc[mask & (True if tier == "ALL_ROWS" else master["primary_or_sensitivity"].eq(tier))]
            if selected.empty:
                continue
            traits = selected.drop_duplicates("external_trait")
            summary_rows.append({
                "comparison_group": name, "analysis_tier": tier,
                "external_trait_count": selected["external_trait"].nunique(),
                "pair_count": len(selected),
                "fdr_significant_pair_count": int(selected["fdr"].le(LOCKED_FDR).sum()),
                "mean_rg": selected["rg"].mean(), "median_rg": selected["rg"].median(),
                "mean_abs_rg": selected["abs_rg"].mean(), "median_abs_rg": selected["abs_rg"].median(),
                "mean_gcov_intercept": selected["gcov_intercept"].mean(),
                "median_gcov_intercept": selected["gcov_intercept"].median(),
                "mean_snp_overlap": selected["SNP_overlap"].mean(),
                "mean_external_h2_intercept": traits["external_h2_intercept"].mean(),
                "interpretation_limit": "Descriptive source-group comparison; endpoint and ancestry differences prevent causal attribution to FinnGen",
            })
    return pd.DataFrame(summary_rows), pair_rows


def cancellation_candidates(master: pd.DataFrame, hubs: pd.DataFrame) -> pd.DataFrame:
    primary = master.loc[master["primary_or_sensitivity"].eq(PRIMARY)].copy()
    hub_counts = hubs.set_index("sleep_trait")["fdr_significant_connections"].to_dict()
    hub_cutoff = float(hubs["fdr_significant_connections"].quantile(.75))
    external_counts = primary.loc[primary["fdr"].le(LOCKED_FDR)].groupby("external_trait").size().to_dict()
    family = primary.groupby("sleep_phenotype_family")["sleep_trait"].unique().to_dict()
    rows = []
    for row in primary.loc[primary["fdr"].gt(LOCKED_FDR) & primary["abs_rg"].le(.10)].itertuples(index=False):
        siblings = primary.loc[
            primary["external_trait"].eq(row.external_trait)
            & primary["sleep_trait"].isin(family[row.sleep_phenotype_family])
            & ~primary["sleep_trait"].eq(row.sleep_trait)
        ]
        strong_siblings = siblings.loc[siblings["fdr"].le(LOCKED_FDR) & siblings["abs_rg"].ge(.15)]
        reasons = []
        if len(strong_siblings) >= 1:
            reasons.append("RELATED_SLEEP_PHENOTYPE_HAS_SIGNIFICANT_ABS_RG_GE_0.15")
        if hub_counts[row.sleep_trait] >= hub_cutoff and external_counts.get(row.external_trait, 0) >= 3:
            reasons.append("NULL_EDGE_CONNECTS_TWO_BROAD_GLOBAL_HUBS")
        if row.abs_rg <= .05:
            reasons.append("GLOBAL_ABS_RG_LE_0.05")
        if (len(strong_siblings) or len(reasons) >= 3) and "GLOBAL_ABS_RG_LE_0.05" in reasons:
            rows.append({
                "sleep_trait": row.sleep_trait, "external_trait": row.external_trait,
                "external_domain": row.external_domain, "rg": row.rg, "se": row.se,
                "p": row.p, "fdr": row.fdr, "candidate_reasons": ";".join(reasons),
                "related_strong_pairs": ";".join(
                    f"{item.sleep_trait}:{item.rg:.4g}" for item in strong_siblings.itertuples()
                ) or "NA",
                "global_cancellation_status": "HYPOTHESIS_ONLY_NOT_DEMONSTRATED",
                "required_test": "LAVA or HDL-L with signed local estimates",
            })
    return pd.DataFrame(rows).sort_values(["external_domain", "external_trait", "sleep_trait"])


def write_outputs(root: Path, outputs: dict[str, object], validate_only: bool) -> None:
    for relative, value in outputs.items():
        path = root / relative
        if isinstance(value, pd.DataFrame):
            payload = frame_text(value)
        elif isinstance(value, str):
            payload = value if value.endswith("\n") else value + "\n"
        else:
            payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
        if validate_only:
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"deep-analysis artifact drifted: {relative}")
        else:
            atomic_text(path, payload)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--permutations", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if args.permutations < 1000:
        fail("at least 1,000 whole-profile domain permutations are required")
    root = Path(args.root).resolve()
    diagnostics, paths, log_paths = validate_inputs(root)
    master, panel, _h2, modes = build_master(root, paths)
    summary_json, summary_md = global_summary(master)
    sleep = sleep_hubs(master, panel)
    external = external_hubs(master, modes, panel)
    domains = domain_architecture(master, args.permutations, args.seed)
    measurement = measurement_comparison(master, modes)
    similarity, sleep_linkage, sleep_order, pca, external_clusters = similarity_and_clustering(master, panel)
    effects = effect_tiers(master)
    confidence = result_confidence(master)
    network_nodes, network_edges = network_tables(master, panel, modes)
    bridge = bridges(sleep, external_clusters, network_nodes)
    surprise = surprising_pairs(master, measurement)
    negatives = negative_rg(master)
    testing = multiple_testing(master)
    leaveout = leave_one_out(master)
    substitutions, substitution_pairs = source_substitution(master, paths)
    cancellation = cancellation_candidates(master, sleep)

    input_files = [path for name, path in paths.items() if path.is_file() and name != "requested_phase1_completion"] + log_paths
    hashes = pd.DataFrame([
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(set(input_files))
    ])
    diagnostics["analyzed_input_file_count"] = len(hashes)
    diagnostics["input_hash_manifest_sha256"] = hashlib.sha256(frame_text(hashes).encode()).hexdigest()
    diagnostics["analysis_parameters"] = {
        "locked_fdr_threshold": LOCKED_FDR,
        "domain_permutations": args.permutations,
        "random_seed": args.seed,
        "domain_permutation_unit": "whole external-trait 12-rg profile",
        "sleep_profile_scope": "31 primary external traits; T2D and melanoma excluded",
    }
    outputs: dict[str, object] = {
        "results/analysis/input_verification.json": diagnostics,
        "results/analysis/input_hashes.tsv": hashes,
        "results/analysis/phase1_master_analysis.tsv": master,
        "results/analysis/global_summary.json": summary_json,
        "results/analysis/global_summary.md": summary_md,
        "results/analysis/sleep_trait_hubs.tsv": sleep,
        "results/analysis/external_trait_hubs.tsv": external,
        "results/analysis/domain_architecture.tsv": domains,
        "results/analysis/measurement_mode_comparison.tsv": measurement,
        "results/analysis/sleep_profile_similarity.tsv": similarity,
        "results/analysis/sleep_profile_linkage.tsv": sleep_linkage,
        "results/analysis/sleep_profile_dendrogram_order.tsv": sleep_order,
        "results/analysis/sleep_profile_pca.tsv": pca,
        "results/analysis/external_trait_clusters.tsv": external_clusters,
        "results/analysis/cross_domain_bridges.tsv": bridge,
        "results/analysis/effect_size_tiers.tsv": effects,
        "results/analysis/negative_rg.tsv": negatives,
        "results/analysis/surprising_pairs.tsv": surprise,
        "results/analysis/result_confidence.tsv": confidence,
        "results/analysis/source_substitution_sensitivity.tsv": substitutions,
        "results/analysis/source_substitution_pairs.tsv": substitution_pairs,
        "results/analysis/network_nodes.tsv": network_nodes,
        "results/analysis/network_edges.tsv": network_edges,
        "results/analysis/multiple_testing_sensitivity.tsv": testing,
        "results/analysis/leave_one_out_sensitivity.tsv": leaveout,
        "results/analysis/global_cancellation_candidates.tsv": cancellation,
    }
    write_outputs(root, outputs, args.validate_only)
    print(
        f"PHASE1_DEEP_ANALYSIS_{'VALID' if args.validate_only else 'BUILT'} "
        f"pairs={len(master)} primary_discoveries={int(master['locked_primary_significant'].sum())} "
        f"negative_discoveries={len(negatives)} inputs={len(hashes)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
