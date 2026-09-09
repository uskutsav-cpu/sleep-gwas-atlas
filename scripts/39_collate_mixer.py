#!/usr/bin/env python3
"""Collate complete 20-replicate MiXeR summaries into canonical atlas tables."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import mixer_tasks


UNIVARIATE_FIELDS = [
    "trait_id", "mixer_version", "container_digest", "reference_commit", "input_scope",
    "replicates", "pi_mean", "pi_std", "sig2_beta_mean", "sig2_beta_std",
    "sig2_zero_mean", "sig2_zero_std", "h2_mean", "h2_std", "n_causal_p9_mean",
    "n_causal_p9_std", "AIC", "BIC", "power_interpretation", "bivariate_eligibility",
    "source_summary",
]
BIVARIATE_FIELDS = [
    "sleep_trait", "non_sleep_trait", "mixer_version", "container_digest",
    "reference_commit", "replicates", "dice_mean", "dice_std", "pi1_mean", "pi1_std",
    "pi2_mean", "pi2_std", "pi12_mean", "pi12_std", "n_shared_p9_mean",
    "n_shared_p9_std", "rho_beta_mean", "rho_beta_std", "rg_mean", "rg_std",
    "fraction_concordant_mean", "fraction_concordant_std", "best_vs_min_AIC",
    "best_vs_min_BIC", "best_vs_max_AIC", "best_vs_max_BIC", "source_summary",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_summary(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or not path.stat().st_size:
        raise SystemExit(f"ERROR: missing combined MiXeR summary: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"ERROR: empty combined MiXeR summary: {path}")
    return rows


def number(row: dict[str, str], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, ValueError) as exc:
        raise SystemExit(f"ERROR: missing/non-numeric MiXeR field {key}") from exc
    if not math.isfinite(value):
        raise SystemExit(f"ERROR: non-finite MiXeR field {key}")
    return value


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    if any(list(row) != fields for row in rows):
        raise SystemExit(f"ERROR: MiXeR rows differ from the canonical schema: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def publish_phase(
    root: Path, policy: dict[str, object], panel_path: Path, phase: str,
    path: Path, fields: list[str], rows: list[dict[str, object]],
) -> None:
    provenance_path = root / str(policy[f"{phase}_result_provenance_path"])
    if path.exists() != provenance_path.exists():
        raise SystemExit(f"ERROR: partial immutable MiXeR {phase} publication exists")
    payload = mixer_tasks.table_text(fields, rows)
    if path.exists():
        if path.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"ERROR: existing MiXeR {phase} table differs from recomputation")
        mixer_tasks.validate_phase_results(
            root, policy, panel_path, phase, path, len(rows),
        )
        return
    write_tsv(path, fields, rows)
    mixer_tasks.seal_phase_results(
        root, policy, panel_path, phase, path, len(rows),
    )


def univariate(root: Path, policy: dict[str, object], panel: list[dict[str, str]]) -> list[dict[str, object]]:
    rows = []
    for trait in (row["trait_id"] for row in panel):
        replicate_files = [
            root / f"results/mixer/univariate/{trait}.{stage}.rep{replicate}.json"
            for replicate in range(1, policy["fit_replicates"] + 1)
            for stage in ("fit", "test")
        ]
        if not all(path.is_file() and path.stat().st_size for path in replicate_files):
            raise SystemExit(f"ERROR: incomplete 20-replicate fit1/test1 family for {trait}")
        path = root / f"results/mixer/univariate/{trait}.fit.summary.csv"
        summary = read_summary(path)
        if len(summary) != 1:
            raise SystemExit(f"ERROR: expected one combined univariate row for {trait}")
        source = summary[0]
        aic = number(source, "AIC")
        bic = number(source, "BIC")
        eligibility = "ELIGIBLE" if aic > policy["univariate_aic_threshold"] else "INELIGIBLE_LOW_POWER"
        power = "SUPPORTED_AIC_BIC" if aic > 0 and bic > 0 else (
            "BORDERLINE_POWER_AIC_ONLY" if aic > 0 else "INSUFFICIENT_POWER"
        )
        rows.append({
            "trait_id": trait,
            "mixer_version": policy["mixer_release"],
            "container_digest": policy["container_amd64_manifest_digest"],
            "reference_commit": policy["mixer_reference_commit_at_lock"],
            "input_scope": "FULL_POST_QC_AUTOSOMAL",
            "replicates": policy["fit_replicates"],
            "pi_mean": number(source, "pi (mean)"),
            "pi_std": number(source, "pi (std)"),
            "sig2_beta_mean": number(source, "sig2_beta (mean)"),
            "sig2_beta_std": number(source, "sig2_beta (std)"),
            "sig2_zero_mean": number(source, "sig2_zero (mean)"),
            "sig2_zero_std": number(source, "sig2_zero (std)"),
            "h2_mean": number(source, "h2 (mean)"),
            "h2_std": number(source, "h2 (std)"),
            "n_causal_p9_mean": number(source, "nc@p9 (mean)"),
            "n_causal_p9_std": number(source, "nc@p9 (std)"),
            "AIC": aic,
            "BIC": bic,
            "power_interpretation": power,
            "bivariate_eligibility": eligibility,
            "source_summary": str(path.relative_to(root)),
        })
    return rows


def bivariate(
    root: Path,
    policy: dict[str, object],
    panel: list[dict[str, str]],
    univ: list[dict[str, object]],
) -> list[dict[str, object]]:
    eligible = {row["trait_id"] for row in univ if row["bivariate_eligibility"] == "ELIGIBLE"}
    sleeps = [row["trait_id"] for row in panel if row["domain"] == "sleep" and row["trait_id"] in eligible]
    diseases = [row["trait_id"] for row in panel if row["domain"] != "sleep" and row["trait_id"] in eligible]
    rows = []
    for sleep in sleeps:
        for disease in diseases:
            replicate_files = [
                root / f"results/mixer/bivariate/{sleep}_vs_{disease}.{stage}.rep{replicate}.json"
                for replicate in range(1, policy["fit_replicates"] + 1)
                for stage in ("fit", "test")
            ]
            if not all(path.is_file() and path.stat().st_size for path in replicate_files):
                raise SystemExit(f"ERROR: incomplete 20-replicate fit2/test2 family for {sleep}/{disease}")
            path = root / f"results/mixer/bivariate/{sleep}_vs_{disease}.csv"
            summary = read_summary(path)
            fit_rows = [row for row in summary if ".fit.json" in row.get("fname", "")]
            if len(fit_rows) != 1:
                raise SystemExit(f"ERROR: expected one combined fit row for {sleep}/{disease}")
            source = fit_rows[0]
            rows.append({
                "sleep_trait": sleep,
                "non_sleep_trait": disease,
                "mixer_version": policy["mixer_release"],
                "container_digest": policy["container_amd64_manifest_digest"],
                "reference_commit": policy["mixer_reference_commit_at_lock"],
                "replicates": policy["fit_replicates"],
                "dice_mean": number(source, "dice (mean)"),
                "dice_std": number(source, "dice (std)"),
                "pi1_mean": number(source, "pi1 (mean)"),
                "pi1_std": number(source, "pi1 (std)"),
                "pi2_mean": number(source, "pi2 (mean)"),
                "pi2_std": number(source, "pi2 (std)"),
                "pi12_mean": number(source, "pi12 (mean)"),
                "pi12_std": number(source, "pi12 (std)"),
                "n_shared_p9_mean": number(source, "nc12@p9 (mean)"),
                "n_shared_p9_std": number(source, "nc12@p9 (std)"),
                "rho_beta_mean": number(source, "rho_beta (mean)"),
                "rho_beta_std": number(source, "rho_beta (std)"),
                "rg_mean": number(source, "rg (mean)"),
                "rg_std": number(source, "rg (std)"),
                "fraction_concordant_mean": number(source, "fraction_concordant_within_shared (mean)"),
                "fraction_concordant_std": number(source, "fraction_concordant_within_shared (std)"),
                "best_vs_min_AIC": number(source, "best_vs_min_AIC"),
                "best_vs_min_BIC": number(source, "best_vs_min_BIC"),
                "best_vs_max_AIC": number(source, "best_vs_max_AIC"),
                "best_vs_max_BIC": number(source, "best_vs_max_BIC"),
                "source_summary": str(path.relative_to(root)),
            })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--univariate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads((root / "config/mixer_analysis_policy.json").read_text(encoding="utf-8"))
    panel_path = root / "config/analysis_panel.tsv"
    panel = read_tsv(panel_path)
    if len(panel) != 45:
        raise SystemExit("ERROR: expected locked 45-trait panel")
    univ = univariate(root, policy, panel)
    publish_phase(
        root, policy, panel_path, "univariate",
        root / "results/tables/mixer_univariate.tsv", UNIVARIATE_FIELDS, univ,
    )
    if args.univariate_only:
        eligible = sum(row["bivariate_eligibility"] == "ELIGIBLE" for row in univ)
        print(f"Published 45 univariate MiXeR summaries; {eligible} eligible")
        return 0
    pair_rows = bivariate(root, policy, panel, univ)
    publish_phase(
        root, policy, panel_path, "bivariate",
        root / "results/tables/mixer_bivariate.tsv", BIVARIATE_FIELDS, pair_rows,
    )
    print(f"Published MiXeR: 45 univariate traits and {len(pair_rows)} eligible sleep-by-non-sleep pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
