#!/usr/bin/env python3
"""Run and collate the frozen 12 sleep x 7 latent-frailty secondary rg family."""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
SLEEP_TRAITS = [
    "insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype",
    "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency",
    "accel_sleep_duration", "sleep_timing",
]
FACTOR_TRAITS = [
    "frailty_general", "frailty_factor_1", "frailty_factor_2",
    "frailty_factor_3", "frailty_factor_4", "frailty_factor_5",
    "frailty_factor_6",
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def run(command: list[str], *, cwd: Path = REPO) -> None:
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_rg_log(path: Path, factors: list[str]) -> bool:
    if not path.is_file():
        return False
    lines = path.read_text(errors="replace").splitlines()
    try:
        table_start = next(i for i, line in enumerate(lines)
                           if "Summary of Genetic Correlation Results" in line)
    except StopIteration:
        return False
    body = []
    for line in lines[table_start + 2:]:
        if not line.strip():
            break
        if line.startswith("Analysis finished"):
            break
        body.append(line)
    return len(body) == len(factors) and all(
        any(f"{trait}.sumstats.gz" in line for line in body) for trait in factors
    )


def run_family(external_root: Path, repo: Path = REPO) -> Path:
    external_root = external_root.resolve()
    workspace = external_root / "analysis-workspace/frailty_v1"
    sleep_munged = workspace / "munged"
    factor_munged = workspace / "munged_secondary"
    ref = external_root / "analysis-workspace/reference/eur_w_ld_chr"
    logs = repo / "frailty_paper/results/frailty_v1/logs_rg_latent_factors"
    logs.mkdir(parents=True, exist_ok=True)

    h2_sleep = read_tsv(repo / "frailty_paper/results/frailty_v1/h2_sleep_panel.tsv").set_index("trait")
    h2_factors = read_tsv(repo / "frailty_paper/results/frailty_v1/h2_latent_factors.tsv").set_index("trait")
    if set(SLEEP_TRAITS).difference(h2_sleep.index) or not h2_sleep.loc[SLEEP_TRAITS, "verdict"].eq("PASS").all():
        raise RuntimeError("the locked sleep h2 gate is incomplete or failed")
    if set(FACTOR_TRAITS).difference(h2_factors.index) or not h2_factors.loc[FACTOR_TRAITS, "verdict"].eq("PASS").all():
        raise RuntimeError("the locked latent-factor h2 gate is incomplete or failed")

    munged_paths: dict[str, Path] = {}
    for trait, directory in [(t, sleep_munged) for t in SLEEP_TRAITS] + [(t, factor_munged) for t in FACTOR_TRAITS]:
        path = directory / f"{trait}.sumstats.gz"
        if not path.is_file():
            raise RuntimeError(f"munged file missing after h2 gate: {path}")
        # The h2 jobs already consumed every file; avoid another full decompression here.
        munged_paths[trait] = path

    python = repo / ".ldsc-env/bin/python"
    ldsc = repo / "ldsc/ldsc.py"
    if not python.is_file() or not ldsc.is_file():
        raise RuntimeError("the pinned LDSC runtime or source script is missing")
    failures: list[str] = []
    for sleep in SLEEP_TRAITS:
        log_prefix = logs / f"rg_{sleep}"
        log_file = Path(str(log_prefix) + ".log")
        if valid_rg_log(log_file, FACTOR_TRAITS):
            print(f"RESUME: valid 7-factor rg log for {sleep}", flush=True)
            continue
        factors = ",".join(str(munged_paths[t]) for t in FACTOR_TRAITS)
        command = [str(python), str(ldsc), "--rg", f"{munged_paths[sleep]},{factors}",
                   "--ref-ld-chr", str(ref) + "/", "--w-ld-chr", str(ref) + "/",
                   "--out", str(log_prefix)]
        try:
            run(command, cwd=repo)
            if not valid_rg_log(log_file, FACTOR_TRAITS):
                raise RuntimeError("LDSC log lacks one or more expected factor result rows")
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            failures.append(f"{sleep}: {exc}")
            print(f"FAILED: {sleep}: {exc}", file=sys.stderr, flush=True)

    if failures:
        raise RuntimeError("not publishing a reduced rg family; repair/retry these runs: " + "; ".join(failures))

    combined_config = repo / "frailty_paper/config/secondary_sleep_factor_panel_v1.tsv"
    inclusion = repo / "frailty_paper/config/secondary_sleep_factor_inclusion_v1.tsv"
    raw = repo / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors_collated.tsv"
    run([sys.executable, "scripts/05_collate.py", "--mode", "rg", "--config", str(combined_config),
         "--inclusion", str(inclusion), "--logdir", str(logs), "--out", str(raw)], cwd=repo)
    results = pd.read_csv(raw, sep="\t")
    if len(results) != 84 or results.duplicated(["sleep_trait", "disease_trait"]).any():
        raise RuntimeError(f"expected exactly 84 unique secondary rg pairs, got {len(results)}")

    # Cohort provenance is deliberately conservative: the latent study used
    # UKB-containing source GWAS, but exact pairwise intersections are unknown.
    results["analysis_family"] = "secondary_sleep_x_latent_frailty"
    results["family_denominator"] = 84
    results["cohort_overlap_status"] = "EXPECTED_OR_POSSIBLE; exact overlap UNKNOWN"
    results["interpretation_status"] = "SENSITIVITY_ONLY"
    results["claim_limit"] = "global genetic correlation; not independent replication or causality"
    results = results.rename(columns={"p": "p_value", "fdr": "q_value"})
    out = repo / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv"
    results.to_csv(out, sep="\t", index=False, float_format="%.12g")
    script_hash = sha256(Path(__file__))
    config_hash = sha256(repo / "frailty_paper/config/secondary_frailty_factors_v1.yaml")
    execution_hash_file = repo / "frailty_paper/config/secondary_sleep_factor_rg_execution_v1.sha256"
    execution_hash = sha256(execution_hash_file)
    panel_hash = sha256(repo / "frailty_paper/config/secondary_sleep_factor_panel_v1.tsv")
    inclusion_hash = sha256(repo / "frailty_paper/config/secondary_sleep_factor_inclusion_v1.tsv")
    ldsc_hash = sha256(ldsc)
    input_hashes: dict[Path, str] = {}
    manifest_rows = []
    for sleep in SLEEP_TRAITS:
        log_file = logs / f"rg_{sleep}.log"
        for factor in FACTOR_TRAITS:
            left, right = munged_paths[sleep], munged_paths[factor]
            for path in (left, right):
                if path not in input_hashes:
                    input_hashes[path] = sha256(path)
            manifest_rows.append({
                "sleep_trait": sleep,
                "frailty_factor": factor,
                "sleep_munged_file": str(left),
                "sleep_munged_sha256": input_hashes[left],
                "factor_munged_file": str(right),
                "factor_munged_sha256": input_hashes[right],
                "ldsc_log": str(log_file),
                "ldsc_log_sha256": sha256(log_file),
                "secondary_family_config_sha256": config_hash,
                "execution_hash_file_sha256": execution_hash,
                "combined_trait_panel_sha256": panel_hash,
                "family_inclusion_manifest_sha256": inclusion_hash,
                "ldsc_source_sha256": ldsc_hash,
                "runner_sha256": script_hash,
                "analysis_status": "SENSITIVITY_ONLY_COMPLETE",
            })
    manifest = repo / "frailty_paper/results/frailty_v1/rg_sleep_latent_factors_run_manifest.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest, sep="\t", index=False)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path,
                        default=Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1"))
    args = parser.parse_args()
    output = run_family(args.external_root)
    print(f"PASS: complete 84-pair sensitivity-only family; wrote {output}", flush=True)


if __name__ == "__main__":
    main()
