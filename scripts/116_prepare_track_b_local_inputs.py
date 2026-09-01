#!/usr/bin/env python3
"""Prepare and freeze the three-pair Track B local-analysis input family."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import shutil
from pathlib import Path

import numpy as np


POLICY = Path("config/track_b_local_analysis_policy.json")
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
PANEL = Path("config/analysis_panel.tsv")
INTERCEPTS = Path("results/tables/ldsc_intercept_45x45.tsv")
LAVA_LOCUS = Path("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
INPUT_INFO = Path("results/track_b/lava_input_info.tsv")
PAIR_INPUT = Path("results/track_b/lava_pair_manifest.tsv")
OVERLAP = Path("results/track_b/lava_sample_overlap.txt")
PROVENANCE = Path("results/track_b/lava_input_provenance.tsv")
CONDITIONAL = Path("results/track_b/local_conditional_manifest.tsv")
ROBUSTNESS = Path("results/track_b/local_method_robustness_plan.tsv")
HANDOFF = Path("results/track_b/LOCAL_ANALYSIS_HANDOFF.md")
LOCK = Path("results/track_b/local_analysis_input.lock.json")
REFERENCE_PROVENANCE = Path("ref/lava/ukb_v1.1/reference.provenance.json")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def load_policy() -> dict[str, object]:
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    if (
        policy["lava_version"] != "0.1.5"
        or policy["expected_loci"] != 2495
        or policy["expected_discovery_pairs"] != 3
        or policy["expected_analysis_traits"] != 8
        or policy["planned_univariate_tests"] != 2495 * 8
        or policy["planned_bivariate_pair_locus_family_max"] != 2495 * 3
    ):
        raise SystemExit("ERROR: Track B local policy scope drifted")
    expected = policy["univariate_alpha"] / policy["planned_univariate_tests"]
    if not math.isclose(expected, policy["univariate_p_threshold"], rel_tol=1e-12):
        raise SystemExit("ERROR: Track B local univariate threshold is inconsistent")
    return policy


def validate_locus_file(policy: dict[str, object]) -> None:
    if not LAVA_LOCUS.is_file():
        raise SystemExit("ERROR: LAVA 2,495-locus definition is missing")
    with LAVA_LOCUS.open(encoding="utf-8") as handle:
        lines = sum(1 for _ in handle)
    if lines != int(policy["expected_loci"]) + 1:
        raise SystemExit(f"ERROR: expected 2,495 LAVA loci, found {lines-1}")


def validate_sumstats(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing LAVA input: {path}")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    required = {"SNP", "A1", "A2", "Z", "N"}
    if not required.issubset(fields):
        raise SystemExit(f"ERROR: {path} lacks LAVA fields: {sorted(required-set(fields))}")
    return ",".join(fields)


def input_rows(policy: dict[str, object], panel: dict[str, dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    info: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    for trait in policy["trait_order"]:
        row = panel[str(trait)]
        path = Path("data/munged") / f"{trait}.sumstats.gz"
        schema = validate_sumstats(path)
        if row["type"] == "binary":
            cases, controls, prevalence = row["ncase"], row["ncontrol"], row["pop_prev"]
            if not cases.isdigit() or not controls.isdigit() or not 0 < float(prevalence) < 1:
                raise SystemExit(f"ERROR: unresolved binary metadata for {trait}")
        else:
            cases = controls = prevalence = "NA"
        info.append({
            "phenotype": trait, "cases": cases, "controls": controls,
            "prevalence": prevalence, "filename": str(path),
        })
        provenance.append({
            "phenotype": trait, "filename": str(path), "bytes": path.stat().st_size,
            "sha256": sha256(path), "observed_columns": schema,
            "validation_status": "VALIDATED_FOR_LAVA_LOCAL_NOT_FINE_MAPPING",
        })
    return info, provenance


def pair_rows(policy: dict[str, object]) -> list[dict[str, object]]:
    manifest = read_tsv(PAIR_MANIFEST)
    by_id = {row["pair_id"]: row for row in manifest}
    if set(by_id) != {"A", "B", "CONTROL"}:
        raise SystemExit("ERROR: frozen Track B pair family drifted")
    rows = []
    for pair_id in policy["pair_order"]:
        row = by_id[str(pair_id)]
        rows.append({
            "pair_id": pair_id,
            "trait1": row["sleep_trait"],
            "trait2": row["external_trait"],
            "discovery_rg": row["discovery_rg"],
            "discovery_SE": row["SE"],
            "discovery_P": row["P"],
            "discovery_FDR": row["FDR"],
            "planned_loci": policy["expected_loci"],
            "bivariate_family": "BH_FDR_ACROSS_ALL_ACTUALLY_TESTED_PAIR_LOCUS_ROWS_FOR_THREE_FROZEN_PAIRS",
            "pair_replacement": "FORBIDDEN",
        })
    return rows


def overlap_text(policy: dict[str, object]) -> tuple[str, dict[str, float]]:
    all_rows = read_tsv(INTERCEPTS)
    all_traits = [row["trait_id"] for row in all_rows]
    if len(all_rows) != 45 or list(all_rows[0].keys())[1:] != all_traits:
        raise SystemExit("ERROR: canonical 45-trait intercept matrix drifted")
    index = {trait: row for trait, row in zip(all_traits, all_rows)}
    traits = [str(value) for value in policy["trait_order"]]
    matrix = np.array([[float(index[first][second]) for second in traits] for first in traits])
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, atol=1e-10):
        raise SystemExit("ERROR: selected LDSC intercept matrix is invalid")
    diagonal = np.diag(matrix)
    correlation = matrix / np.sqrt(np.outer(diagonal, diagonal))
    correlation = np.round((correlation + correlation.T) / 2, 5)
    np.fill_diagonal(correlation, 1.0)
    eigenvalues = np.linalg.eigvalsh(correlation)
    if eigenvalues[0] <= 0:
        raise SystemExit("ERROR: Track B overlap matrix is not positive definite")
    buffer = io.StringIO()
    buffer.write(" " + " ".join(traits) + "\n")
    for trait, values in zip(traits, correlation):
        buffer.write(trait + " " + " ".join(f"{value:.5f}" for value in values) + "\n")
    return buffer.getvalue(), {
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "condition_number": float(eigenvalues[-1] / eigenvalues[0]),
    }


def conditional_rows(policy: dict[str, object]) -> list[dict[str, object]]:
    pairs = {row["pair_id"]: row for row in read_tsv(PAIR_MANIFEST)}
    output = []
    for pair_id in policy["pair_order"]:
        covariates = policy["conditional_covariates"][pair_id]
        if not covariates:
            output.append({
                "pair_id": pair_id, "trait1": pairs[pair_id]["sleep_trait"],
                "trait2": pairs[pair_id]["external_trait"], "covariate": "NONE",
                "rationale": policy["conditional_covariate_rationale"][pair_id],
                "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS", "execution_status": "NOT_PLANNED_FOR_CONTROL",
            })
            continue
        for covariate in covariates:
            output.append({
                "pair_id": pair_id, "trait1": pairs[pair_id]["sleep_trait"],
                "trait2": pairs[pair_id]["external_trait"], "covariate": covariate,
                "rationale": policy["conditional_covariate_rationale"][pair_id],
                "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS", "execution_status": "PENDING_LAVA_REFERENCE",
            })
    return output


def robustness_rows() -> list[dict[str, object]]:
    return [
        {
            "method": "HDL-L", "status": "BLOCKED_BY_SOFTWARE_AND_REFERENCE",
            "reference": "NO_PINNED_ANCESTRY_MATCHED_HDL_L_REFERENCE_PRESENT",
            "comparison_rule": "Compare locus overlap, sign, approximate magnitude and corrected significance; preserve discordance.",
            "allowed_labels": "CONCORDANT;PARTIAL_CONCORDANCE;DISCORDANT",
        },
        {
            "method": "rho-HESS", "status": "FALLBACK_REQUIRES_SEPARATE_PINNED_REFERENCE_AND_CONTRACT",
            "reference": "NO_TRACK_B_SPECIFIC_SEALED_RHO_HESS_REFERENCE_PRESENT",
            "comparison_rule": "May be used only as an independently implemented local-correlation sensitivity, not as a proxy for LAVA.",
            "allowed_labels": "CONCORDANT;PARTIAL_CONCORDANCE;DISCORDANT",
        },
    ]


def handoff_text(policy: dict[str, object], free_bytes: int) -> str:
    return f"""# Track B local-analysis production handoff

Current status: `BLOCKED_BY_DATA` and `BLOCKED_BY_COMPUTE`.

The three frozen pairs, eight required analysis traits, 2,495 loci, LDSC overlap submatrix, conditional covariates, correction families, and input hashes are prepared and locked. The official LAVA UK Biobank European v1.1 LD payload is absent.

Current free space at preparation: {free_bytes} bytes ({free_bytes/1024**3:.3f} GiB).
Policy minimum for checksum-ledgered acquisition and extraction: {policy['reference_minimum_free_bytes']} bytes (35 GiB).

## Minimum production environment

- x86_64 Linux recommended
- R 4.3.x with LAVA 0.1.5
- at least 16 GiB RAM; 32 GiB recommended for parallel or conditional work
- at least 35 GiB free for reference acquisition; 60 GiB recommended for outputs/checkpoints
- the exact repository commit and ignored dense/munged inputs whose hashes are in `local_analysis_input.lock.json`

## Exact execution order

```bash
python3 scripts/110_track_b_checkpoint.py --verify
python3 scripts/111_freeze_track_b_pairs.py --verify
python3 scripts/113_freeze_track_b_replication_sources.py --verify
python3 scripts/115_build_track_b_dense_qc.py --verify
python3 scripts/116_prepare_track_b_local_inputs.py --verify
bash scripts/32_download_lava_reference.sh --download
python3 scripts/lava_contract.py --verify-reference --rehash
```

Do not run the existing 45-trait/396-pair LAVA result family as a substitute for Track B. A Track B runner must consume `results/track_b/lava_pair_manifest.tsv`, use all 2,495 loci, run univariate h2 for all eight predeclared traits, and apply BH FDR across every actually tested locus row for the three frozen pairs. Conditional tests must use only the covariates in `local_conditional_manifest.tsv`.

No local result exists yet. No empty table or synthetic output is presented as science.
"""


def build() -> dict[Path, str]:
    policy = load_policy()
    validate_locus_file(policy)
    panel = {row["trait_id"]: row for row in read_tsv(PANEL)}
    if [trait for trait in policy["trait_order"] if trait not in panel]:
        raise SystemExit("ERROR: Track B local trait absent from locked panel")
    info, provenance = input_rows(policy, panel)
    pairs = pair_rows(policy)
    overlap, overlap_diagnostics = overlap_text(policy)
    conditional = conditional_rows(policy)
    robustness = robustness_rows()
    free_bytes = shutil.disk_usage(Path(".")).free
    outputs = {
        INPUT_INFO: tsv_text(["phenotype", "cases", "controls", "prevalence", "filename"], info),
        PAIR_INPUT: tsv_text([
            "pair_id", "trait1", "trait2", "discovery_rg", "discovery_SE", "discovery_P",
            "discovery_FDR", "planned_loci", "bivariate_family", "pair_replacement",
        ], pairs),
        OVERLAP: overlap,
        PROVENANCE: tsv_text([
            "phenotype", "filename", "bytes", "sha256", "observed_columns", "validation_status",
        ], provenance),
        CONDITIONAL: tsv_text([
            "pair_id", "trait1", "trait2", "covariate", "rationale", "selection_timing", "execution_status",
        ], conditional),
        ROBUSTNESS: tsv_text(["method", "status", "reference", "comparison_rule", "allowed_labels"], robustness),
        HANDOFF: handoff_text(policy, free_bytes),
    }
    immutable_hashes = {str(path): hashlib.sha256(value.encode()).hexdigest() for path, value in outputs.items() if path != HANDOFF}
    lock = {
        "schema_version": 1,
        "analysis_id": policy["analysis_id"],
        "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS",
        "local_results_accessed_before_input_freeze": False,
        "policy_sha256": sha256(POLICY),
        "pair_manifest_sha256": sha256(PAIR_MANIFEST),
        "locus_file_sha256": sha256(LAVA_LOCUS),
        "input_artifact_sha256": immutable_hashes,
        "overlap_diagnostics": overlap_diagnostics,
        "reference_status": "READY" if REFERENCE_PROVENANCE.is_file() else "BLOCKED_BY_DATA",
        "required_reference_provenance": str(REFERENCE_PROVENANCE),
        "reference_minimum_free_bytes": policy["reference_minimum_free_bytes"],
        "current_compute_snapshot": {
            "free_bytes": free_bytes,
            "status": "PASS" if free_bytes >= int(policy["reference_minimum_free_bytes"]) else "BLOCKED_BY_COMPUTE",
            "note": "Operational snapshot excluded from immutable artifact hashes; rerun to refresh after cleanup or handoff.",
        },
        "result_substitution_policy": "FORBIDDEN_SYNTHETIC_OR_PARTIAL_OUTPUTS_CANNOT_BE_PROMOTED_TO_REAL_LOCAL_RESULTS",
    }
    outputs[LOCK] = json.dumps(lock, indent=2, sort_keys=True) + "\n"
    return outputs


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build()
    if args.verify:
        # Operational free-space and handoff snapshots may change. Refresh them, then
        # require every immutable table to match the current scientific inputs.
        for path, value in expected.items():
            if path in {HANDOFF, LOCK}:
                continue
            if not path.is_file() or path.read_text(encoding="utf-8") != value:
                raise SystemExit(f"ERROR: Track B local input artifact missing or drifted: {path}")
        observed = json.loads(LOCK.read_text(encoding="utf-8"))
        for key in (
            "analysis_id", "selection_timing", "local_results_accessed_before_input_freeze",
            "policy_sha256", "pair_manifest_sha256", "locus_file_sha256",
            "input_artifact_sha256", "overlap_diagnostics", "result_substitution_policy",
        ):
            expected_lock = json.loads(expected[LOCK])
            if observed.get(key) != expected_lock.get(key):
                raise SystemExit(f"ERROR: Track B local input lock drifted: {key}")
        print("verified Track B local inputs: 3 pairs, 8 traits, 2495 loci")
        return
    for path, value in expected.items():
        atomic_text(path, value)
    print("wrote Track B local inputs: 3 pairs, 8 traits, 2495 loci; reference status blocked")


if __name__ == "__main__":
    main()
