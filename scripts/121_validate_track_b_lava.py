#!/usr/bin/env python3
"""Validate, summarize, and transactionally publish the exact Track B LAVA family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import math
import os
import shutil
import stat
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("track_b_lava_contract", ROOT / "scripts/119_track_b_lava_contract.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load Track B LAVA contract")
contract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contract)

NA_VALUES = {"", "NA", "NaN", "nan"}
RHO_Z_975 = 1.959963984540054
LOCAL_FIELDS = [
    "pair_id", "chr", "start", "end", "locus_id", "trait1", "trait2",
    "local_h2_trait1", "local_h2_trait2", "local_h2_scale", "local_covariance", "local_rg", "SE",
    "SE_method",
    "local_rg_CI_lower", "local_rg_CI_upper", "P", "FDR", "direction", "SNP_count",
    "QC", "interpretation_status",
]
STATUS_FIELDS = [
    "LOC", "CHR", "START", "STOP", "status", "n_snps", "K", "univariate_tested",
    "eligible_bivariate_pairs", "elapsed_seconds", "analysis_fingerprint", "lava_version", "reference_prefix",
]
UNIV_FIELDS = [
    "LOC", "CHR", "START", "STOP", "phen", "h2.obs", "h2.latent", "ascertained", "p",
    "analysis_status", "error", "n_snps", "K", "univariate_test_family_n",
    "univariate_p_threshold", "p_bonferroni", "p_fdr",
]
BIVAR_FIELDS = [
    "LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2", "discovery_rg",
    "discovery_SE", "discovery_P", "discovery_FDR", "local_covariance", "rho", "rho.lower",
    "rho.upper", "r2", "r2.lower", "r2.upper", "p", "analysis_status", "error",
    "bivariate_test_family_n", "p_fdr", "fdr_significant",
]
RAW_CONDITIONAL_FIELDS = [
    "LOC", "CHR", "START", "STOP", "pair_id", "trait1", "trait2", "conditional_model_id",
    "covariates", "pcor", "ci.lower", "ci.upper", "p", "r2.trait1_z", "r2.trait2_z",
    "analysis_status", "error", "conditional_test_family_n", "p_fdr", "fdr_significant",
]
CONDITIONAL_FIELDS = [
    "pair_id", "chr", "start", "end", "locus_id", "trait1", "trait2",
    "conditional_model_id", "covariates", "unconditioned_local_rg", "conditioned_local_rg",
    "conditioned_local_rg_CI_lower", "conditioned_local_rg_CI_upper",
    "absolute_attenuation_fraction", "P", "FDR", "direction",
    "r2_trait1_conditioners", "r2_trait2_conditioners", "QC", "interpretation_status",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        try:
            relative = path.relative_to(ROOT)
        except ValueError:
            relative = path
        fail(f"missing non-empty Track B LAVA result: {relative}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def atomic_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists():
            if artifact_identity(path) != artifact_identity(temporary):
                fail(f"existing staged summary differs; overwrite is forbidden: {path.relative_to(ROOT)}")
            return
        try:
            os.link(temporary, path)
        except FileExistsError:
            if artifact_identity(path) != artifact_identity(temporary):
                fail(f"staged summary appeared with different content: {path.relative_to(ROOT)}")
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)
        fsync_directory(path.parent)


def artifact_identity(path: Path) -> tuple[int, str]:
    try:
        label = path.relative_to(ROOT)
    except ValueError:
        label = path
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size == 0:
                fail(f"missing non-empty staged artifact: {label}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        fail(f"could not hash staged artifact {label}: {error}")

    def identity_fields(value: os.stat_result) -> tuple[int, int, int, int, int]:
        return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns

    if identity_fields(before) != identity_fields(after) or identity_fields(after) != identity_fields(current):
        fail(f"staged artifact changed while hashing: {label}")
    return after.st_size, digest.hexdigest()


def number(value: str, label: str) -> float:
    try:
        observed = float(value)
    except ValueError:
        fail(f"invalid numeric {label}: {value}")
    if not math.isfinite(observed):
        fail(f"non-finite numeric {label}: {value}")
    return observed


def nullable_number(value: str, label: str) -> float | None:
    if value in NA_VALUES:
        return None
    return number(value, label)


def require_na(value: str, label: str) -> None:
    if value not in NA_VALUES:
        fail(f"{label} must be NA when no valid test exists")


def probability(value: str, label: str = "probability") -> float:
    observed = number(value, label)
    if not 0 <= observed <= 1:
        fail(f"{label} outside [0,1]: {value}")
    return observed


def bh(values: list[float]) -> list[float]:
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [0.0] * count
    running = 1.0
    for rank0 in range(count - 1, -1, -1):
        index = order[rank0]
        running = min(running, min(1.0, values[index] * count / (rank0 + 1)))
        adjusted[index] = running
    return adjusted


def boolean(value: str, label: str = "logical value") -> bool:
    if value.lower() not in {"true", "false"}:
        fail(f"invalid {label}: {value}")
    return value.lower() == "true"


def close(observed: str, wanted: float, label: str) -> None:
    if not math.isclose(number(observed, label), wanted, rel_tol=1e-8, abs_tol=1e-12):
        fail(f"incorrect {label}: observed={observed} expected={wanted}")


def fmt(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.12g}"


def direction(value: float | None) -> str:
    if value is None:
        return "NOT_ESTIMATED"
    if value > 0:
        return "POSITIVE"
    if value < 0:
        return "NEGATIVE"
    return "ZERO"


def load_loci(policy: dict[str, Any]) -> tuple[list[str], dict[str, dict[str, str]]]:
    path = contract.LOCUS_FILE
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not lines:
        fail("LAVA locus definition is empty")
    header = lines[0]
    if not {"LOC", "CHR", "START", "STOP"}.issubset(header):
        fail("LAVA locus definition lacks required coordinates")
    rows = [dict(zip(header, values, strict=True)) for values in lines[1:]]
    loci = [row["LOC"] for row in rows]
    if len(loci) != policy["expected_loci"] or len(set(loci)) != len(loci):
        fail("LAVA locus family drifted")
    return loci, {row["LOC"]: row for row in rows}


def exact_fields(fields: list[str], expected: list[str], label: str) -> None:
    if fields != expected:
        fail(f"{label} result schema drifted: observed={fields} expected={expected}")


def validate_coordinates(row: dict[str, str], locus_index: dict[str, dict[str, str]], label: str) -> None:
    if row["LOC"] not in locus_index:
        fail(f"unknown locus in {label}: {row['LOC']}")
    expected = locus_index[row["LOC"]]
    if (row["CHR"], row["START"], row["STOP"]) != (expected["CHR"], expected["START"], expected["STOP"]):
        fail(f"{label} coordinates drifted: {row['LOC']}")


def validate_low_level(paths: list[Path]) -> dict[str, Any]:
    policy = contract.load_policy()
    _, input_info = read_tsv(ROOT / "results/track_b/lava_input_info.tsv")
    _, pairs = read_tsv(ROOT / "results/track_b/lava_pair_manifest.tsv")
    _, models = read_tsv(ROOT / "results/track_b/local_conditional_manifest.tsv")
    traits = [row["phenotype"] for row in input_info]
    binary_traits = {row["phenotype"] for row in input_info if row["cases"] not in NA_VALUES}
    if traits != policy["trait_order"] or [row["pair_id"] for row in pairs] != policy["pair_order"]:
        fail("Track B LAVA result inputs differ from policy order")
    expected_model_ids = [
        model["conditional_model_id"]
        for pair_id in policy["pair_order"]
        for model in policy["conditional_models"][pair_id]
    ] + ["NONE"]
    if [row["conditional_model_id"] for row in models] != expected_model_ids:
        fail("Track B conditional model manifest differs from policy")

    loci, locus_index = load_loci(policy)
    status_fields, status = read_tsv(paths[0])
    univ_fields, univ = read_tsv(paths[1])
    bivar_fields, bivar = read_tsv(paths[2])
    conditional_fields, conditional = read_tsv(paths[3])
    exact_fields(status_fields, STATUS_FIELDS, "status")
    exact_fields(univ_fields, UNIV_FIELDS, "univariate")
    exact_fields(bivar_fields, BIVAR_FIELDS, "bivariate")
    exact_fields(conditional_fields, RAW_CONDITIONAL_FIELDS, "conditional")

    fingerprint = contract.run_fingerprint()
    if len(status) != len(loci) or [row["LOC"] for row in status] != loci:
        fail("locus status does not cover the ordered 2,495-locus family")
    if {row["analysis_fingerprint"] for row in status} != {fingerprint}:
        fail("locus status fingerprint differs from current Track B contract")
    if {row["lava_version"] for row in status} != {policy["lava_version"]}:
        fail("locus status has wrong LAVA version")
    if {row["reference_prefix"] for row in status} != {policy["reference_prefix"]}:
        fail("locus status has wrong LAVA reference prefix")
    allowed_locus_status = {"PROCESSED", "PROCESS_FAILED", "UNIVARIATE_FAILED"}
    for row in status:
        if row["status"] not in allowed_locus_status:
            fail(f"invalid locus status: {row['status']}")
        validate_coordinates(row, locus_index, "locus status")
    failed_loci = sum(row["status"] != "PROCESSED" for row in status)
    if failed_loci / len(status) > policy["maximum_locus_failure_fraction"]:
        fail("locus failure fraction exceeds policy")

    expected_univ = [(locus, trait) for locus in loci for trait in traits]
    observed_univ = [(row["LOC"], row["phen"]) for row in univ]
    if len(univ) != policy["planned_univariate_tests"] or observed_univ != expected_univ:
        fail("univariate result is not the exact 2,495 x 8 family")
    if {int(row["univariate_test_family_n"]) for row in univ} != {policy["planned_univariate_tests"]}:
        fail("univariate family-size provenance drifted")
    threshold = float(policy["univariate_p_threshold"])
    allowed_univ_status = {"TESTED", "PHENOTYPE_DROPPED", "LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED"}
    univ_family_p: list[float] = []
    tested_univ: list[dict[str, str]] = []
    for row in univ:
        validate_coordinates(row, locus_index, "univariate result")
        if row["analysis_status"] not in allowed_univ_status:
            fail(f"invalid univariate status: {row['analysis_status']}")
        close(row["univariate_p_threshold"], threshold, "univariate threshold")
        if row["analysis_status"] == "TESTED":
            p_value = probability(row["p"], "univariate P")
            h2_observed = number(row["h2.obs"], "local h2 observed")
            if h2_observed < 0:
                fail(f"negative observed-scale local h2: {row['LOC']}/{row['phen']}")
            if row["phen"] in binary_traits:
                if number(row["h2.latent"], "local h2 latent") < 0:
                    fail(f"negative latent-scale local h2: {row['LOC']}/{row['phen']}")
            else:
                require_na(row["h2.latent"], "continuous-trait latent-scale local h2")
            close(row["p_bonferroni"], min(1.0, p_value * len(univ)), "univariate Bonferroni value")
            univ_family_p.append(p_value)
            tested_univ.append(row)
            require_na(row["error"], "tested univariate error")
        else:
            for field in ("h2.obs", "h2.latent", "p"):
                require_na(row[field], f"untested univariate {field}")
            if not row.get("error", "").strip() or row["ascertained"] not in NA_VALUES:
                fail("untested univariate row lacks an error or contains an ascertainment result")
            require_na(row["p_fdr"], "untested univariate FDR")
            require_na(row["p_bonferroni"], "untested univariate Bonferroni value")
            univ_family_p.append(1.0)
    if (len(univ) - len(tested_univ)) / len(univ) > policy["maximum_univariate_untested_fraction"]:
        fail("univariate untested fraction exceeds policy")
    univ_adjusted = bh(univ_family_p)
    for index, row in enumerate(univ):
        if row["analysis_status"] == "TESTED":
            close(row["p_fdr"], univ_adjusted[index], f"univariate FDR {row['LOC']}/{row['phen']}")

    significant_traits: dict[str, set[str]] = {locus: set() for locus in loci}
    for row in tested_univ:
        if probability(row["p"], "univariate P") <= threshold:
            significant_traits[row["LOC"]].add(row["phen"])
    pair_defs = [(row["pair_id"], row["trait1"], row["trait2"]) for row in pairs]
    expected_bivar = [
        (locus, pair_id, first, second)
        for locus in loci
        for pair_id, first, second in pair_defs
        if first in significant_traits[locus] and second in significant_traits[locus]
    ]
    observed_bivar = [(row["LOC"], row["pair_id"], row["trait1"], row["trait2"]) for row in bivar]
    if observed_bivar != expected_bivar:
        fail("bivariate result differs from the locally eligible frozen three-pair family")
    univ_by_locus = {locus: [row for row in univ if row["LOC"] == locus] for locus in loci}
    bivar_count_by_locus = {locus: sum(row["LOC"] == locus for row in bivar) for locus in loci}
    for status_row in status:
        locus = status_row["LOC"]
        locus_univ = univ_by_locus[locus]
        if int(status_row["univariate_tested"]) != sum(row["analysis_status"] == "TESTED" for row in locus_univ):
            fail(f"locus univariate-tested count drifted: {locus}")
        if int(status_row["eligible_bivariate_pairs"]) != bivar_count_by_locus[locus]:
            fail(f"locus eligible-bivariate count drifted: {locus}")
        if status_row["status"] == "PROCESSED":
            n_snps, components = int(status_row["n_snps"]), int(status_row["K"])
            if n_snps <= 0 or components <= 0:
                fail(f"processed locus has nonpositive SNP/component count: {locus}")
            if any(row["n_snps"] != status_row["n_snps"] or row["K"] != status_row["K"] for row in locus_univ):
                fail(f"univariate SNP/component counts differ from locus status: {locus}")
            if any(row["analysis_status"] in {"LOCUS_PROCESS_FAILED", "UNIVARIATE_FAILED"} for row in locus_univ):
                fail(f"processed locus contains locus-level failure rows: {locus}")
        elif status_row["status"] == "PROCESS_FAILED":
            if (
                any(row["analysis_status"] != "LOCUS_PROCESS_FAILED" for row in locus_univ)
                or bivar_count_by_locus[locus] != 0
                or status_row["n_snps"] not in NA_VALUES
                or status_row["K"] not in NA_VALUES
            ):
                fail(f"process-failed locus has incompatible scientific rows: {locus}")
        else:
            if (
                any(row["analysis_status"] != "UNIVARIATE_FAILED" for row in locus_univ)
                or bivar_count_by_locus[locus] != 0
                or int(status_row["n_snps"]) <= 0
                or int(status_row["K"]) <= 0
                or any(row["n_snps"] != status_row["n_snps"] or row["K"] != status_row["K"] for row in locus_univ)
            ):
                fail(f"univariate-failed locus has incompatible scientific rows: {locus}")
    if bivar and {int(row["bivariate_test_family_n"]) for row in bivar} != {len(bivar)}:
        fail("bivariate test-family size drifted")
    bivar_family_p: list[float] = []
    tested_bivar: list[dict[str, str]] = []
    pair_index = {row["pair_id"]: row for row in pairs}
    for row in bivar:
        validate_coordinates(row, locus_index, "bivariate result")
        frozen_pair = pair_index[row["pair_id"]]
        for observed_field, frozen_field in (
            ("discovery_rg", "discovery_rg"), ("discovery_SE", "discovery_SE"),
            ("discovery_P", "discovery_P"), ("discovery_FDR", "discovery_FDR"),
        ):
            if not math.isclose(
                number(row[observed_field], observed_field), number(frozen_pair[frozen_field], frozen_field),
                rel_tol=1e-12, abs_tol=1e-15,
            ):
                fail(f"bivariate discovery statistic drifted: {row['LOC']}/{row['pair_id']}/{observed_field}")
        status_value = row["analysis_status"]
        if status_value == "TESTED":
            require_na(row["error"], "tested bivariate error")
            rho = number(row["rho"], "local rg")
            rho_lower = number(row["rho.lower"], "local rg CI lower")
            rho_upper = number(row["rho.upper"], "local rg CI upper")
            r2 = number(row["r2"], "local r2")
            r2_lower = number(row["r2.lower"], "local r2 CI lower")
            r2_upper = number(row["r2.upper"], "local r2 CI upper")
            number(row["local_covariance"], "local covariance")
            if not -1 <= rho <= 1 or rho_lower > rho or rho_upper < rho:
                fail(f"invalid local rg or CI: {row['LOC']}/{row['pair_id']}")
            if not 0 <= r2 <= 1 or r2_lower > r2 or r2_upper < r2:
                fail(f"invalid local r2 or CI: {row['LOC']}/{row['pair_id']}")
            p_value = probability(row["p"], "bivariate P")
            bivar_family_p.append(p_value)
            tested_bivar.append(row)
        elif status_value == "BIVARIATE_FAILED":
            diagnostics = {
                field: nullable_number(row[field], f"failed bivariate {field}")
                for field in ("local_covariance", "rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper", "p")
            }
            if diagnostics["rho"] is not None and not -1 <= diagnostics["rho"] <= 1:
                fail("failed bivariate rho is outside [-1,1]")
            if diagnostics["r2"] is not None and not 0 <= diagnostics["r2"] <= 1:
                fail("failed bivariate r2 is outside [0,1]")
            if diagnostics["p"] is not None and not 0 <= diagnostics["p"] <= 1:
                fail("failed bivariate P is outside [0,1]")
            if all(diagnostics[field] is not None for field in ("rho", "rho.lower", "rho.upper")) and not diagnostics["rho.lower"] <= diagnostics["rho"] <= diagnostics["rho.upper"]:
                fail("failed bivariate retained rho/CI diagnostics are inconsistent")
            if all(diagnostics[field] is not None for field in ("r2", "r2.lower", "r2.upper")) and not diagnostics["r2.lower"] <= diagnostics["r2"] <= diagnostics["r2.upper"]:
                fail("failed bivariate retained r2/CI diagnostics are inconsistent")
            require_na(row["p_fdr"], "failed bivariate FDR")
            if not row.get("error", "").strip() or boolean(row["fdr_significant"], "bivariate significance"):
                fail("failed bivariate row lacks an error or is marked significant")
            bivar_family_p.append(1.0)
        else:
            fail(f"invalid bivariate status: {status_value}")
    if bivar and (len(bivar) - len(tested_bivar)) / len(bivar) > policy["maximum_bivariate_failure_fraction"]:
        fail("bivariate failure fraction exceeds policy")
    bivar_adjusted = bh(bivar_family_p)
    for index, row in enumerate(bivar):
        if row["analysis_status"] == "TESTED":
            close(row["p_fdr"], bivar_adjusted[index], f"bivariate FDR {row['LOC']}/{row['pair_id']}")
            if boolean(row["fdr_significant"], "bivariate significance") != (bivar_adjusted[index] <= policy["bivariate_fdr_alpha"]):
                fail("bivariate significance label differs from FDR")

    model_defs: dict[str, list[dict[str, str]]] = {
        pair_id: [row for row in models if row["pair_id"] == pair_id and row["conditional_model_id"] != "NONE"]
        for pair_id in ("A", "B")
    }
    expected_conditional = [
        (row["LOC"], row["pair_id"], row["trait1"], row["trait2"], model["conditional_model_id"])
        for row in bivar
        if row["analysis_status"] == "TESTED" and boolean(row["fdr_significant"], "bivariate significance") and row["pair_id"] in {"A", "B"}
        for model in model_defs[row["pair_id"]]
    ]
    observed_conditional = [
        (row["LOC"], row["pair_id"], row["trait1"], row["trait2"], row["conditional_model_id"])
        for row in conditional
    ]
    if observed_conditional != expected_conditional:
        fail("conditional result does not exactly cover FDR-supported A/B loci and frozen models")
    univ_index = {(row["LOC"], row["phen"]): row for row in univ}
    model_index = {row["conditional_model_id"]: row for row in models if row["conditional_model_id"] != "NONE"}
    tested_conditional: list[dict[str, str]] = []
    failed_conditional: list[dict[str, str]] = []
    eligible_conditional: list[dict[str, str]] = []
    conditional_family_p: list[float] = []
    max_r2 = float(policy["runtime"]["conditional_max_r2"])
    for row in conditional:
        validate_coordinates(row, locus_index, "conditional result")
        model = model_index.get(row["conditional_model_id"])
        if model is None or row["covariates"] != model["covariates"] or row["pair_id"] != model["pair_id"]:
            fail("conditional model identity or covariate set drifted")
        covariates = model["covariates"].split(";")
        eligible = all(
            (source := univ_index[(row["LOC"], covariate)])["analysis_status"] == "TESTED"
            and probability(source["p"], "conditioner univariate P") <= threshold
            for covariate in covariates
        )
        status_value = row["analysis_status"]
        if not eligible:
            if status_value != "CONDITIONER_LOCAL_H2_INELIGIBLE":
                fail("ineligible conditioner model was analyzed")
            if not row.get("error", "").strip():
                fail("conditioner-ineligible conditional row lacks an explanation")
            for field in ("pcor", "ci.lower", "ci.upper", "p", "r2.trait1_z", "r2.trait2_z", "p_fdr"):
                require_na(row[field], f"conditioner-ineligible conditional {field}")
        elif status_value == "TESTED":
            require_na(row["error"], "tested conditional error")
            pcor = number(row["pcor"], "partial local rg")
            ci_lower = number(row["ci.lower"], "partial local rg CI lower")
            ci_upper = number(row["ci.upper"], "partial local rg CI upper")
            r2_first = number(row["r2.trait1_z"], "conditional trait1 r2")
            r2_second = number(row["r2.trait2_z"], "conditional trait2 r2")
            if not -1 <= pcor <= 1 or ci_lower > pcor or ci_upper < pcor:
                fail(f"invalid partial local rg or CI: {row['LOC']}/{row['conditional_model_id']}")
            if not 0 <= r2_first < max_r2 or not 0 <= r2_second < max_r2:
                fail(f"invalid conditional r2: {row['LOC']}/{row['conditional_model_id']}")
            p_value = probability(row["p"], "conditional P")
            tested_conditional.append(row)
            eligible_conditional.append(row)
            conditional_family_p.append(p_value)
        elif status_value in {"CONDITIONAL_FAILED", "CONDITIONAL_UNSTABLE_MAX_R2"}:
            if not row.get("error", "").strip():
                fail("failed conditional row lacks an error")
            diagnostics = {
                field: nullable_number(row[field], f"failed conditional {field}")
                for field in ("pcor", "ci.lower", "ci.upper", "p", "r2.trait1_z", "r2.trait2_z")
            }
            if diagnostics["pcor"] is not None and not -1 <= diagnostics["pcor"] <= 1:
                fail("failed conditional pcor is outside [-1,1]")
            if diagnostics["p"] is not None and not 0 <= diagnostics["p"] <= 1:
                fail("failed conditional P is outside [0,1]")
            for field in ("r2.trait1_z", "r2.trait2_z"):
                if diagnostics[field] is not None and not 0 <= diagnostics[field] <= 1:
                    fail(f"failed conditional {field} is outside [0,1]")
            if status_value == "CONDITIONAL_UNSTABLE_MAX_R2":
                required = [diagnostics[field] for field in ("pcor", "ci.lower", "ci.upper", "r2.trait1_z", "r2.trait2_z")]
                if any(value is None for value in required) or max(diagnostics["r2.trait1_z"], diagnostics["r2.trait2_z"]) < max_r2:
                    fail("max-r2 unstable conditional row lacks the required retained diagnostics")
                require_na(row["p"], "max-r2 unstable conditional P")
            failed_conditional.append(row)
            eligible_conditional.append(row)
            conditional_family_p.append(1.0)
        else:
            fail("eligible conditional locus has invalid terminal status")
        if status_value != "TESTED":
            require_na(row["p_fdr"], "untested conditional FDR")
            if boolean(row["fdr_significant"], "conditional significance"):
                fail("untested conditional row is marked significant")
    denominator = len(eligible_conditional)
    if denominator and len(failed_conditional) / denominator > policy["maximum_conditional_failure_fraction"]:
        fail("conditional failure fraction exceeds policy")
    conditional_adjusted = bh(conditional_family_p)
    adjusted_by_identity = {
        (row["LOC"], row["pair_id"], row["conditional_model_id"]): adjusted
        for row, adjusted in zip(eligible_conditional, conditional_adjusted, strict=True)
    }
    for row in conditional:
        if int(row["conditional_test_family_n"]) != len(eligible_conditional):
            fail("conditional test-family size drifted")
        if row["analysis_status"] == "TESTED":
            wanted = adjusted_by_identity[(row["LOC"], row["pair_id"], row["conditional_model_id"])]
            close(row["p_fdr"], wanted, f"conditional FDR {row['LOC']}/{row['conditional_model_id']}")
            if boolean(row["fdr_significant"], "conditional significance") != (wanted <= policy["bivariate_fdr_alpha"]):
                fail("conditional significance label differs from FDR")

    return {
        "policy": policy, "loci": loci, "locus_index": locus_index, "traits": traits,
        "pairs": pairs, "models": models, "status": status, "univ": univ, "bivar": bivar,
        "conditional": conditional, "univ_index": univ_index,
        "bivar_index": {(row["LOC"], row["pair_id"]): row for row in bivar},
        "status_index": {row["LOC"]: row for row in status},
    }


def local_interpretation(row: dict[str, str] | None, significant_by_pair: dict[str, list[dict[str, str]]], pair_id: str) -> str:
    if row is None or row["analysis_status"] != "TESTED":
        return "UNDERPOWERED"
    if not boolean(row["fdr_significant"], "bivariate significance"):
        return "NO_LOCAL_SHARING"
    significant = significant_by_pair[pair_id]
    signs = {direction(number(item["rho"], "local rg")) for item in significant}
    if {"POSITIVE", "NEGATIVE"}.issubset(signs):
        return "OPPOSING_LOCAL_EFFECTS"
    if len(significant) >= 2:
        return "MULTIPLE_LOCAL_LOCI"
    return "STRONG_LOCAL_SHARING"


def conditional_interpretation(
    status_value: str, rho: float | None, pcor: float | None, fdr_significant: bool,
    partial_threshold: float, full_threshold: float,
) -> tuple[str, float | None, str]:
    if status_value == "CONDITIONER_LOCAL_H2_INELIGIBLE":
        return "UNDERPOWERED", None, "CONDITIONER_LOCAL_H2_INELIGIBLE"
    if status_value != "TESTED" or rho is None or pcor is None or rho == 0:
        attenuation = None if rho is None or pcor is None or rho == 0 else 1.0 - abs(pcor) / abs(rho)
        return "UNSTABLE", attenuation, status_value
    attenuation = 1.0 - abs(pcor) / abs(rho)
    if rho * pcor < 0:
        return "UNSTABLE", attenuation, "SIGN_REVERSAL"
    if fdr_significant and attenuation < partial_threshold:
        return "INDEPENDENT_LOCAL_COMPONENT", attenuation, "PASS"
    if not fdr_significant and attenuation >= full_threshold:
        return "FULLY_ATTENUATED", attenuation, "PASS"
    return "PARTIALLY_ATTENUATED", attenuation, "PASS"


def build_summaries(context: dict[str, Any]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    policy = context["policy"]
    significant_by_pair = {
        pair["pair_id"]: [
            row for row in context["bivar"]
            if row["pair_id"] == pair["pair_id"] and row["analysis_status"] == "TESTED"
            and boolean(row["fdr_significant"], "bivariate significance")
        ]
        for pair in context["pairs"]
    }
    local_rows: list[dict[str, str]] = []
    for pair in context["pairs"]:
        pair_id, trait1, trait2 = pair["pair_id"], pair["trait1"], pair["trait2"]
        for locus in context["loci"]:
            coordinates = context["locus_index"][locus]
            first = context["univ_index"][(locus, trait1)]
            second = context["univ_index"][(locus, trait2)]
            status = context["status_index"][locus]
            local = context["bivar_index"].get((locus, pair_id))
            tested = local is not None and local["analysis_status"] == "TESTED"
            rho = number(local["rho"], "local rg") if tested else None
            rho_lower = number(local["rho.lower"], "local rg CI lower") if tested else None
            rho_upper = number(local["rho.upper"], "local rg CI upper") if tested else None
            standard_error = (rho_upper - rho_lower) / (2 * RHO_Z_975) if tested else None
            if status["status"] != "PROCESSED":
                qc = f"LOCUS_{status['status']}"
            elif local is None:
                qc = "INSUFFICIENT_LOCAL_H2"
            elif not tested:
                qc = "BIVARIATE_FAILED"
            else:
                qc = "PASS"
            local_rows.append({
                "pair_id": pair_id, "chr": coordinates["CHR"], "start": coordinates["START"],
                "end": coordinates["STOP"], "locus_id": locus, "trait1": trait1, "trait2": trait2,
                "local_h2_trait1": first["h2.obs"] if first["analysis_status"] == "TESTED" else "NA",
                "local_h2_trait2": second["h2.obs"] if second["analysis_status"] == "TESTED" else "NA",
                "local_h2_scale": "OBSERVED",
                "local_covariance": local["local_covariance"] if tested else "NA",
                "local_rg": local["rho"] if tested else "NA", "SE": fmt(standard_error),
                "SE_method": "CI_EQUIVALENT_FROM_LAVA_95CI_WIDTH_DIV_3.919927969080108",
                "local_rg_CI_lower": local["rho.lower"] if tested else "NA",
                "local_rg_CI_upper": local["rho.upper"] if tested else "NA",
                "P": local["p"] if tested else "NA", "FDR": local["p_fdr"] if tested else "NA",
                "direction": direction(rho), "SNP_count": status["n_snps"], "QC": qc,
                "interpretation_status": local_interpretation(local, significant_by_pair, pair_id),
            })

    partial_threshold = float(policy["interpretation_policy"]["conditional_attenuation_fraction_partial"])
    full_threshold = float(policy["interpretation_policy"]["conditional_attenuation_fraction_full"])
    conditional_rows: list[dict[str, str]] = []
    for row in context["conditional"]:
        coordinates = context["locus_index"][row["LOC"]]
        local = context["bivar_index"][(row["LOC"], row["pair_id"])]
        rho = number(local["rho"], "unconditioned local rg")
        pcor = nullable_number(row["pcor"], "conditioned local rg")
        fdr_significant = boolean(row["fdr_significant"], "conditional significance")
        interpretation, attenuation, qc = conditional_interpretation(
            row["analysis_status"], rho, pcor, fdr_significant, partial_threshold, full_threshold,
        )
        conditional_rows.append({
            "pair_id": row["pair_id"], "chr": coordinates["CHR"], "start": coordinates["START"],
            "end": coordinates["STOP"], "locus_id": row["LOC"], "trait1": row["trait1"],
            "trait2": row["trait2"], "conditional_model_id": row["conditional_model_id"],
            "covariates": row["covariates"], "unconditioned_local_rg": local["rho"],
            "conditioned_local_rg": row["pcor"], "conditioned_local_rg_CI_lower": row["ci.lower"],
            "conditioned_local_rg_CI_upper": row["ci.upper"],
            "absolute_attenuation_fraction": fmt(attenuation), "P": row["p"], "FDR": row["p_fdr"],
            "direction": direction(pcor), "r2_trait1_conditioners": row["r2.trait1_z"],
            "r2_trait2_conditioners": row["r2.trait2_z"], "QC": qc,
            "interpretation_status": interpretation,
        })
    return local_rows, conditional_rows


def validate_summaries(paths: list[Path], expected: tuple[list[dict[str, str]], list[dict[str, str]]]) -> None:
    local_fields, local_rows = read_tsv(paths[0])
    conditional_fields, conditional_rows = read_tsv(paths[1])
    if local_fields != LOCAL_FIELDS or local_rows != expected[0]:
        fail("04_lava_local_results.tsv differs from the exact validated 3 x 2,495 summary")
    if conditional_fields != CONDITIONAL_FIELDS or conditional_rows != expected[1]:
        fail("05_local_conditional_results.tsv differs from the validated conditional model family")


def publish(staged: list[Path], validated_identity: dict[Path, tuple[int, str]] | None = None) -> None:
    for source in staged:
        if not source.is_file() or source.stat().st_size == 0:
            fail(f"staged Track B LAVA result is missing: {source.relative_to(ROOT)}")
    if validated_identity is None:
        validated_identity = {path: artifact_identity(path) for path in staged}
    temporaries: list[Path] = []
    try:
        for index, (source, destination) in enumerate(zip(staged, contract.RESULTS, strict=True), start=1):
            if artifact_identity(source) != validated_identity[source]:
                fail(f"staged Track B LAVA result changed after validation: {source.relative_to(ROOT)}")
            if destination.exists():
                if artifact_identity(destination) != validated_identity[source]:
                    fail(
                        "unsealed or sealed canonical output differs from validated staging; "
                        f"overwrite is forbidden: {destination.relative_to(ROOT)}"
                    )
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{index}.publishing")
            with source.open("rb") as source_handle, temporary.open("xb") as temporary_handle:
                shutil.copyfileobj(source_handle, temporary_handle, length=1024 * 1024)
                temporary_handle.flush()
                os.fsync(temporary_handle.fileno())
            temporaries.append(temporary)
            if artifact_identity(temporary) != validated_identity[source]:
                fail(f"Track B LAVA publication copy differs from validated stage: {source.relative_to(ROOT)}")
            temporary_stat = temporary.stat()
            try:
                os.link(temporary, destination)
            except FileExistsError:
                fail(f"canonical Track B LAVA output appeared during publication: {destination.relative_to(ROOT)}")
            destination_stat = destination.stat()
            if (destination_stat.st_dev, destination_stat.st_ino) != (temporary_stat.st_dev, temporary_stat.st_ino):
                fail(f"Track B LAVA no-replace link identity mismatch: {destination.relative_to(ROOT)}")
            fsync_directory(destination.parent)
            temporary.unlink()
            fsync_directory(destination.parent)
        for source, destination in zip(staged, contract.RESULTS, strict=True):
            if artifact_identity(source) != validated_identity[source] or artifact_identity(destination) != validated_identity[source]:
                fail("Track B LAVA artifact changed between validation and result sealing")
            fsync_directory(destination.parent)
        # The long execution may span hours or days. Fully revalidate the exact
        # chromosome projections and rehash the reference at the last possible point.
        contract.fully_verify_chromosome_inputs()
        contract.validate_reference(rehash=True)
        expected_canonical = {
            str(destination.relative_to(ROOT)): validated_identity[source]
            for source, destination in zip(staged, contract.RESULTS, strict=True)
        }
        contract.seal_results(expected_canonical)
    except BaseException:
        # Matching canonical files are intentionally retained. They are immutable,
        # independently hash-checked restart progress, not an ambiguous partial result.
        for path in temporaries:
            path.unlink(missing_ok=True)
            fsync_directory(path.parent)
        raise


def validate(paths: list[Path], *, build_summary_files: bool) -> dict[Path, tuple[int, str]]:
    low_level_before = {path: artifact_identity(path) for path in paths[:4]}
    context = validate_low_level(paths[:4])
    summaries = build_summaries(context)
    if build_summary_files:
        atomic_tsv(paths[4], LOCAL_FIELDS, summaries[0])
        atomic_tsv(paths[5], CONDITIONAL_FIELDS, summaries[1])
    validate_summaries(paths[4:], summaries)
    low_level_after = {path: artifact_identity(path) for path in paths[:4]}
    if low_level_before != low_level_after:
        fail("staged low-level Track B LAVA output changed during validation")
    return {path: artifact_identity(path) for path in paths}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--validate-staged-results", action="store_true")
    mode.add_argument("--verify-staged-results", action="store_true")
    mode.add_argument("--publish-staged-results", action="store_true")
    parser.add_argument("--staging-dir", type=Path)
    parser.add_argument("--fingerprint")
    args = parser.parse_args()
    if args.validate_staged_results or args.verify_staged_results or args.publish_staged_results:
        current_fingerprint = contract.run_fingerprint()
        fingerprint = args.fingerprint or current_fingerprint
        if fingerprint != current_fingerprint:
            fail("supplied staging fingerprint differs from the live Track B contract")
        if args.staging_dir is not None:
            directory = args.staging_dir.resolve(strict=True)
            attempts = (contract.CHECKPOINT_DIR / fingerprint / "finalize" / ".attempts").resolve()
            if directory.parent != attempts or directory.is_symlink() or not directory.name.startswith("unit_all."):
                fail("explicit staging directory is not the current supervised finalize attempt")
            staged = [directory / path.name for path in contract.RESULTS]
        else:
            staged = contract.staged_results(fingerprint)
        validated_identity = validate(staged, build_summary_files=args.validate_staged_results)
        if args.publish_staged_results:
            publish(staged, validated_identity)
            print("TRACK_B_LAVA_STAGED_RESULTS_VALIDATED_AND_PUBLISHED")
        elif args.verify_staged_results:
            print("TRACK_B_LAVA_STAGED_RESULTS_READ_ONLY_VERIFIED")
        else:
            print("TRACK_B_LAVA_STAGED_RESULTS_SEMANTICALLY_VALIDATED")
    else:
        validate(contract.RESULTS, build_summary_files=False)
        contract.validate_results()
        print("TRACK_B_LAVA_RESULTS_SEMANTICALLY_VALIDATED")


if __name__ == "__main__":
    main()
