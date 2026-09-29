#!/usr/bin/env python3
"""Materialize full LAVA h2/CI details after the frozen family has completed.

The production runner intentionally stores compact locus checkpoints. This
post-run step restores local h2 fields from their receipt-bound worker outputs
and replays only loci with a tested bivariate result to capture LAVA's native
confidence intervals. It refuses partial families and checks replayed rho/P
against the frozen result before publishing any detail table.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "extensions/brain6"
sys.path.insert(0, str(EXTENSION))
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from brain6.artifacts import verify_artifact  # noqa: E402
from brain6.io import read_json, sha256  # noqa: E402
from run_lava_family import EXPECTED_PAIRS, load_loci, validate_inputs, verify_unit  # noqa: E402

FAMILY_PATH = ROOT / "brain6/config/lava_family_v2.json"
EXECUTION_PATH = ROOT / "brain6/config/lava_execution_v1.json"
R_REPLAY = ROOT / "brain6/scripts/replay_lava_tested_bivariate.R"
ROUND_OFF_EXECUTION_PATH = ROOT / "brain6/config/lava_execution_roundoff_v1.json"
ROUND_OFF_RUNNER = ROOT / "brain6/scripts/run_lava_family_roundoff_v1.py"
ROUND_OFF_WORKER = ROOT / "brain6/scripts/run_lava_family_locus_roundoff_v1.R"
ROUND_OFF_REPLAY = ROOT / "brain6/scripts/replay_lava_tested_bivariate_roundoff_v1.R"
H2_FIELDS = ("h2.obs", "h2.latent", "ascertained", "p")


def resolve_run_variant(identity: dict) -> dict[str, Any]:
    """Bind post-run processing to one of the explicitly supported run/code locks."""
    common = {
        "family_lock_sha256": sha256(FAMILY_PATH),
        "family_aggregator_source_sha256": sha256(EXTENSION / "brain6/family24.py"),
        "runtime_validator_sha256": sha256(ROOT / "scripts/133_validate_track_b_lava_runtime.R"),
    }
    variants = [
        {"name": "baseline", "execution_path": EXECUTION_PATH,
         "runner_path": ROOT / "brain6/scripts/run_lava_family.py",
         "worker_path": ROOT / "brain6/scripts/run_lava_family_locus.R",
         "replay_path": R_REPLAY, "roundoff_patch": False,
         "default_output_dir": ROOT / "brain6/results/local/lava_family_v2_detailed"},
        {"name": "roundoff_v1", "execution_path": ROUND_OFF_EXECUTION_PATH,
         "runner_path": ROUND_OFF_RUNNER, "worker_path": ROUND_OFF_WORKER,
         "replay_path": ROUND_OFF_REPLAY, "roundoff_patch": True,
         "default_output_dir": ROOT / "brain6/results/local/lava_family_roundoff_v1_detailed"},
    ]
    for variant in variants:
        execution = read_json(variant["execution_path"])
        if variant["roundoff_patch"]:
            patch = execution.get("numerical_stability_patch", {})
            if (patch.get("id") != "LAVA015_BLOCK_REDUCED_SYMMETRY_V1"
                    or patch.get("worker_script_sha256") != sha256(variant["worker_path"])
                    or patch.get("scientific_parameters_changed") is not False):
                continue
        expected = {
            **common,
            "execution_lock_sha256": sha256(variant["execution_path"]),
            "orchestrator_script_sha256": sha256(variant["runner_path"]),
            "worker_sha256": sha256(variant["worker_path"]),
        }
        if all(identity.get(key) == value for key, value in expected.items()):
            return variant
    raise ValueError("LAVA run identity does not match a supported, checksum-locked worker/execution variant")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_gzip_tsv(path: Path, fields: list[str], data: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.partial")
    try:
        with tmp.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as zipped:
                import io
                text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(data)
                text.flush()
                text.detach()
        payload = tmp.read_bytes()
        if path.exists():
            if path.read_bytes() != payload:
                raise FileExistsError(f"Refusing to replace a different LAVA detail artifact: {path}")
            tmp.unlink()
        else:
            os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def discover_run(output_root: Path, run_id: str | None) -> Path:
    if run_id:
        candidates = [output_root / run_id]
    else:
        candidates = sorted(p for p in output_root.iterdir()
                            if p.is_dir() and (p / "family_decision.json").is_file())
    if len(candidates) != 1 or not candidates[0].is_dir():
        raise ValueError("Specify --run-id for exactly one completed LAVA run")
    return candidates[0]


def canonical_univariate(run_root: Path, loci: list[dict[str, str]], traits: list[str]) -> list[dict[str, Any]]:
    aggregate_rows = read_rows(run_root / "univariate_results.tsv")
    aggregate = {(r["trait_id"], r["locus_id"]): r for r in aggregate_rows}
    expected = {(trait, str(locus["LOC"])) for trait in traits for locus in loci}
    if len(aggregate_rows) != len(expected) or set(aggregate) != expected:
        raise ValueError("Aggregated unique local-h2 table does not cover the frozen family")
    observed: dict[tuple[str, str], dict[str, str]] = {}
    sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    for locus in loci:
        locus_id = str(locus["LOC"])
        checkpoint = run_root / "loci" / locus_id
        for raw in read_rows(checkpoint / "worker_output/univariate.tsv"):
            trait = raw.get("phen", "")
            key = (trait, locus_id)
            if trait not in traits:
                raise ValueError(f"Unexpected univariate phenotype: {key}")
            pair_id = raw.get("pair_id", "")
            if pair_id in sources[key]:
                raise ValueError(f"Duplicate local-h2 test for pair and trait/locus: {pair_id}/{key}")
            if key in observed:
                for field in H2_FIELDS:
                    left, right = observed[key].get(field, "NA"), raw.get(field, "NA")
                    if left != right:
                        try:
                            equal = math.isclose(float(left), float(right), rel_tol=0, abs_tol=1e-12)
                        except (TypeError, ValueError):
                            equal = left == right
                        if not equal:
                            raise ValueError(f"Pair-specific local-h2 details disagree at {key}/{field}")
            else:
                observed[key] = raw
            sources[key].add(pair_id)

    output = []
    for key in sorted(expected, key=lambda x: (int(x[1]), traits.index(x[0]))):
        row = aggregate[key]
        detail = observed.get(key)
        is_tested = row["status"] == "TESTED"
        if is_tested != (detail is not None):
            raise ValueError(f"Raw local-h2 output and aggregated status disagree at {key}")
        if detail is not None and not math.isclose(float(row["p"]), float(detail["p"]),
                                                   rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Aggregated local-h2 P differs from receipt output at {key}")
        output.append({
            "trait_id": key[0], "locus_id": key[1], "status": row["status"],
            "local_h2_observed": detail.get("h2.obs", "NA") if detail else "NA",
            "local_h2_latent": detail.get("h2.latent", "NA") if detail else "NA",
            "ascertained": detail.get("ascertained", "NA") if detail else "NA",
            "p_h2": detail.get("p", "NA") if detail else "NA",
            "se_local_h2": "NA", "se_status": "NOT_RETURNED_BY_PINNED_LAVA_RUN_UNIV",
            "source_pair_ids": ";".join(sorted(sources[key])),
            "reason": row.get("reason", ""),
        })
    return output


def verify_completed_run(run_root: Path, variant: dict[str, Any]) -> tuple[dict, dict, dict, list[dict[str, str]], list[dict[str, str]]]:
    family, execution = read_json(FAMILY_PATH), read_json(variant["execution_path"])
    identity = read_json(run_root / "run_identity.json")
    run_id = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if run_id != run_root.name:
        raise ValueError("Run directory does not match its source-bound identity")
    if (identity.get("family_lock_sha256") != sha256(FAMILY_PATH)
            or identity.get("execution_lock_sha256") != sha256(variant["execution_path"])
            or identity.get("orchestrator_script_sha256") != sha256(variant["runner_path"])
            or identity.get("worker_sha256") != sha256(variant["worker_path"])
            or identity.get("family_aggregator_source_sha256") != sha256(EXTENSION / "brain6/family24.py")
            or identity.get("runtime_validator_sha256") != sha256(ROOT / "scripts/133_validate_track_b_lava_runtime.R")):
        raise ValueError("LAVA run identity no longer matches the frozen family/execution locks")
    if list(family["pairs"]) != list(EXPECTED_PAIRS):
        raise ValueError("Frozen LAVA pair order differs from the production runner")
    decision_path = run_root / "family_decision.json"
    manifest_path = run_root / "local_family_manifest.json"
    if not decision_path.is_file() or not manifest_path.is_file():
        raise ValueError("LAVA family is not complete; detail materialization is post-run only")
    decision = read_json(decision_path)
    if decision.get("analysis_id") != family["analysis_id"]:
        raise ValueError("Unexpected LAVA family decision identity")
    manifest = read_json(manifest_path)
    result_path = Path(manifest["results"])
    planned_path = Path(manifest["planned_units"])
    if not result_path.is_file() or not planned_path.is_file():
        raise ValueError("LAVA aggregate input is missing")
    final_artifact = Path(decision["local_family_artifact"])
    artifact_receipt = verify_artifact(final_artifact)
    if artifact_receipt.get("scientific_status") != decision.get("local_family_status"):
        raise ValueError("LAVA local-family receipt status differs from decision")
    if sha256(final_artifact / "receipt.json") != decision.get("local_family_receipt_sha256"):
        raise ValueError("LAVA local-family receipt checksum differs from family decision")
    bound_inputs = {Path(r["path"]).resolve(): r for r in artifact_receipt.get("inputs", [])}
    for path in (manifest_path, result_path, planned_path):
        record = bound_inputs.get(path.resolve())
        if record is None or record["sha256"] != sha256(path) or record["bytes"] != path.stat().st_size:
            raise ValueError(f"LAVA family artifact does not bind source input: {path}")
    if (sha256(Path(manifest["worker_receipts_path"])) != manifest["worker_receipts_sha256"]
            or sha256(Path(manifest["univariate_results_path"])) != manifest["univariate_results_sha256"]):
        raise ValueError("LAVA receipt or univariate aggregate differs from its manifest checksum")
    correction_path = Path(decision["family_correction_path"])
    if not correction_path.is_file() or sha256(correction_path) != decision.get("family_correction_sha256"):
        raise ValueError("LAVA correction file differs from the family decision checksum")
    loci_path = ROOT / family["locus_definition"]["path"]
    loci = load_loci(loci_path, family["locus_definition"]["sha256"], int(family["locus_definition"]["n_loci"]))
    plan = read_rows(planned_path)
    expected_plan = [(pair, str(locus["LOC"])) for locus in loci for pair in family["pairs"]]
    if [(r["pair_id"], r["locus_id"]) for r in plan] != expected_plan:
        raise ValueError("LAVA planned-unit ledger differs from the exact frozen family")
    results = read_rows(result_path)
    expected_keys = set(expected_plan)
    if len(results) != len(expected_plan) or {(r["pair_id"], r["locus_id"]) for r in results} != expected_keys:
        raise ValueError("LAVA results do not preserve every frozen pair-by-locus slot")
    receipts = read_rows(Path(manifest["worker_receipts_path"]))
    if len(receipts) != len(loci):
        raise ValueError("Completed family does not have one verified receipt per frozen locus")
    receipt_by_locus = {r["locus_id"]: r for r in receipts}
    if set(receipt_by_locus) != {str(r["LOC"]) for r in loci}:
        raise ValueError("LAVA receipt index does not cover all frozen loci")
    for locus in loci:
        locus_id = str(locus["LOC"])
        unit = run_root / "loci" / locus_id
        receipt = verify_unit(unit, locus_id, set(family["pairs"]), run_id)
        indexed = receipt_by_locus[locus_id]
        if receipt is None or indexed["receipt_sha256"] != sha256(unit / "receipt.json"):
            raise ValueError(f"LAVA receipt index mismatch at locus {locus_id}")
    aggregate_by_key = {(r["pair_id"], r["locus_id"]): r for r in results}
    for locus in loci:
        locus_id = str(locus["LOC"])
        unit_rows = read_rows(run_root / "loci" / locus_id / "worker_output/pair_results.tsv")
        for unit_row in unit_rows:
            key = unit_row["pair_id"], locus_id
            aggregate_row = aggregate_by_key[key]
            if any(unit_row.get(field, "") != aggregate_row.get(field, "")
                   for field in ("pair_id", "locus_id", "status", "p", "local_rg", "reason")):
                raise ValueError(f"Aggregated LAVA result differs from receipt at {key}")
    corrections = read_rows(correction_path)
    if {(r["pair_id"], r["locus_id"]) for r in corrections} != expected_keys:
        raise ValueError("LAVA correction table does not cover every prespecified slot")
    raw_by_key = {(r["pair_id"], r["locus_id"]): r for r in results}
    local_family_rows = read_rows(final_artifact / "local_family.tsv")
    if {(r["pair_id"], r["locus_id"]) for r in local_family_rows} != expected_keys:
        raise ValueError("LAVA corrected artifact does not cover every prespecified slot")
    local_by_key = {(r["pair_id"], r["locus_id"]): r for r in local_family_rows}
    for row in corrections:
        key = row["pair_id"], row["locus_id"]
        source = raw_by_key[key]
        corrected = local_by_key[key]
        if row["status"] != source["status"]:
            raise ValueError("LAVA correction status differs from the source result table")
        tested = source["status"] == "TESTED"
        expected_p = float(source["p"]) if tested else 1.0
        if not math.isclose(float(row["p_for_family_correction"]), expected_p,
                            rel_tol=0, abs_tol=1e-12):
            raise ValueError("LAVA correction does not apply the frozen tested/non-tested P rule")
        expected_q = corrected["family_fdr"] if tested else "1.0"
        if not math.isclose(float(row["family_fdr"]), float(expected_q), rel_tol=0, abs_tol=1e-12):
            raise ValueError("LAVA exported family correction differs from the corrected artifact")
        if tested and (not math.isclose(float(corrected["p"]), float(source["p"]), rel_tol=0, abs_tol=1e-12)
                       or not math.isclose(float(corrected["local_rg"]), float(source["local_rg"]), rel_tol=0, abs_tol=1e-12)):
            raise ValueError("LAVA corrected artifact estimate differs from the raw result")
    return family, execution, decision, loci, results


def replay_bivariate(run_root: Path, family: dict, execution: dict,
                     results: list[dict[str, str]], input_root: Path,
                     reference: Path, rscript: Path, output_path: Path,
                     variant: dict[str, Any]) -> list[dict[str, str]]:
    tested = [r for r in results if r["status"] == "TESTED"]
    if not tested:
        return []
    pair_lock = read_json(ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json")
    pair_lock_body = {k: v for k, v in pair_lock.items() if k != "lock_sha256"}
    pair_lock_sha = hashlib.sha256(json.dumps(pair_lock_body, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()
    if pair_lock_sha != family["pair_lock_sha256"]:
        raise ValueError("Reviewed deep-pair lock differs from the frozen LAVA family")
    pair_map = {r["pair_id"]: r for r in pair_lock["pairs"]}
    materialized = read_json(input_root / "provenance.json")
    if sha256(input_root / "provenance.json") != read_json(run_root / "run_identity.json")["materialized_inputs_sha256"]:
        raise ValueError("Post-run detail inputs differ from the LAVA run's materialized-input identity")
    ref_provenance = reference / "reference.provenance.json"
    if materialized.get("official_reference_provenance_sha256") != sha256(ref_provenance):
        raise ValueError("Post-run detail replay reference differs from materialized LAVA inputs")
    coverage_path = Path(materialized["pair_coverage"]["path"])
    coverage_meta = materialized["pair_coverage"]
    if coverage_path.stat().st_size != int(coverage_meta["bytes"]) or sha256(coverage_path) != coverage_meta["sha256"]:
        raise ValueError("Pair-by-locus reference coverage differs from its frozen input manifest")
    coverage = {(r["pair_id"], r["locus_id"]): int(r["shared_reference_variants"])
                for r in read_rows(coverage_path)}
    locus_chromosome = {}
    with (ROOT / family["locus_definition"]["path"]).open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) == 4:
                locus_chromosome[fields[0]] = int(fields[1])
    pairs = []
    for result in tested:
        spec = pair_map[result["pair_id"]]
        matrix = next(x for x in family["overlap_gate"]["pairwise_sample_overlap_file_manifest"]["matrices"]
                      if x["pair_id"] == result["pair_id"])
        chrom = locus_chromosome.get(result["locus_id"])
        if chrom is None:
            raise ValueError(f"Locus chromosome is absent from frozen locus file: {result['locus_id']}")
        pairs.append({
            **result, "phenotypes": [spec["sleep_trait"], spec["disease_trait"]],
            "input_info": str(input_root / f"locus_{result['locus_id']}" / "input_info.tsv"),
            "sample_overlap_file": str(ROOT / matrix["path"]),
            "reference_prefix": str(reference / "lava-ukb-v1.1") + f"_chr{chrom}",
            "shared_reference_variants": coverage[(result["pair_id"], result["locus_id"])],
        })
    config = {
        "analysis_id": family["analysis_id"], "family_lock_sha256": sha256(FAMILY_PATH),
        "execution_lock_sha256": sha256(variant["execution_path"]),
        "roundoff_patch": variant["roundoff_patch"],
        "loci_file": str((ROOT / family["locus_definition"]["path"]).resolve()),
        "input_root": str(input_root.resolve()), "random_seed": int(family["execution"]["random_seed"]),
        "locus_processing": execution["locus_processing"],
        "univariate": {**execution["univariate"],
                       "gate_p_strictly_less_than": family["univariate_gate"]["p_threshold_strictly_less_than"]},
        "bivariate": execution["bivariate"],
        "pairs": pairs, "output_path": str(output_path.resolve()),
    }
    fd, temp = tempfile.mkstemp(prefix=".lava_detail_replay_", suffix=".json", dir=output_path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(config, f, sort_keys=True); f.write("\n")
        helper = variant["replay_path"]
        proc = subprocess.run([str(rscript), str(helper), temp], cwd=ROOT,
                              capture_output=True, text=True, check=False)
        if proc.returncode:
            raise RuntimeError("LAVA bivariate replay failed: " + (proc.stdout + proc.stderr)[-4000:])
    finally:
        Path(temp).unlink(missing_ok=True)
    return read_rows(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-id")
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=ROOT / "ref/lava/ukb_v1.1")
    parser.add_argument("--rscript", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path,
                        help="Optional release directory; defaults to a run-variant-specific path")
    args = parser.parse_args()
    run_root = discover_run(args.output_root.resolve(), args.run_id)
    identity = read_json(run_root / "run_identity.json")
    variant = resolve_run_variant(identity)
    family, execution, decision, loci, results = verify_completed_run(run_root, variant)
    input_root = args.input_root.resolve()
    input_provenance = input_root / "provenance.json"
    if not input_provenance.is_file() or sha256(input_provenance) != identity.get("materialized_inputs_sha256"):
        raise ValueError("Post-run detail inputs differ from the completed run's materialized-input identity")
    reference = args.reference.resolve()
    expected_reference = (ROOT / "ref/lava/ukb_v1.1").resolve()
    if reference != expected_reference:
        raise ValueError("Post-run replay must use the exact reference path used by the production runner")
    reference_provenance = reference / "reference.provenance.json"
    materialized_input_record = read_json(input_provenance)
    if (not reference_provenance.is_file()
            or materialized_input_record.get("official_reference_provenance_sha256") != sha256(reference_provenance)):
        raise ValueError("Post-run detail reference differs from the frozen materialized-input provenance")
    # Reuse the production validator so every locus-specific summary-stat file,
    # input-info record, overlap matrix, and coverage row is rechecked here too.
    validate_inputs(FAMILY_PATH, variant["execution_path"], input_root, reference)
    if not args.rscript.is_file() or sha256(args.rscript.resolve()) != identity.get("rscript_binary_sha256"):
        raise ValueError("Post-run LAVA replay runtime differs from the frozen binary")
    validator = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
    runtime = subprocess.run([str(args.rscript.resolve()), str(validator)], cwd=ROOT,
                             capture_output=True, text=True, check=False)
    runtime_sha256 = hashlib.sha256(runtime.stdout.encode("utf-8")).hexdigest()
    if (runtime.returncode or runtime.stdout.count("TRACK_B_LAVA_RUNTIME\tstatus=PASS\t") != 1
            or runtime_sha256 != identity.get("runtime_validation_output_sha256")):
        raise ValueError("Pinned LAVA runtime validation no longer matches the completed run")
    reference_check = subprocess.run([sys.executable, "scripts/lava_contract.py", "--verify-reference"],
                                     cwd=ROOT, capture_output=True, text=True, check=False)
    if reference_check.returncode or "LAVA_REFERENCE_VALIDATED" not in reference_check.stdout:
        raise ValueError("Official LAVA reference failed verification: " +
                         (reference_check.stdout + reference_check.stderr)[-4000:])
    if decision.get("overall_status") not in {"PASS", "INSUFFICIENT_EVIDENCE", "FAILED_QC_NOT_CONSUMED"}:
        raise ValueError("Unexpected LAVA family decision status")
    h2_rows = canonical_univariate(run_root, loci, list(family["trait_ids"]))
    output_dir = (args.output_dir.resolve() if args.output_dir else variant["default_output_dir"])
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    h2_path = output_dir / "local_h2_tests.tsv.gz"
    rg_path = output_dir / "local_rg_detailed.tsv.gz"
    replay_work = Path(tempfile.mkdtemp(prefix="brain6-lava-detail-replay-"))
    replay_path = replay_work / "bivariate.tsv"
    replay_sha256 = "NO_TESTED_SLOT_REPLAY_REQUIRED"
    try:
        replayed = replay_bivariate(run_root, family, execution, results, input_root,
                                    reference, args.rscript.resolve(), replay_path, variant)
        if replayed:
            replay_sha256 = sha256(replay_path)
    finally:
        import shutil
        shutil.rmtree(replay_work, ignore_errors=True)
    replay_by_key = {(r["pair_id"], r["locus_id"]): r for r in replayed}
    tested_keys = {(r["pair_id"], r["locus_id"]) for r in results if r["status"] == "TESTED"}
    if set(replay_by_key) != tested_keys:
        raise ValueError("LAVA CI replay does not cover exactly the frozen TESTED slots")
    correction_path = Path(read_json(run_root / "family_decision.json")["family_correction_path"])
    correction = {(r["pair_id"], r["locus_id"]): r for r in read_rows(correction_path)}
    detailed = []
    for row in results:
        key = row["pair_id"], row["locus_id"]
        replay = replay_by_key.get(key)
        if replay:
            for field in ("rho", "p"):
                if not math.isclose(float(row["local_rg"] if field == "rho" else row["p"]),
                                    float(replay[field]), rel_tol=0, abs_tol=1e-10):
                    raise ValueError(f"Rerun LAVA {field} differs from frozen production result at {key}")
        adj = correction[key]
        detailed.append({
            "pair_id": row["pair_id"], "locus_id": row["locus_id"], "status": row["status"],
            "local_rg": row["local_rg"], "se_local_rg": "NA",
            "se_status": "NOT_RETURNED_BY_PINNED_LAVA_RUN_BIVAR",
            "ci_lower": replay.get("rho.lower", "NA") if replay else "NA",
            "ci_upper": replay.get("rho.upper", "NA") if replay else "NA",
            "p": row["p"], "family_fdr": adj["family_fdr"],
            "p_for_family_correction": adj["p_for_family_correction"],
            "r2": replay.get("r2", "NA") if replay else "NA",
            "r2_ci_lower": replay.get("r2.lower", "NA") if replay else "NA",
            "r2_ci_upper": replay.get("r2.upper", "NA") if replay else "NA",
            "reason": row["reason"],
        })
    staging_dir = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", suffix=".partial",
                                        dir=output_dir.parent))
    staged_h2_path = staging_dir / h2_path.name
    staged_rg_path = staging_dir / rg_path.name
    published = False
    try:
        write_gzip_tsv(staged_h2_path, ["trait_id", "locus_id", "status", "local_h2_observed",
            "local_h2_latent", "ascertained", "p_h2", "se_local_h2", "se_status",
            "source_pair_ids", "reason"], h2_rows)
        write_gzip_tsv(staged_rg_path, ["pair_id", "locus_id", "status", "local_rg", "se_local_rg",
            "se_status", "ci_lower", "ci_upper", "p", "family_fdr", "p_for_family_correction",
            "r2", "r2_ci_lower", "r2_ci_upper", "reason"], detailed)
        provenance = {
            "schema_version": 1,
            "analysis_id": family["analysis_id"],
            "run_id": run_root.name,
            "run_variant": variant["name"],
            "overall_family_status": decision["overall_status"],
            "local_family_artifact": decision["local_family_artifact"],
            "local_family_receipt_sha256": decision["local_family_receipt_sha256"],
            "family_lock_sha256": sha256(FAMILY_PATH),
            "execution_lock_sha256": sha256(variant["execution_path"]),
            "roundoff_patch_applied_during_replay": variant["roundoff_patch"],
            "run_identity_sha256": sha256(run_root / "run_identity.json"),
            "production_materializer_script_sha256": identity.get("materializer_script_sha256"),
            "run_results_sha256": sha256(Path(read_json(run_root / "local_family_manifest.json")["results"])),
            "worker_receipts_sha256": sha256(Path(read_json(run_root / "local_family_manifest.json")["worker_receipts_path"])),
            "univariate_results_sha256": sha256(Path(read_json(run_root / "local_family_manifest.json")["univariate_results_path"])),
            "planned_units_sha256": sha256(Path(read_json(run_root / "local_family_manifest.json")["planned_units"])),
            "materialized_inputs_provenance_sha256": sha256(input_provenance),
            "official_reference_provenance_sha256": sha256(reference_provenance),
            "replay_script_sha256": sha256(variant["replay_path"]),
            "bivariate_replay_tsv_sha256": replay_sha256,
            "python_materializer_sha256": sha256(Path(__file__)),
            "rscript_binary_sha256": sha256(args.rscript.resolve()),
            "runtime_validator_sha256": sha256(validator),
            "runtime_validation_output_sha256": runtime_sha256,
            "reference_validation_stdout_sha256": hashlib.sha256(reference_check.stdout.encode("utf-8")).hexdigest(),
            "rscript_version": "4.3.3",
            "lava_version": "0.1.5",
            "n_loci": len(loci),
            "n_pair_locus_slots": len(detailed),
            "n_univariate_tests": len(h2_rows),
            "n_bivariate_tested_replayed": len(replayed),
            "bivariate_replay_match": True,
            "standard_error_availability": "LAVA 0.1.5 run.univ/run.bivar outputs do not expose SE; no SE values are inferred",
            "outputs": {p.name: sha256(p) for p in (staged_h2_path, staged_rg_path)},
        }
        provenance_name = ("lava_family_roundoff_v1_detail_provenance.json"
                           if variant["roundoff_patch"] else "lava_family_v2_detail_provenance.json")
        prov_path = staging_dir / provenance_name
        payload = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode()
        prov_path.write_bytes(payload)
        if output_dir.exists():
            existing = {p.name: sha256(p) for p in output_dir.iterdir() if p.is_file()}
            staged = {p.name: sha256(p) for p in staging_dir.iterdir() if p.is_file()}
            if existing != staged:
                raise FileExistsError(f"Refusing to replace a different LAVA detail release: {output_dir}")
            import shutil
            shutil.rmtree(staging_dir)
        else:
            os.replace(staging_dir, output_dir)
        published = True
    finally:
        if not published and staging_dir.exists():
            import shutil
            shutil.rmtree(staging_dir, ignore_errors=True)
    print(json.dumps({"status": "PASS_DETAIL_MATERIALIZATION", "run_id": run_root.name,
                      "n_pair_locus_slots": len(detailed), "n_univariate_rows": len(h2_rows),
                      "n_bivariate_replayed": len(replayed), "output_dir": str(output_dir)}, indent=2))


if __name__ == "__main__":
    main()
