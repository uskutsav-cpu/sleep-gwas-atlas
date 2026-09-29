#!/usr/bin/env python3
"""Validate and combine the four Brain6-only continuous-duration LAVA batches.

This aggregator is intentionally separate from canonical LAVA v3. It refuses
partial, duplicate, failed, or configuration-mismatched batches and never
overwrites an existing aggregate output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REL_DIR = Path("brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1")
COMMON_FIELDS = (
    "analysis_id", "trait_id", "strict_gate_p", "reference_prefix",
    "input_root", "loci_file", "random_seed", "execution_policy", "scope",
)
ROW_FIELDS = (
    "trait_id", "locus_id", "status", "n_snps", "n_components",
    "h2_obs", "p", "strict_gate_pass", "reason",
)
SUMMARY_FIELDS = (
    "analysis_id", "scope", "rows", "tested", "not_run", "failed",
    "strict_gate_pass", "threshold",
)
NOT_RUN_REASONS = {
    "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS",
    "LOW_LOCAL_H2_UNDERPOWERED",
    "LAVA_PROCESS_LOCUS_RETURNED_NULL",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"Expected a JSON object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ROW_FIELDS:
            raise ValueError(f"Unexpected columns in {path}")
        return list(reader)


def _finite_float(value: str, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"Non-finite {label}: {value!r}")
    return number


def _is_true(value: str) -> bool:
    if value in ("TRUE", "True", "true", "1"):
        return True
    if value in ("FALSE", "False", "false", "0"):
        return False
    raise ValueError(f"Invalid Boolean value: {value!r}")


def _write_new(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())
        # Hard-link makes the create-only condition atomic on the same volume.
        os.link(temp_name, path)
    finally:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass


def aggregate(root: Path, *, write: bool = False) -> dict[str, Any]:
    root = root.resolve()
    run_dir = root / REL_DIR
    base_path = run_dir / "config.json"
    base = read_json(base_path)
    if base.get("analysis_id") != "brain6_power_optimized_sensitivity_v1_longsleep_trait_only_screen":
        raise ValueError("Unexpected analysis ID")
    if not str(base.get("scope", "")).startswith("Brain6-only"):
        raise ValueError("Screen is not explicitly Brain6-only")
    if int(base.get("execution_policy", {}).get("worker_count", -1)) != 4:
        raise ValueError("Aggregate config does not freeze four workers")

    audit_provenance_path = root / "brain6/results/power_optimized_sensitivity_v1/provenance.json"
    audit_provenance = read_json(audit_provenance_path)
    audit_root = root / "brain6/results/power_optimized_sensitivity_v1"
    for rel, expected_hash in audit_provenance.get("inputs", {}).items():
        path = root / rel
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"Inherited audit input hash mismatch: {rel}")
    for rel, expected_hash in audit_provenance.get("outputs", {}).items():
        path = audit_root / rel
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"Inherited audit output hash mismatch: {rel}")
    decision_rule = audit_provenance.get("decision_rule", {})
    minimum_improvement = int(decision_rule.get("minimum_additional_loci", -1))
    improvement_points = float(
        decision_rule.get("material_improvement_strict_gate_eligible_loci_percentage_points", -1)
    )
    if minimum_improvement != 250 or improvement_points != 10:
        raise ValueError("Predeclared promotion rule differs from the frozen Brain6 audit")
    if decision_rule.get("downstream_association_results_consulted") is not False:
        raise ValueError("Promotion decision was not outcome-blinded")

    base_ids = [str(x) for x in base.get("locus_ids", [])]
    expected_total = len(base_ids)
    if expected_total != 2495 or len(set(base_ids)) != expected_total:
        raise ValueError("Aggregate config must contain 2,495 unique frozen loci")
    threshold = float(base["strict_gate_p"])
    workers: list[dict[str, Any]] = []
    worker_rows: list[dict[str, str]] = []
    worker_input_hashes: dict[str, str] = {}
    output_paths: list[str] = []

    for index in range(1, 5):
        config_path = run_dir / f"worker_{index}.config.json"
        cfg = read_json(config_path)
        worker_input_hashes[str(config_path.relative_to(root))] = sha256(config_path)
        for field in COMMON_FIELDS:
            if cfg.get(field) != base.get(field):
                raise ValueError(f"Worker {index} configuration mismatch: {field}")
        ids = [str(x) for x in cfg.get("locus_ids", [])]
        if len(ids) != int(cfg.get("expected_loci", -1)) or len(set(ids)) != len(ids):
            raise ValueError(f"Worker {index} has an invalid locus assignment")
        result_path = (root / cfg["output_tsv"]).resolve()
        summary_path = (root / cfg["summary_json"]).resolve()
        if result_path.parent != run_dir or summary_path.parent != run_dir:
            raise ValueError(f"Worker {index} output escapes the screen directory")
        output_paths.extend((str(result_path), str(summary_path)))
        rows = read_tsv(result_path)
        summary = read_json(summary_path)
        if len(rows) != len(ids):
            raise ValueError(f"Worker {index} row count does not match its assignment")
        if [r["locus_id"] for r in rows] != ids:
            raise ValueError(f"Worker {index} row order/identity differs from its assignment")
        if any(r["trait_id"] != base["trait_id"] for r in rows):
            raise ValueError(f"Worker {index} contains a different trait")
        counts = {status: sum(r["status"] == status for r in rows)
                  for status in ("TESTED", "NOT_RUN", "FAILED")}
        gate_pass = 0
        for row in rows:
            status = row["status"]
            passed = _is_true(row["strict_gate_pass"])
            if status == "TESTED":
                p = _finite_float(row["p"], "p-value")
                _finite_float(row["h2_obs"], "observed local h2")
                if not row["n_snps"].isdigit() or not row["n_components"].isdigit():
                    raise ValueError("TESTED row has invalid SNP/component counts")
                if int(row["n_snps"]) <= 0 or int(row["n_components"]) <= 0:
                    raise ValueError("TESTED row must have positive SNP/component counts")
                if not 0 <= p <= 1 or passed != (p < threshold):
                    raise ValueError("Strict-gate flag does not match the frozen p threshold")
                gate_pass += int(passed)
                if row["reason"]:
                    raise ValueError("TESTED row has an unexpected failure reason")
            elif status == "NOT_RUN":
                if passed or row["reason"] not in NOT_RUN_REASONS:
                    raise ValueError("NOT_RUN row has an invalid reason or gate pass")
                if row["p"] not in ("", "NA") or row["h2_obs"] not in ("", "NA"):
                    raise ValueError("NOT_RUN row unexpectedly contains an estimate")
                if row["n_snps"] not in ("", "NA") or row["n_components"] not in ("", "NA"):
                    raise ValueError("NOT_RUN row unexpectedly contains SNP/component counts")
            else:
                raise ValueError(f"Unacceptable worker status {status!r}; do not aggregate")
        if counts["FAILED"]:
            raise ValueError(f"Worker {index} contains failed rows; do not aggregate")
        expected_summary = {
            "analysis_id": base["analysis_id"], "scope": base["scope"],
            "rows": len(rows), "tested": counts["TESTED"],
            "not_run": counts["NOT_RUN"], "failed": 0,
            "strict_gate_pass": gate_pass,
        }
        for field, expected in expected_summary.items():
            if summary.get(field) != expected:
                raise ValueError(f"Worker {index} summary mismatch: {field}")
        summary_threshold = float(summary.get("threshold", math.nan))
        # The R worker serializes numeric JSON values to five significant
        # digits. Accept that exact serialization as well as full precision,
        # while still rejecting a threshold from any other config.
        serialized_threshold = float(f"{threshold:.5g}")
        if summary_threshold not in (threshold, serialized_threshold):
            raise ValueError(f"Worker {index} threshold mismatch")
        worker_input_hashes[str(result_path.relative_to(root))] = sha256(result_path)
        worker_input_hashes[str(summary_path.relative_to(root))] = sha256(summary_path)
        workers.append({"worker": index, "rows": len(rows), **counts,
                        "strict_gate_pass": gate_pass})
        worker_rows.extend(rows)

    if len(set(output_paths)) != 8:
        raise ValueError("Worker output paths are not unique")
    worker_ids = [row["locus_id"] for row in worker_rows]
    if len(worker_ids) != expected_total or len(set(worker_ids)) != expected_total:
        raise ValueError("Combined worker rows are incomplete or duplicated")
    if set(worker_ids) != set(base_ids):
        raise ValueError("Combined worker rows do not match the frozen locus set")
    locus_order = {locus_id: index for index, locus_id in enumerate(base_ids)}
    worker_rows.sort(key=lambda row: locus_order[row["locus_id"]])
    if [row["locus_id"] for row in worker_rows] != base_ids:
        raise ValueError("Combined output could not be restored to frozen locus order")

    baseline_path = root / "brain6/results/power_optimized_sensitivity_v1/current_sleep_trait_audit.tsv"
    worker_input_hashes[str(audit_provenance_path.relative_to(root))] = sha256(audit_provenance_path)
    worker_input_hashes[str(baseline_path.relative_to(root))] = sha256(baseline_path)
    with baseline_path.open(newline="", encoding="utf-8") as f:
        baseline_rows = {r["trait"]: r for r in csv.DictReader(f, delimiter="\t")}
    baseline = baseline_rows["longsleep"]
    baseline_gate_pass = int(baseline["strict_univariate_gate_pass"])
    total_tested = sum(r["status"] == "TESTED" for r in worker_rows)
    total_not_run = sum(r["status"] == "NOT_RUN" for r in worker_rows)
    total_gate_pass = sum(_is_true(r["strict_gate_pass"]) for r in worker_rows)
    additional_gate_pass = total_gate_pass - baseline_gate_pass
    eligible_for_full_sensitivity = additional_gate_pass >= minimum_improvement
    summary_out = {
        "analysis_id": base["analysis_id"], "scope": base["scope"],
        "status": "COMPLETE_SCREEN_NO_FAILED_ROWS", "rows": expected_total,
        "tested": total_tested, "not_run": total_not_run, "failed": 0,
        "strict_gate_p": threshold, "strict_gate_pass": total_gate_pass,
        "canonical_longsleep_strict_gate_pass": baseline_gate_pass,
        "additional_strict_gate_pass": additional_gate_pass,
        "minimum_additional_for_promotion": minimum_improvement,
        "eligible_for_full_sensitivity_screen_gate": eligible_for_full_sensitivity,
        "canonical_decision_unchanged": "FAILED_QC_NOT_PROMOTED",
        "worker_summaries": workers,
        "interpretation": "A passing screen gate only qualifies a separate full sensitivity run; it does not promote results or change canonical QC.",
    }
    output_tsv = (root / base["output_tsv"]).resolve()
    output_json = (root / base["summary_json"]).resolve()
    provenance_path = run_dir / "aggregate.provenance.json"
    if output_tsv.parent != run_dir or output_json.parent != run_dir:
        raise ValueError("Aggregate outputs must remain inside the candidate screen directory")
    if len({output_tsv, output_json, provenance_path}) != 3:
        raise ValueError("Aggregate output paths collide")
    if any(path.exists() for path in (output_tsv, output_json, provenance_path)):
        raise FileExistsError("Refusing to overwrite an existing aggregate output")
    tsv_content = "\t".join(ROW_FIELDS) + "\n" + "".join(
        "\t".join(row[field] for field in ROW_FIELDS) + "\n" for row in worker_rows
    )
    json_content = json.dumps(summary_out, indent=2, sort_keys=True) + "\n"
    provenance = {
        "schema_version": 1, "analysis_id": base["analysis_id"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Brain6-only; separate outcome-blinded trait-only sensitivity screen",
        "worker_count": 4, "base_config_sha256": sha256(base_path),
        "aggregator_sha256": sha256(Path(__file__).resolve()),
        "upstream_audit_provenance_sha256": sha256(audit_provenance_path),
        "upstream_audit_inputs": audit_provenance.get("inputs", {}),
        "upstream_audit_outputs": audit_provenance.get("outputs", {}),
        "inputs": worker_input_hashes,
        "outputs": {
            str(output_tsv.relative_to(root)): hashlib.sha256(tsv_content.encode()).hexdigest(),
            str(output_json.relative_to(root)): hashlib.sha256(json_content.encode()).hexdigest(),
        },
        "status": summary_out["status"],
    }
    provenance_content = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    if write:
        # Recheck before writing so concurrent reruns fail closed.
        if any(path.exists() for path in (output_tsv, output_json, provenance_path)):
            raise FileExistsError("Aggregate output appeared during validation")
        _write_new(output_tsv, tsv_content)
        # Never remove prior outputs; a partial aggregate remains visible for
        # explicit recovery rather than being silently rolled back.
        _write_new(output_json, json_content)
        _write_new(provenance_path, provenance_content)
    summary_out["write_requested"] = write
    return summary_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--write", action="store_true",
                        help="write configured aggregate files; refuses any overwrite")
    args = parser.parse_args()
    print(json.dumps(aggregate(args.root, write=args.write), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
