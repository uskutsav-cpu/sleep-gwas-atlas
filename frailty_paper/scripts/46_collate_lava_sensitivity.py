#!/usr/bin/env python3
"""Validate and collate the complete locked FI x sleep LAVA sensitivity family."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCK = ROOT / "frailty_paper/config/lava_frailty_sensitivity_v1.yaml"
OUTPUT_ROOT = Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1")
TRAITS = ("accel_sleep_duration", "chronotype", "insomnia", "longsleep", "napping", "shortsleep",
          "sleep_apnea", "sleep_efficiency", "sleep_timing", "sleepdur", "sleepiness", "snoring")
LOCI = 2495
THRESHOLD = 4.45335114673792e-7
BH_ALPHA = 0.05
MAX_UNIVARIATE_UNTESTED = 0.05
MAX_LOCUS_FAILURE = 0.01
MAX_BIVARIATE_FAILURE = 0.05


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh_adjust(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original] * n / rank)
        adjusted[original] = min(1.0, running)
    return adjusted


def atomic_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    temp = path.with_name(path.name + f".{os.getpid()}.partial")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def normalize_p(value) -> float:
    return float(value) if isinstance(value, (int, float)) and math.isfinite(value) and 0 <= value <= 1 else 1.0


def main() -> int:
    lock_text = LOCK.read_text(encoding="utf-8")
    analysis_match = re.search(r"^analysis_id:\s*([^\s#]+)\s*$", lock_text, re.MULTILINE)
    runtime_block = re.search(r"^runtime:\s*\n((?:^[ \t]+.*\n?)*)", lock_text, re.MULTILINE)
    reference_match = re.search(r"^[ \t]+reference_root:\s*(.+?)\s*$", runtime_block.group(1), re.MULTILINE) if runtime_block else None
    if not analysis_match or not reference_match:
        raise RuntimeError("analysis ID or pinned reference root is absent from frozen lock")
    reference_root = Path(reference_match.group(1))
    reference_manifest = reference_root / "extracted_manifest.tsv"
    if not reference_manifest.is_file():
        raise RuntimeError("pinned external LAVA reference manifest is missing")
    lock_hash = sha256(LOCK)
    input_path = OUTPUT_ROOT / "input_manifest.json"
    if not input_path.is_file():
        raise RuntimeError("prepared input manifest is absent")
    input_hash = sha256(input_path)
    inputs = json.loads(input_path.read_text(encoding="utf-8"))
    if inputs.get("analysis_id") != analysis_match.group(1) or inputs.get("lock_sha256") != lock_hash:
        raise RuntimeError("prepared inputs do not match the current lock")
    if inputs.get("candidate_pair_locus_slots") != 29940 or set(inputs.get("pairs", {})) != set(TRAITS):
        raise RuntimeError("prepared family dimensions or trait set drifted")
    loci_path = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
    with loci_path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        loci = [line.split() for line in handle if line.strip()]
    if header != ["LOC", "CHR", "START", "STOP"] or len(loci) != LOCI:
        raise RuntimeError("frozen locus file is invalid")
    runner_script = ROOT / "frailty_paper/scripts/45_run_lava_sensitivity.py"
    worker_script = ROOT / "frailty_paper/scripts/44_run_lava_sensitivity_locus.R"
    prepare_script = ROOT / "frailty_paper/scripts/43_prepare_lava_sensitivity_inputs.py"
    locus_by_index = {i: row for i, row in enumerate(loci, 1)}

    pair_rows: list[dict] = []
    bivar_rows: list[dict] = []
    per_locus_fi: dict[int, list[dict]] = defaultdict(list)
    fi_context_status: dict[int, list[str]] = defaultdict(list)
    univ_rows: list[dict] = []
    missing: list[str] = []
    for pair_order, trait in enumerate(TRAITS, 1):
        pair_root = OUTPUT_ROOT / "pairs" / trait
        if set(inputs["pairs"].get(trait, {}).keys()) < {"overlap_file", "input_info_files"}:
            raise RuntimeError(f"pair input manifest incomplete: {trait}")
        for index in range(1, LOCI + 1):
            path = pair_root / "loci" / f"locus_{index:04d}.json"
            if not path.is_file():
                missing.append(str(path))
                continue
            rec = json.loads(path.read_text(encoding="utf-8"))
            expected = locus_by_index[index]
            identity = (rec.get("schema_version"), rec.get("trait"), rec.get("pair_order"), rec.get("locus_index"))
            if identity != (1, trait, pair_order, index):
                raise RuntimeError(f"receipt identity mismatch: {path}")
            if [str(rec.get(k)) for k in ("LOC", "CHR", "START", "STOP")] != expected:
                raise RuntimeError(f"locus coordinates mismatch: {path}")
            univ_status = rec.get("univ_status")
            process_status = rec.get("process_status")
            if process_status not in {"PROCESSED", "PROCESS_FAILED"}:
                raise RuntimeError(f"unknown process status in {path}")
            if univ_status not in {"TESTED", "PARTIAL", "UNIVARIATE_FAILED", "UNIVARIATE_INVALID", "PHENOTYPE_DROPPED", "NOT_RUN"}:
                raise RuntimeError(f"unknown univariate status in {path}")
            fi_status = rec.get("fi_status")
            sleep_status = rec.get("sleep_status")
            per_trait_statuses = {"TESTED", "PROCESS_FAILED", "UNIVARIATE_FAILED", "UNIVARIATE_INVALID", "PHENOTYPE_DROPPED", "NOT_RUN"}
            if fi_status not in per_trait_statuses or sleep_status not in per_trait_statuses:
                raise RuntimeError(f"unknown per-trait univariate status in {path}")
            if rec.get("bivar_status") not in {"TESTED", "NOT_ELIGIBLE", "BIVARIATE_FAILED", "BIVARIATE_INVALID"}:
                raise RuntimeError(f"unknown bivariate status in {path}")
            ps = []
            for status_key, h2_key, p_key in (("fi_status", "fi_h2_obs", "fi_p"), ("sleep_status", "sleep_h2_obs", "sleep_p")):
                if rec[status_key] == "TESTED":
                    if not isinstance(rec.get(h2_key), (int, float)) or not math.isfinite(rec[h2_key]) or rec[h2_key] < 0:
                        raise RuntimeError(f"invalid tested univariate h2 in {path}")
                    if not isinstance(rec.get(p_key), (int, float)) or not math.isfinite(rec[p_key]) or not 0 <= rec[p_key] <= 1:
                        raise RuntimeError(f"invalid tested univariate p-value in {path}")
                    ps.append(float(rec[p_key]))
                else:
                    ps.append(1.0)
            eligible = fi_status == "TESTED" and sleep_status == "TESTED" and ps[0] <= THRESHOLD and ps[1] <= THRESHOLD
            b_status = rec["bivar_status"]
            if eligible and b_status not in {"TESTED", "BIVARIATE_FAILED", "BIVARIATE_INVALID"}:
                raise RuntimeError(f"eligible bivariate test missing explicit result or failure in {path}")
            if not eligible and b_status == "TESTED":
                raise RuntimeError(f"ineligible slot contains a tested bivariate result in {path}")
            pair_rows.append({
                "trait": trait, "pair_order": pair_order, "locus_index": index,
                "LOC": rec["LOC"], "CHR": rec["CHR"], "START": rec["START"], "STOP": rec["STOP"],
                "process_status": process_status, "univariate_status": univ_status,
                "fi_status": fi_status, "sleep_status": sleep_status,
                "fi_local_p": rec.get("fi_p"), "sleep_local_p": rec.get("sleep_p"),
                "bivariate_eligible": str(eligible).lower(), "bivariate_status": b_status,
                "n_snps": rec.get("n_snps"), "K": rec.get("K"),
                "elapsed_seconds": rec.get("elapsed_seconds"), "error": rec.get("error")})
            fi_context_status[index].append(fi_status)
            if fi_status == "TESTED":
                per_locus_fi[index].append({"trait": trait, "h2": rec["fi_h2_obs"], "p": ps[0], "status": univ_status})
            if sleep_status == "TESTED":
                univ_rows.append({"phenotype": trait, "locus_index": index, "LOC": rec["LOC"],
                                  "CHR": rec["CHR"], "START": rec["START"], "STOP": rec["STOP"],
                                  "h2_obs": rec["sleep_h2_obs"], "p": ps[1], "status": "TESTED",
                                  "source_pair": trait})
            else:
                univ_rows.append({"phenotype": trait, "locus_index": index, "LOC": rec["LOC"],
                                  "CHR": rec["CHR"], "START": rec["START"], "STOP": rec["STOP"],
                                  "h2_obs": None, "p": None, "status": sleep_status, "source_pair": trait})
            bivar_rows.append({
                "trait": trait, "locus_index": index, "LOC": rec["LOC"], "CHR": rec["CHR"],
                "START": rec["START"], "STOP": rec["STOP"], "eligible": str(eligible).lower(),
                "status": b_status if eligible else "INELIGIBLE",
                "rho": rec.get("rho") if eligible and b_status == "TESTED" else None,
                "rho_lower": rec.get("rho_lower") if eligible and b_status == "TESTED" else None,
                "rho_upper": rec.get("rho_upper") if eligible and b_status == "TESTED" else None,
                "r2": rec.get("r2") if eligible and b_status == "TESTED" else None,
                "r2_lower": rec.get("r2_lower") if eligible and b_status == "TESTED" else None,
                "r2_upper": rec.get("r2_upper") if eligible and b_status == "TESTED" else None,
                "p": normalize_p(rec.get("bivar_p")) if eligible and b_status == "TESTED" else 1.0,
                "q": None})

    if missing:
        raise RuntimeError(f"incomplete LAVA family: {len(missing)} receipts missing; first={missing[0]}")
    if len(pair_rows) != 29940 or len(bivar_rows) != 29940:
        raise RuntimeError("receipt counts do not equal the locked 29,940 slots")
    fi_tested = 0
    fi_conflicts = []
    for index in range(1, LOCI + 1):
        values = per_locus_fi.get(index, [])
        if values:
            reference = values[0]
            if any(row["p"] != reference["p"] or row["h2"] != reference["h2"] for row in values[1:]):
                fi_conflicts.append(index)
            elif len(values) == 12 and len(fi_context_status[index]) == 12:
                fi_tested += 1
            if len(values) == 12 and len(fi_context_status[index]) == 12:
                fi_status = "TESTED"
                fi_h2, fi_p = reference["h2"], reference["p"]
            else:
                fi_status = "INCOMPLETE_PAIR_CONTEXTS"
                fi_h2 = fi_p = None
        else:
            fi_status = "UNTESTED"
            fi_h2 = fi_p = None
        loc = loci[index - 1]
        univ_rows.append({"phenotype": "frailty", "locus_index": index, "LOC": loc[0],
                          "CHR": loc[1], "START": loc[2], "STOP": loc[3],
                          "h2_obs": fi_h2, "p": fi_p, "status": fi_status,
                          "source_pair": "all_12_pair_contexts"})
    if fi_conflicts:
        raise RuntimeError(f"FI local-univariate estimates conflict across pair contexts at {len(fi_conflicts)} loci; first={fi_conflicts[0]}")
    univ_rows = [r for r in univ_rows if r["phenotype"] != "frailty"] + [r for r in univ_rows if r["phenotype"] == "frailty"]
    if len(univ_rows) != 32435:
        raise RuntimeError(f"unique local-univariate table has {len(univ_rows)} rows, expected 32,435")
    for row in univ_rows:
        p = normalize_p(row.get("p"))
        row["p"] = p if row["status"] == "TESTED" else None
        row["eligible_at_locked_threshold"] = str(row["status"] == "TESTED" and p <= THRESHOLD).lower()

    qvalues = bh_adjust([float(row["p"]) for row in bivar_rows])
    for row, q in zip(bivar_rows, qvalues):
        row["q"] = q
        row["q_le_0_05"] = str(q <= BH_ALPHA).lower()
    eligible = [r for r in pair_rows if r["bivariate_eligible"] == "true"]
    bivar_fail = sum(r["bivariate_status"] in {"BIVARIATE_FAILED", "BIVARIATE_INVALID"} for r in eligible)
    locus_fail_indices = {r["locus_index"] for r in pair_rows
                          if r["process_status"] == "PROCESS_FAILED" or r["fi_status"] != "TESTED" or r["sleep_status"] != "TESTED"}
    locus_fail = len(locus_fail_indices)
    univ_untested = sum(r["status"] != "TESTED" for r in univ_rows)
    family = [
        {"gate": "candidate_pair_locus_slots", "observed": len(pair_rows), "expected": 29940, "limit": "exact", "pass": str(len(pair_rows) == 29940).lower()},
        {"gate": "local_univariate_untested_fraction", "observed": univ_untested / 32435, "expected": 32435, "limit": MAX_UNIVARIATE_UNTESTED, "pass": str(univ_untested / 32435 <= MAX_UNIVARIATE_UNTESTED).lower()},
        {"gate": "locus_failure_fraction", "observed": locus_fail / LOCI, "expected": LOCI, "limit": MAX_LOCUS_FAILURE, "pass": str(locus_fail / LOCI <= MAX_LOCUS_FAILURE).lower()},
        {"gate": "eligible_bivariate_failure_fraction", "observed": bivar_fail / max(1, len(eligible)), "expected": len(eligible), "limit": MAX_BIVARIATE_FAILURE, "pass": str(bivar_fail / max(1, len(eligible)) <= MAX_BIVARIATE_FAILURE).lower()},
    ]
    family_complete = all(r["pass"] == "true" for r in family)
    qc = {"analysis_id": analysis_match.group(1), "lock_sha256": lock_hash,
          "input_manifest_sha256": input_hash, "candidate_pair_locus_slots": 29940,
          "source_file_checksums": inputs["traits"],
          "locus_file_sha256": sha256(loci_path),
          "reference_root": str(reference_root), "reference_manifest_sha256": sha256(reference_manifest),
          "script_sha256": {"prepare": sha256(prepare_script), "worker": sha256(worker_script),
                            "runner": sha256(runner_script), "collator": sha256(Path(__file__))},
          "locked_runtime": {"R": "4.3.3", "LAVA": "0.1.5", "LAVA_commit": "e729a245f7b6923967a96804fbf5246eadf2d6c6"},
          "command_template": ["python3", "frailty_paper/scripts/45_run_lava_sensitivity.py",
                               "<trait> <pair_order> <locus_index> <input_info> <overlap_file> <locus_file> <reference_root> <input_alias> <receipt_json>"],
          "unique_local_univariate_tests": 32435, "local_univariate_tested": sum(r["status"] == "TESTED" for r in univ_rows),
          "local_univariate_untested": univ_untested, "fi_loci_tested_in_all_12_contexts": fi_tested,
          "eligible_bivariate_slots": len(eligible), "bivariate_tested": len(eligible) - bivar_fail,
          "bivariate_failed_eligible": bivar_fail, "pair_locus_failures": locus_fail,
          "family_gates_pass": family_complete, "sensitivity_only": True,
          "generated_utc": datetime.now(timezone.utc).isoformat(), "gates": family}

    univ_cols = ["phenotype", "locus_index", "LOC", "CHR", "START", "STOP", "h2_obs", "p", "status", "eligible_at_locked_threshold", "source_pair"]
    pair_cols = ["trait", "pair_order", "locus_index", "LOC", "CHR", "START", "STOP", "process_status", "univariate_status", "fi_status", "sleep_status", "fi_local_p", "sleep_local_p", "bivariate_eligible", "bivariate_status", "n_snps", "K", "elapsed_seconds", "error"]
    biv_cols = ["trait", "locus_index", "LOC", "CHR", "START", "STOP", "eligible", "status", "rho", "rho_lower", "rho_upper", "r2", "r2_lower", "r2_upper", "p", "q", "q_le_0_05"]
    atomic_tsv(OUTPUT_ROOT / "local_univariate.tsv", univ_rows, univ_cols)
    atomic_tsv(OUTPUT_ROOT / "local_bivariate.tsv", bivar_rows, biv_cols)
    atomic_tsv(OUTPUT_ROOT / "pair_locus_status.tsv", pair_rows, pair_cols)
    atomic_tsv(OUTPUT_ROOT / "family_qc.tsv", family, ["gate", "observed", "expected", "limit", "pass"])
    atomic_json(OUTPUT_ROOT / "run_manifest.json", qc)
    print(f"COLLATION_COMPLETE univariate={len(univ_rows)} pair_locus={len(pair_rows)} eligible_bivariate={len(eligible)} gates_pass={family_complete}")
    return 0 if family_complete else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
