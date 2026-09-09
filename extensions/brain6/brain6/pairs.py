"""Select from measured atlas rows, not biological plausibility or old proposals."""
from __future__ import annotations
import math
from pathlib import Path
from .io import read_json, read_tsv, write_json, write_tsv, require, file_record, json_hash, safe_id

BRAIN = ["adhd", "mdd", "scz", "bipolar", "alz", "parkinson"]
SLEEP = ["insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype", "sleepiness",
         "napping", "snoring", "sleep_apnea", "sleep_efficiency", "accel_sleep_duration", "sleep_timing"]


def rank_pairs(matrix: str | Path, config: dict) -> list[dict]:
    """Preserve existing family FDR; do not recalculate after choosing six diseases."""
    names = {k: k for k in ["sleep_trait", "disease_trait", "rg", "se", "p", "fdr"]}
    names.update(config.get("matrix_columns", {}))
    qc_col = config.get("matrix_qc_column", "analysis_status")
    required = list(names.values()) + [qc_col]
    rows = list(read_tsv(matrix, required))
    expected = config.get("expected_matrix_rows", 396)
    require(len(rows) == expected, f"Expected {expected} atlas rows; got {len(rows)}")
    valid_sleep = config.get("sleep_traits", SLEEP)
    brains = config.get("brain_traits", BRAIN)
    require(len(valid_sleep) == len(set(valid_sleep)) == 12, "Expected 12 distinct sleep traits")
    require(set(brains) == set(BRAIN), "Six chosen disorder identities changed")
    seen = set()
    selected = []
    for raw in rows:
        r = {key: raw[col] for key, col in names.items()}
        pair = (r["sleep_trait"], r["disease_trait"])
        require(pair not in seen, f"Duplicate atlas pair: {pair}")
        seen.add(pair)
        if r["disease_trait"] not in brains:
            continue
        require(r["sleep_trait"] in valid_sleep, f"Unknown sleep trait: {pair}")
        status = raw[qc_col]
        flags = []
        for k in ["rg", "se", "p", "fdr"]:
            try:
                r[k] = float(r[k])
            except (ValueError, TypeError):
                r[k] = None
        numeric = all(r[k] is not None and math.isfinite(r[k]) for k in ["rg", "se", "p", "fdr"])
        if not numeric:
            flags.append("NONFINITE_RESULT")
        else:
            if abs(r["rg"]) > 1 or r["se"] <= 0 or not 0 <= r["p"] <= 1 or not 0 <= r["fdr"] <= 1:
                flags.append("INVALID_ESTIMATE")
            if r["fdr"] >= config.get("global_fdr_max", .05):
                flags.append("NOT_GLOBAL_FDR_SIGNIFICANT")
        if status not in config.get("primary_qc_values", ["PRIMARY"]):
            flags.append("NOT_PRIMARY_QC")
        r.update(analysis_status=status, eligible=not flags, reason=";".join(flags) or "ELIGIBLE",
                 selection_rank="", pair_id=safe_id(f"{pair[0]}__{pair[1]}"))
        selected.append(r)
    require({(r["sleep_trait"], r["disease_trait"]) for r in selected} ==
            {(s, b) for s in valid_sleep for b in brains}, "Incomplete 12 x 6 brain subset")
    # Sort by existing family FDR first, then effect magnitude, then stable ID.
    # This is a transparent prioritization rule, NOT an optimized weighted score.
    output = []
    for b in brains:
        group = [r for r in selected if r["disease_trait"] == b]
        group.sort(key=lambda r: (not r["eligible"], r["fdr"] if r["fdr"] is not None else 1,
                                  -abs(r["rg"]) if r["rg"] is not None else 0, r["pair_id"]))
        rank = 0
        for r in group:
            if r["eligible"]:
                rank += 1; r["selection_rank"] = rank
        output.extend(group)
    return output


def write_ranking(matrix, config_path, out):
    rows = rank_pairs(matrix, read_json(config_path))
    write_tsv(out, list(rows[0]), rows)
    return rows


def freeze_pairs(matrix, config_path, decisions_path, out, reviewer: str, *, synthetic=False):
    """Explicit pair decision table avoids freezing by proposal example alone."""
    require(bool(reviewer.strip()) and reviewer.lower() not in {"todo", "unresolved"}, "Named scientific review required")
    cfg = read_json(config_path)
    require(cfg.get("synthetic", False) == synthetic, "Synthetic lock must be explicitly marked")
    ranked = {r["pair_id"]: r for r in rank_pairs(matrix, cfg)}
    decisions = list(read_tsv(decisions_path, ["disease_trait", "sleep_trait", "role", "reason"]))
    primary = set(); included = []; keys = set()
    for d in decisions:
        require(d["disease_trait"] in BRAIN, "Unknown disorder")
        require(d["role"] in {"PRIMARY", "SECONDARY", "NO_ELIGIBLE_PAIR"}, "Unknown pair role")
        require(d["reason"].strip(), "Each decision needs a justification")
        if d["role"] == "NO_ELIGIBLE_PAIR":
            require(not any(r["eligible"] for r in ranked.values() if r["disease_trait"] == d["disease_trait"]),
                    "Cannot call NO_ELIGIBLE_PAIR when eligible evidence exists")
            require(d["disease_trait"] not in primary, "Duplicate disorder decision")
            primary.add(d["disease_trait"])
            continue
        key = f"{d['sleep_trait']}__{d['disease_trait']}"
        require(key in ranked and ranked[key]["eligible"], f"Ineligible/unverified selected pair: {key}")
        require(key not in keys, "Duplicate selected pair")
        keys.add(key)
        if d["role"] == "PRIMARY":
            require(d["disease_trait"] not in primary, "Multiple primary partners for one disorder")
            primary.add(d["disease_trait"])
        included.append({**ranked[key], "role": d["role"], "decision_reason": d["reason"]})
    require(primary == set(BRAIN), "Each disorder needs a primary pair or explicit NO_ELIGIBLE_PAIR")
    require(0 <= len(included) <= 12, "Choose at most 12 pairs; a fully null panel remains reportable")
    for b in BRAIN:
        require(sum(r["disease_trait"] == b and r["role"] == "SECONDARY" for r in included) <= 1,
                "At most one secondary per disorder")
    native = cfg.get("native_policy", {})
    require(native.get("pleiotropy") in {"PLACO_PLUS", "PLACO"}, "Declare the new family method")
    if native["pleiotropy"] == "PLACO":
        require(native.get("uncorrelated_traits_justification"), "Original PLACO requires uncorrelated inputs justification")
    lock = {"schema_version": 1, "experiment": cfg.get("experiment", "brain6-v1"),
            "synthetic": synthetic, "reviewed_by": reviewer, "selection_is_post_global_screen": True,
            "atlas": file_record(matrix), "config": file_record(config_path),
            "decisions": file_record(decisions_path), "pairs": included,
            "n_pairs": len(included),
            "disorder_decisions": decisions,
            "genomewide_family_threshold": 5e-8 / len(included) if included else None,
            "threshold_interpretation": "conventional genomewide threshold with additional Bonferroni across selected pairs; not a selective-inference correction",
            "native_policy": native, "protected_legacy_track": "A/B/CONTROL immutable",
            "note": "Discovery-selected follow-up is not independent replication. Preserve old core FDR."}
    lock["lock_sha256"] = json_hash(lock)
    write_json(out, lock)
    return lock


def verify_lock(path):
    lock = read_json(path)
    digest = lock.pop("lock_sha256", None)
    require(digest == json_hash(lock), "Pair lock was edited/corrupted")
    lock["lock_sha256"] = digest
    require(lock["n_pairs"] == len(lock["pairs"]), "Broken family size")
    require(len({r["pair_id"] for r in lock["pairs"]}) == lock["n_pairs"], "Duplicate locked pair")
    expected = 5e-8 / lock["n_pairs"] if lock["n_pairs"] else None
    require(lock.get("genomewide_family_threshold") == expected, "Broken family threshold")
    return lock
