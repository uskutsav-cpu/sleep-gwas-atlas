#!/usr/bin/env python3
"""Lock all 45-trait by all-context S-PrediXcan tests before any TWAS result."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from datetime import datetime, timezone
from pathlib import Path


MAPPING_FIELDS = [
    "mapping_id", "trait_id", "model_family", "modality", "full_gwas_path", "full_gwas_sha256",
    "model_variant_path", "model_variant_sha256", "mapped_gwas_path", "mapped_gwas_lock_path",
    "results_accessed_before_lock",
]
ELIGIBILITY_FIELDS = [
    "trait_id", "trait_domain", "trait_type", "gwas_N", "gwas_h2", "gwas_h2_scale",
    "analysis_status", "reason", "results_accessed_before_lock",
]
RUN_FIELDS = [
    "run_id", "trait_id", "trait_domain", "model_id", "model_family", "modality", "context",
    "mapping_id", "mapped_gwas_path", "mapped_gwas_lock_path", "model_db_path", "model_db_sha256",
    "covariance_path", "covariance_sha256", "model_snp_key", "gwas_N", "gwas_h2",
    "gwas_h2_scale", "variance_control_status", "output_path",
    "results_accessed_before_lock",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def qc_prefiltered(path: Path) -> bool:
    return any(line.startswith("prefilter_strategy\t") for line in path.read_text(encoding="utf-8").splitlines())


def full_input(root: Path, trait: str) -> Path | None:
    for directory in ("data/harmonized_mixer_full", "data/harmonized"):
        data = root / directory / f"{trait}.harmonized.tsv.gz"
        qc = root / directory / f"{trait}.qc.txt"
        if data.is_file() and qc.is_file() and not qc_prefiltered(qc):
            return data
    return None


def trait_eligibility(trait: dict[str, str], evidence: dict[str, str]) -> dict[str, object]:
    reason = "PASS"
    try:
        n_total = int(evidence["n_total"])
        h2 = float(evidence["h2"])
    except (KeyError, ValueError):
        n_total, h2, reason = 0, float("nan"), "MISSING_OR_NONNUMERIC_N_OR_H2"
    if reason == "PASS" and n_total <= 0:
        reason = "NONPOSITIVE_N_TOTAL"
    if reason == "PASS" and (not math.isfinite(h2) or not 0 < h2 <= 1):
        reason = "H2_OUTSIDE_OPEN_CLOSED_0_1"
    status = "ELIGIBLE" if reason == "PASS" else "NOT_APPLICABLE"
    return {
        "trait_id": trait["trait_id"], "trait_domain": trait["domain"], "trait_type": trait["type"],
        "gwas_N": n_total if n_total > 0 else "NA", "gwas_h2": format(h2, ".15g") if math.isfinite(h2) else "NA",
        "gwas_h2_scale": evidence.get("h2_scale", "") or "NA", "analysis_status": status,
        "reason": reason, "results_accessed_before_lock": "NO",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--preflight", default="results/tables/molecular_preflight.json")
    parser.add_argument("--models", default="results/tables/twas_model_registry.tsv")
    parser.add_argument("--models-lock", default="results/tables/twas_model_registry.lock.json")
    parser.add_argument("--mapping-out", default="results/tables/twas_gwas_mapping_manifest.tsv")
    parser.add_argument("--eligibility-out", default="results/tables/twas_trait_eligibility.tsv")
    parser.add_argument("--run-out", default="results/tables/twas_run_manifest.tsv")
    parser.add_argument("--lock-out", default="results/tables/twas_run_manifest.lock.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, preflight_path = root / args.policy, root / args.preflight
    model_path, model_lock_path = root / args.models, root / args.models_lock
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    if preflight.get("policy_sha256") != sha256(policy_path) or not preflight.get("metaxcan_runtime_ready"):
        fail("pinned molecular preflight does not have a ready MetaXcan runtime")
    if not preflight.get("upstream_fine_mapping_ready"):
        fail("TWAS manifest remains downstream of locked fine-mapping outputs")
    model_lock = json.loads(model_lock_path.read_text(encoding="utf-8"))
    model_fields, models = read_tsv(model_path)
    del model_fields
    if model_lock.get("registry_sha256") != sha256(model_path) or model_lock.get("model_ids_in_locked_order") != [row["model_id"] for row in models]:
        fail("TWAS model registry differs from lock")
    panel_fields, panel = read_tsv(root / "config/analysis_panel.tsv")
    del panel_fields
    if len(panel) != 45 or len({row["trait_id"] for row in panel}) != 45:
        fail("TWAS trait family is not the locked 45-trait panel")
    full = {row["trait_id"]: full_input(root, row["trait_id"]) for row in panel}
    missing = [trait for trait, path in full.items() if path is None]
    if missing:
        fail(f"full dense TWAS inputs are absent for {len(missing)} traits: {','.join(missing)}")
    atlas_traits_path = root / "results/atlas/traits.tsv"
    _, atlas_traits = read_tsv(atlas_traits_path)
    atlas_by_id = {row["trait_id"]: row for row in atlas_traits}
    if len(atlas_traits) != 45 or set(atlas_by_id) != {row["trait_id"] for row in panel}:
        fail("atlas trait evidence does not cover the exact locked panel")
    eligibility_rows: list[dict[str, object]] = []
    eligible: dict[str, dict[str, object]] = {}
    for trait in panel:
        evidence = atlas_by_id[trait["trait_id"]]
        assessment = trait_eligibility(trait, evidence)
        eligibility_rows.append(assessment)
        if assessment["analysis_status"] == "ELIGIBLE":
            eligible[trait["trait_id"]] = {
                "gwas_N": assessment["gwas_N"], "gwas_h2": assessment["gwas_h2"],
                "gwas_h2_scale": assessment["gwas_h2_scale"],
            }
    mapping_rows: list[dict[str, object]] = []
    modality_by_family: dict[str, str] = {}
    for model in models:
        family, modality = model["model_family"], model["modality"]
        if family in modality_by_family and modality_by_family[family] != modality:
            fail(f"TWAS model family has multiple modalities: {family}")
        modality_by_family[family] = modality
    if set(modality_by_family) != set(policy["twas"]["model_families"]):
        fail("TWAS model registry does not cover the locked model families")
    for trait in panel:
        if trait["trait_id"] not in eligible:
            continue
        for family in policy["twas"]["model_families"]:
            modality = modality_by_family[family]
            variant_path = root / "results/tables/twas_model_variants" / f"{family}.tsv.gz"
            if not variant_path.is_file():
                fail(f"TWAS model variant registry is absent: {family}")
            mapping_id = f"TWMAP__{trait['trait_id']}__{family}"
            mapped = root / "data/twas" / family / f"{trait['trait_id']}.tsv.gz"
            mapping_rows.append({
                "mapping_id": mapping_id, "trait_id": trait["trait_id"], "model_family": family,
                "modality": modality, "full_gwas_path": str(full[trait["trait_id"]].relative_to(root)),
                "full_gwas_sha256": sha256(full[trait["trait_id"]]),
                "model_variant_path": str(variant_path.relative_to(root)), "model_variant_sha256": sha256(variant_path),
                "mapped_gwas_path": str(mapped.relative_to(root)),
                "mapped_gwas_lock_path": str((mapped.with_suffix(".lock.json")).relative_to(root)),
                "results_accessed_before_lock": "NO",
            })
    mapping_by_key = {(row["trait_id"], row["model_family"]): row for row in mapping_rows}
    run_rows: list[dict[str, object]] = []
    for trait in panel:
        if trait["trait_id"] not in eligible:
            continue
        for model in models:
            mapping = mapping_by_key[(trait["trait_id"], model["model_family"])]
            calibration = eligible[trait["trait_id"]]
            run_id = f"TWAS__{trait['trait_id']}__{model['model_id']}"
            run_rows.append({
                "run_id": run_id, "trait_id": trait["trait_id"], "trait_domain": trait["domain"],
                "model_id": model["model_id"], "model_family": model["model_family"],
                "modality": model["modality"], "context": model["context"],
                "mapping_id": mapping["mapping_id"], "mapped_gwas_path": mapping["mapped_gwas_path"],
                "mapped_gwas_lock_path": mapping["mapped_gwas_lock_path"],
                "model_db_path": model["model_db_path"], "model_db_sha256": model["model_db_sha256"],
                "covariance_path": model["covariance_path"], "covariance_sha256": model["covariance_sha256"],
                "model_snp_key": model["model_snp_key"],
                "gwas_N": calibration["gwas_N"], "gwas_h2": format(float(calibration["gwas_h2"]), ".15g"),
                "gwas_h2_scale": calibration["gwas_h2_scale"], "variance_control_status": "REQUIRED",
                "output_path": f"results/twas/runs/{run_id}/spredixcan.csv", "results_accessed_before_lock": "NO",
            })
    run_ids = [str(row["run_id"]) for row in run_rows]
    if len(run_ids) != len(set(run_ids)) or len(run_rows) != len(eligible) * len(models):
        fail("TWAS run family is incomplete or duplicated")
    mapping_out, eligibility_out = root / args.mapping_out, root / args.eligibility_out
    run_out, lock_out = root / args.run_out, root / args.lock_out
    if any(path.exists() for path in (mapping_out, eligibility_out, run_out, lock_out)):
        fail("TWAS manifests already exist; refusing overwrite")
    mapping_out.parent.mkdir(parents=True, exist_ok=True)
    mapping_out.write_text(table_text(MAPPING_FIELDS, mapping_rows), encoding="utf-8")
    eligibility_out.write_text(table_text(ELIGIBILITY_FIELDS, eligibility_rows), encoding="utf-8")
    run_out.write_text(table_text(RUN_FIELDS, run_rows), encoding="utf-8")
    lock_payload = {
        "schema_version": "atlas-v1.0-twas-runs.1", "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "results_accessed_before_lock": False, "trait_count": 45, "eligible_trait_count": len(eligible),
        "not_applicable_trait_count": 45 - len(eligible), "model_count": len(models),
        "mapping_count": len(mapping_rows), "run_count": len(run_rows), "run_ids_in_locked_order": run_ids,
        "mapping_sha256": sha256(mapping_out), "run_manifest_sha256": sha256(run_out),
        "eligibility_sha256": sha256(eligibility_out), "atlas_traits_sha256": sha256(atlas_traits_path),
        "model_registry_sha256": sha256(model_path), "model_registry_lock_sha256": sha256(model_lock_path),
        "policy_sha256": sha256(policy_path), "preflight_sha256": sha256(preflight_path),
        "claim_limit": policy["claim_limit"],
    }
    lock_out.write_text(json.dumps(lock_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"TWAS_RUN_MANIFEST_OK traits=45 eligible={len(eligible)} models={len(models)} mappings={len(mapping_rows)} runs={len(run_rows)} result_free=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
