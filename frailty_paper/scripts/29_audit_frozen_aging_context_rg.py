#!/usr/bin/env python3
"""Read-only extraction of selected aging-context pairs from the frozen atlas."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
CONTEXT_TRAITS = (
    "healthspan", "parental_lifespan", "longevity", "grip_strength", "alz", "parkinson",
)
PAIR_FIELDS = (
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier", "interpretation_status",
    "global_rg", "global_rg_se", "global_rg_z", "global_rg_p",
    "global_rg_fdr_all_396", "global_rg_fdr_primary_372", "global_rg_primary_significant",
    "effect_direction",
)
TRAIT_FIELDS = (
    "label", "source_id", "dataset_version", "ancestry", "build", "ncase", "ncontrol",
    "n_total", "h2_scale", "h2", "h2_se", "h2_z", "ldsc_intercept", "phase1_verdict",
)
OUTPUT_FIELDS = (*PAIR_FIELDS, "context_label", "source_id", "dataset_version", "ancestry", "build",
                 "ncase", "ncontrol", "n_total", "h2_scale", "h2", "h2_se", "h2_z",
                 "ldsc_intercept", "phase1_verdict", "reuse_status", "fdr_denominator",
                 "cohort_overlap_status", "claim_limit")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_context_rows(
    sleep_traits: list[str], pair_rows: list[dict[str, str]],
    trait_rows: list[dict[str, str]], context_traits: tuple[str, ...] = CONTEXT_TRAITS,
) -> list[dict[str, str]]:
    if len(sleep_traits) != 12 or len(set(sleep_traits)) != 12:
        raise ValueError("expected the locked 12-trait sleep panel")
    traits = {row["trait_id"]: row for row in trait_rows}
    if set(context_traits).difference(traits):
        raise ValueError("one or more locked context traits are absent from atlas traits")
    selected = [row for row in pair_rows if row["non_sleep_trait"] in context_traits]
    expected = {(sleep, context) for sleep in sleep_traits for context in context_traits}
    observed = {(row["sleep_trait"], row["non_sleep_trait"]) for row in selected}
    if len(selected) != len(expected) or observed != expected:
        raise ValueError(f"frozen context family mismatch: expected {len(expected)}, got {len(selected)}")
    output = []
    for row in selected:
        enriched = {field: row[field] for field in PAIR_FIELDS}
        trait = traits[row["non_sleep_trait"]]
        enriched.update({
            "context_label": trait["label"],
            **{field: trait[field] for field in TRAIT_FIELDS if field != "label"},
            "reuse_status": "READ_ONLY_FROZEN_ATLAS_RESULT",
            "fdr_denominator": "396 locked sleep-by-non-sleep pairs",
            "cohort_overlap_status": "Exact participant overlap not reassessed here; consult source-level overlap ledger; unknown is not zero",
            "claim_limit": "Global genetic correlation only; not independent replication, local sharing, colocalization, mechanism, or causality",
        })
        output.append(enriched)
    return output


def publish_idempotently(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise SystemExit(f"REFUSING_NONIDENTICAL_OVERWRITE: {path.relative_to(ROOT)}")
        return
    with path.open("xb") as handle:
        handle.write(payload)


def main() -> int:
    provenance_path = ROOT / "results/atlas/core.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("analysis_id") != "atlas-v1.0-downstream":
        raise SystemExit("UNEXPECTED_CORE_ATLAS_ANALYSIS_ID")
    for relative, expected_hash in provenance["outputs"].items():
        path = ROOT / relative
        if not path.is_file() or sha256(path) != expected_hash:
            raise SystemExit(f"CORE_ATLAS_OUTPUT_HASH_MISMATCH: {relative}")
    rg_matrix = ROOT / "results/tables/rg_matrix.tsv"
    if sha256(rg_matrix) != provenance["rg_matrix_sha256"]:
        raise SystemExit("CORE_ATLAS_RG_MATRIX_HASH_MISMATCH")
    panel_path = ROOT / "config/analysis_panel.tsv"
    lock_path = ROOT / "config/analysis_panel.lock.json"
    if sha256(panel_path) != provenance["panel_sha256"] or sha256(lock_path) != provenance["panel_lock_sha256"]:
        raise SystemExit("CORE_ATLAS_PANEL_LOCK_MISMATCH")

    panel = read_tsv(panel_path)
    sleep_traits = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
    traits_path = ROOT / "results/atlas/traits.tsv"
    pairs_path = ROOT / "results/atlas/trait_pairs.tsv"
    traits = read_tsv(traits_path)
    pairs = read_tsv(pairs_path)
    if len(pairs) != 396:
        raise SystemExit(f"FROZEN_ATLAS_PAIR_COUNT_MISMATCH: {len(pairs)}")
    selected = select_context_rows(sleep_traits, pairs, traits)
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected)
    payload = output.getvalue().encode("utf-8")
    extract_path = ROOT / "frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv"
    publish_idempotently(extract_path, payload)

    contexts = {row["trait_id"]: row for row in traits if row["trait_id"] in CONTEXT_TRAITS}
    manifest = {
        "artifact_role": "READ_ONLY_EXTRACT_OF_FROZEN_CORE_ATLAS_NOT_NEW_ANALYSIS",
        "atlas_analysis_id": provenance["analysis_id"],
        "atlas_source_commit": "660c9b5",
        "source_hashes": {
            "results/atlas/traits.tsv": sha256(traits_path),
            "results/atlas/trait_pairs.tsv": sha256(pairs_path),
            "results/tables/rg_matrix.tsv": sha256(rg_matrix),
            "config/analysis_panel.tsv": sha256(panel_path),
            "config/analysis_panel.lock.json": sha256(lock_path),
            "results/atlas/core.provenance.json": sha256(provenance_path),
        },
        "pair_count": len(selected),
        "context_traits": {trait: {field: contexts[trait][field] for field in TRAIT_FIELDS}
                           for trait in CONTEXT_TRAITS},
        "fdr_denominator": "396 locked sleep-by-non-sleep pairs; original q-values retained",
        "all396_q_le_0_05_count": sum(float(row["global_rg_fdr_all_396"]) <= 0.05 for row in selected),
        "extract": str(extract_path.relative_to(ROOT)),
        "extract_sha256": hashlib.sha256(payload).hexdigest(),
        "claim_limit": "Frozen global rg reuse only; no new replication or mechanistic inference.",
    }
    manifest_path = ROOT / "frailty_paper/manifests/frozen_atlas_aging_context_manifest_v2.json"
    runner_path = Path(__file__).resolve()
    manifest["reproducibility"] = {
        "generator": str(runner_path.relative_to(ROOT)),
        "generator_sha256": sha256(runner_path),
        "command": "python3 frailty_paper/scripts/29_audit_frozen_aging_context_rg.py",
        "python_version": sys.version.split()[0],
    }
    publish_idempotently(manifest_path, (json.dumps(manifest, indent=2) + "\n").encode())
    print(f"FROZEN_AGING_CONTEXT_OK pairs={len(selected)} q_le_0.05={manifest['all396_q_le_0_05_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
