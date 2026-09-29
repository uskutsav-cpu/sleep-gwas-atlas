#!/usr/bin/env python3
"""Audit and idempotently extract the frozen core-atlas frailty rg family."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROV_PATH = ROOT / "results/atlas/core.provenance.json"
TRAITS_PATH = ROOT / "results/atlas/traits.tsv"
PAIRS_PATH = ROOT / "results/atlas/trait_pairs.tsv"
PANEL_PATH = ROOT / "config/analysis_panel.tsv"
PANEL_LOCK_PATH = ROOT / "config/analysis_panel.lock.json"
RAW_FRAILTY_PATH = ROOT / "data/raw/frailty.txt.gz"
EXTRACT_PATH = ROOT / "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv"
MANIFEST_PATH = ROOT / "frailty_paper/manifests/frozen_atlas_frailty_manifest.json"
FIELDS = [
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
    "interpretation_status", "global_rg", "global_rg_se", "global_rg_z",
    "global_rg_p", "global_rg_fdr_all_396", "global_rg_fdr_primary_372",
    "global_rg_primary_significant", "effect_direction",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def publish_idempotently(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise SystemExit(f"REFUSING_NONIDENTICAL_OVERWRITE: {path.relative_to(ROOT)}")
        return
    try:
        with path.open("xb") as handle:
            handle.write(payload)
    except FileExistsError:
        if path.read_bytes() != payload:
            raise SystemExit(f"REFUSING_NONIDENTICAL_OVERWRITE: {path.relative_to(ROOT)}")


def main() -> int:
    provenance = json.loads(PROV_PATH.read_text(encoding="utf-8"))
    if provenance.get("analysis_id") != "atlas-v1.0-downstream":
        raise SystemExit("UNEXPECTED_CORE_ATLAS_ANALYSIS_ID")
    for relative, expected in provenance["outputs"].items():
        path = ROOT / relative
        if not path.is_file() or sha256(path) != expected:
            raise SystemExit(f"CORE_ATLAS_OUTPUT_HASH_MISMATCH: {relative}")
    rg_matrix = ROOT / "results/tables/rg_matrix.tsv"
    if sha256(rg_matrix) != provenance["rg_matrix_sha256"]:
        raise SystemExit("CORE_ATLAS_RG_MATRIX_HASH_MISMATCH")
    if sha256(PANEL_PATH) != provenance["panel_sha256"]:
        raise SystemExit("CORE_ATLAS_PANEL_HASH_MISMATCH")
    if sha256(PANEL_LOCK_PATH) != provenance["panel_lock_sha256"]:
        raise SystemExit("CORE_ATLAS_PANEL_LOCK_HASH_MISMATCH")

    panel = read_tsv(PANEL_PATH)
    sleep_traits = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
    if len(panel) != 45 or len(sleep_traits) != 12:
        raise SystemExit("CORE_ATLAS_PANEL_SHAPE_MISMATCH")

    trait_rows = [row for row in read_tsv(TRAITS_PATH) if row["trait_id"] == "frailty"]
    pair_rows = [
        row for row in read_tsv(PAIRS_PATH)
        if row["non_sleep_trait"] == "frailty"
    ]
    if len(trait_rows) != 1 or len(pair_rows) != 12:
        raise SystemExit("FROZEN_FRAILTY_RESULT_FAMILY_SHAPE_MISMATCH")
    if {row["sleep_trait"] for row in pair_rows} != set(sleep_traits):
        raise SystemExit("FROZEN_FRAILTY_SLEEP_FAMILY_MISMATCH")
    if any(row["analysis_tier"] != "PRIMARY_PHASE1" for row in pair_rows):
        raise SystemExit("FROZEN_FRAILTY_ANALYSIS_TIER_MISMATCH")

    import io
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows([{field: row[field] for field in FIELDS} for row in pair_rows])
    table_payload = output.getvalue().encode("utf-8")
    publish_idempotently(EXTRACT_PATH, table_payload)

    frailty = trait_rows[0]
    raw_sha = sha256(RAW_FRAILTY_PATH)
    with gzip.open(RAW_FRAILTY_PATH, "rb") as handle:
        for _ in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            pass
    manifest = {
        "artifact_role": "READ_ONLY_EXTRACT_OF_FROZEN_CORE_ATLAS_RESULTS_NOT_NEW_ANALYSIS",
        "atlas_analysis_id": provenance["analysis_id"],
        "source_commit": "660c9b5",
        "source_panel_lock_sha256": provenance["panel_lock_sha256"],
        "source_output_sha256": {
            **provenance["outputs"],
            "results/tables/rg_matrix.tsv": provenance["rg_matrix_sha256"],
        },
        "source_h2_summary_sha256_expected": provenance["h2_summary_sha256"],
        "source_h2_summary_present": (ROOT / "results/tables/h2_summary.tsv").is_file(),
        "source_harmonized_frailty_path": frailty["harmonized_path"],
        "source_harmonized_frailty_sha256_expected": frailty["harmonized_sha256"],
        "source_harmonized_frailty_present": (ROOT / frailty["harmonized_path"]).is_file(),
        "source_munged_frailty_path": frailty["munged_path"],
        "source_munged_frailty_present": (ROOT / frailty["munged_path"]).is_file(),
        "raw_frailty_path": str(RAW_FRAILTY_PATH.relative_to(ROOT)),
        "raw_frailty_sha256_recomputed": raw_sha,
        "trait_table_frailty_record": {
            key: frailty[key] for key in (
                "trait_id", "source_id", "dataset_version", "ancestry", "build",
                "n_total", "harmonized_path", "harmonized_sha256", "harmonized_rows",
                "h2", "h2_se", "h2_z", "ldsc_intercept", "ldsc_ratio", "phase1_verdict",
            )
        },
        "frailty_pair_count": len(pair_rows),
        "fdr_all_396_significant_count": sum(
            float(row["global_rg_fdr_all_396"]) <= 0.05 for row in pair_rows
        ),
        "fdr_denominator": "396 locked sleep-by-non-sleep pairs",
        "extracted_table": str(EXTRACT_PATH.relative_to(ROOT)),
        "extracted_table_sha256": hashlib.sha256(table_payload).hexdigest(),
        "claim_limit": (
            "Reused frozen core atlas estimates; not independent replication, local sharing, "
            "pleiotropy, colocalization, mechanism, or causality."
        ),
    }
    manifest_payload = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
    publish_idempotently(MANIFEST_PATH, manifest_payload)
    print(
        "FROZEN_FRAILTY_AUDIT_OK rows=12 all396_fdr_significant="
        f"{manifest['fdr_all_396_significant_count']} raw_sha256={raw_sha}"
    )
    print(
        "SOURCE_REPROCESSING="
        + (
            "AVAILABLE"
            if manifest["source_h2_summary_present"]
            and manifest["source_harmonized_frailty_present"]
            and manifest["source_munged_frailty_present"]
            else "INCOMPLETE_MISSING_H2_OR_HARMONIZED_OR_MUNGED_INPUT"
        ),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
