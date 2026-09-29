#!/usr/bin/env python3
"""Join reportable sleep overlap rows to the authoritative canonical audit.

The historical two-trait table used a nearest-rank quartile and omitted the
locked long-sleep SNP-h2 value. This builder preserves that table and writes a
corrected derivative using the canonical seven-input, receipt-backed audit.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/power_optimized_sensitivity_v1"
INPUT = BASE / "current_sleep_trait_audit_complete.tsv"
CANONICAL = BASE / "canonical_family_trait_power_audit_v1.tsv"
CANONICAL_PROVENANCE = BASE / "canonical_family_trait_power_audit_v1.provenance.json"
OVERLAP_PROVENANCE = BASE / "canonical_sleep_reference_overlap_v1.provenance.json"
OUTPUT = BASE / "current_sleep_trait_audit_complete_v2.tsv"
PROVENANCE = BASE / "current_sleep_trait_audit_complete_v2.provenance.json"
EXPECTED = {
    "insomnia": (6077635, 6061267, 1824, 671, 651),
    "longsleep": (6549769, 6057783, 1204, 1291, 1271),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def build() -> dict[str, object]:
    if OUTPUT.exists() or PROVENANCE.exists():
        raise FileExistsError("Refusing to overwrite the canonical-metrics v2 audit")
    fields, rows = read_tsv(INPUT)
    _, canonical_rows = read_tsv(CANONICAL)
    canonical = {row["trait_id"]: row for row in canonical_rows}
    if {row["trait"] for row in rows} != set(EXPECTED) or set(canonical) & set(EXPECTED) != set(EXPECTED):
        raise ValueError("Expected both canonical sleep traits")

    fields = [*fields, "local_h2_quartile_method", "snp_h2_source"]
    for row in rows:
        trait = row["trait"]
        source = canonical[trait]
        source_tuple = (
            int(row["usable_variant_count"]),
            int(row["frozen_reference_overlap_count"]),
            int(row["processable_tested"]),
            int(row["not_run"]),
            int(row["low_local_h2_not_run"]),
        )
        expected = EXPECTED[trait]
        canonical_tuple = (
            int(source["dense_unique_variants_after_frozen_QC"]),
            int(source["exact_SNP_ID_reference_overlap"]),
            int(source["tested"]),
            int(source["NOT_RUN"]),
            int(source["low_local_h2_NOT_RUN"]),
        )
        if source_tuple != expected or canonical_tuple != expected:
            raise ValueError(f"Sleep trait diagnostic counts disagree for {trait}")
        row["local_h2_median_tested"] = source["local_h2_median_tested"]
        row["local_h2_q25_tested"] = source["local_h2_q25_tested"]
        row["local_h2_q75_tested"] = source["local_h2_q75_tested"]
        row["snp_h2"] = f"{source['SNP_h2_reported_in_locked_global_map']} liability scale"
        row["local_h2_quartile_method"] = "Linear interpolation; TESTED loci only"
        row["snp_h2_source"] = "Locked global-map LDSC result; canonical audit Table S32 source"

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    canonical_prov = json.loads(CANONICAL_PROVENANCE.read_text(encoding="utf-8"))
    overlap_prov = json.loads(OVERLAP_PROVENANCE.read_text(encoding="utf-8"))
    if canonical_prov.get("canonical_run_id") != "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b":
        raise ValueError("Canonical audit is not bound to the immutable LAVA v3 run")
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_current_sleep_trait_audit_canonical_metrics_v2",
        "status": "PASS_CORRECTED_CANONICAL_METRICS",
        "scope": "Corrected reportable diagnostics for the two sleep traits in the seven-input immutable canonical LAVA v3 family",
        "method": "Preserve the exact-ID overlap and receipt counts from the prior report table; replace local-h2 median/quartiles and inherited SNP-h2 with values from the hash-bound seven-input canonical family audit. Quartiles are linear-interpolation quantiles over TESTED loci only.",
        "historical_input_preserved": str(INPUT.relative_to(ROOT)),
        "canonical_run_id": canonical_prov["canonical_run_id"],
        "inputs_sha256": {
            str(INPUT.relative_to(ROOT)): sha256(INPUT),
            str(CANONICAL.relative_to(ROOT)): sha256(CANONICAL),
            str(CANONICAL_PROVENANCE.relative_to(ROOT)): sha256(CANONICAL_PROVENANCE),
            str(OVERLAP_PROVENANCE.relative_to(ROOT)): sha256(OVERLAP_PROVENANCE),
        },
        "outputs_sha256": {OUTPUT.name: sha256(OUTPUT)},
        "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "trait_metrics": {
            row["trait"]: {
                "tested": int(row["processable_tested"]),
                "not_run": int(row["not_run"]),
                "snp_h2": row["snp_h2"],
                "local_h2_median_tested": row["local_h2_median_tested"],
                "local_h2_q25_tested": row["local_h2_q25_tested"],
                "local_h2_q75_tested": row["local_h2_q75_tested"],
            }
            for row in rows
        },
    }
    with PROVENANCE.open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return {"status": provenance["status"], "output": str(OUTPUT.relative_to(ROOT)), "traits": provenance["trait_metrics"]}


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
