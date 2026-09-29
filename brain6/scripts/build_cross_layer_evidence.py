#!/usr/bin/env python3
"""Build a source-bound, conservative evidence table for the five primary pairs."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6"
RULES_PATH = OUT / "config/cross_layer_evidence_rules_v1.json"
RULES_HASH_PATH = RULES_PATH.with_suffix(RULES_PATH.suffix + ".sha256")
OUTPUT_PATH = OUT / "results/evidence/brain6_cross_layer_evidence.tsv"
PROVENANCE_PATH = OUT / "results/evidence/brain6_cross_layer_evidence.provenance.json"
LAVA_RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
TRACKB_PAIR = "insomnia__adhd"
FIELDS = (
    "pair_id", "sleep_trait", "brain_disorder", "global_rg", "global_se", "global_p",
    "original_396_family_q", "original_396_family_status", "replication_class",
    "replication_rg", "replication_se", "replication_p", "replication_direction_concordant",
    "replication_overlap_caveat", "canonical_v3_local_rg_family_status", "canonical_v3_pairwise_status",
    "placo_pair_status", "placo_headline_variants", "placo_independent_loci_status",
    "placo_full_family_status", "fine_mapping_status", "trait_trait_coloc_status",
    "eqtl_sqtl_status", "regulatory_annotation_status", "cell_type_status",
    "overall_evidence_category", "interpretation_boundary",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def assign_category(*, global_significant: bool, local_promoted: bool,
                    replication_admissible: bool, variant_or_molecular_layer: bool,
                    any_evidence: bool) -> str:
    """Apply the frozen precedence without creating a numerical evidence score."""
    if global_significant and replication_admissible and local_promoted and variant_or_molecular_layer:
        return "MULTI_LAYER_SUPPORTED"
    if global_significant:
        return "GENETIC_ONLY"
    if local_promoted:
        return "LOCAL_ONLY"
    if replication_admissible:
        return "REPLICATION_ONLY"
    return "SUGGESTIVE" if any_evidence else "UNRESOLVED"


def immutable_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != content:
            raise FileExistsError(f"Refusing to replace existing cross-layer artifact: {path}")
        return
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_bytes(content)
    tmp.replace(path)


def build_rows() -> tuple[list[dict[str, str]], dict[str, str]]:
    rules = read_json(RULES_PATH)
    rules_sha = sha256(RULES_PATH)
    if RULES_HASH_PATH.read_text(encoding="utf-8").split()[0] != rules_sha:
        raise ValueError("Cross-layer category rules hash is not frozen/current")
    if rules.get("rule_id") != "brain6_cross_layer_evidence_rules_v1" or rules.get("status") != "FROZEN_BEFORE_CROSS_LAYER_TABLE_BUILD":
        raise ValueError("Cross-layer category rules are not in their expected frozen state")

    deep_path = OUT / "config/deep_tracks_v1.tsv"
    deep_hash_path = deep_path.with_suffix(deep_path.suffix + ".sha256")
    if deep_hash_path.read_text(encoding="utf-8").split()[0] != sha256(deep_path):
        raise ValueError("Deep-track selection checksum mismatch")
    selected = tsv(deep_path)
    selected = [row for row in selected if row["eligibility"] == "POST_ATLAS_PRIMARY"]
    if len(selected) != 5:
        raise ValueError(f"Expected five frozen primary pairs, found {len(selected)}")

    global_path = OUT / "results/global/brain6_72_locked.tsv"
    global_rows = {(row["sleep_trait"], row["brain_disorder"]): row for row in tsv(global_path)}
    replication_path = OUT / "results/replication/replication_master.tsv"
    replication_rows = {(row["sleep_trait"], row["brain_disorder"]): row
                        for row in tsv(replication_path)}
    placo_path = OUT / "results/placo/placo_master.tsv"
    placo_rows = {(row["sleep_trait"], row["brain_disorder"]): row for row in tsv(placo_path)}

    family_lock_path = OUT / "config/placo_family_v3/family_lock.json"
    family_lock = read_json(family_lock_path)
    if family_lock.get("protected_legacy_pair") != TRACKB_PAIR:
        raise ValueError("Protected PLACO pair differs from the frozen family lock")
    placo_pairs = set(family_lock["pairs"])
    lava_decision_path = ROOT / "work/lava-canonical-v3-production" / LAVA_RUN_ID / "canonical_family_decision.json"
    lava_decision = read_json(lava_decision_path)
    if (lava_decision.get("run_id") != LAVA_RUN_ID or
            lava_decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or
            lava_decision.get("pairwise_stage_authorized") is not False):
        raise ValueError("Canonical LAVA decision no longer enforces the frozen no-promotion boundary")

    records: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for track in selected:
        sleep, disorder = track["primary_sleep_trait"], track["disorder"]
        pair_id = f"{sleep}__{disorder}"
        key = (sleep, disorder)
        if key in seen or key not in global_rows:
            raise ValueError(f"Missing or duplicate global-map pair: {pair_id}")
        seen.add(key)
        global_row = global_rows[key]
        if global_row["significance_under_original_396_family"].lower() != "true":
            raise ValueError(f"Primary pair is not significant under its inherited family status: {pair_id}")
        replication = replication_rows.get(key)
        placo = placo_rows.get(key)
        if replication is None:
            raise ValueError(f"No replication ledger row for primary pair: {pair_id}")
        if placo:
            if pair_id not in placo_pairs:
                raise ValueError(f"Published PLACO pair is absent from the frozen pair lock: {pair_id}")
            if not placo["stage_status"].startswith("COMPLETE_QC_PASS"):
                raise ValueError(f"PLACO pair result is not pair-QC-passed: {pair_id}")
            result_path = ROOT / placo["output_path"]
            if sha256(result_path) != placo["output_sha256"]:
                raise ValueError(f"PLACO output hash mismatch: {pair_id}")
            if placo["family_lock_sha256"] != sha256(family_lock_path):
                raise ValueError(f"PLACO family lock hash mismatch: {pair_id}")
            placo_status = "PAIR_QC_PASS_PARTIAL_FAMILY"
            variants = placo["headline_variants"]
            loci_status = placo["independent_loci_status"]
        elif pair_id == TRACKB_PAIR:
            if pair_id in placo_pairs:
                raise ValueError("Protected Track B unexpectedly appears in the execution pairs")
            placo_status = "PROTECTED_TRACK_B_ABSENT_NOT_RERUN"
            variants = "NA"
            loci_status = "NOT_AVAILABLE"
        else:
            raise ValueError(f"Unexpected missing PLACO pair: {pair_id}")

        replication_class = replication["replication_class"]
        replication_admissible = replication_class in {
            "DIRECT_REPLICATION", "PROXY_REPLICATION", "DIRECTIONAL_REPLICATION",
        }
        category = assign_category(
            global_significant=True,
            local_promoted=False,
            replication_admissible=replication_admissible,
            variant_or_molecular_layer=False,
            any_evidence=True,
        )
        records.append({
            "pair_id": pair_id,
            "sleep_trait": sleep,
            "brain_disorder": disorder,
            "global_rg": global_row["rg"],
            "global_se": global_row["se"],
            "global_p": global_row["p"],
            "original_396_family_q": global_row["original_BH_FDR_q"],
            "original_396_family_status": "INHERITED_SIGNIFICANT",
            "replication_class": replication_class,
            "replication_rg": replication["replication_rg"],
            "replication_se": replication["replication_SE"],
            "replication_p": replication["replication_P"],
            "replication_direction_concordant": replication["direction_concordant"],
            "replication_overlap_caveat": replication["overlap"],
            "canonical_v3_local_rg_family_status": "FAILED_QC_NOT_PROMOTED",
            "canonical_v3_pairwise_status": "NOT_AVAILABLE_PAIRWISE_STAGE_PROHIBITED",
            "placo_pair_status": placo_status,
            "placo_headline_variants": variants,
            "placo_independent_loci_status": loci_status,
            "placo_full_family_status": "INCOMPLETE_PROTECTED_TRACK_B_ABSENT",
            "fine_mapping_status": "NOT_RUN_NO_COMPLETE_ELIGIBLE_UPSTREAM_FAMILY",
            "trait_trait_coloc_status": "NOT_RUN_NO_COMPLETE_ELIGIBLE_UPSTREAM_FAMILY",
            "eqtl_sqtl_status": "NOT_RUN_NO_PRIORITIZED_QC_PASSED_LOCUS",
            "regulatory_annotation_status": "NOT_RUN_NO_PRIORITIZED_QC_PASSED_LOCUS",
            "cell_type_status": "NOT_RUN_NO_PRIORITIZED_QC_PASSED_LOCUS",
            "overall_evidence_category": category,
            "interpretation_boundary": (
                "Inherited global genetic association; no promoted local-rg or molecular mechanism. "
                "PLACO is partial-family where available; directional replication is not exact replication."
            ),
        })
    return records, {
        str(deep_path.relative_to(ROOT)): sha256(deep_path),
        str(global_path.relative_to(ROOT)): sha256(global_path),
        str(replication_path.relative_to(ROOT)): sha256(replication_path),
        str(placo_path.relative_to(ROOT)): sha256(placo_path),
        str(family_lock_path.relative_to(ROOT)): sha256(family_lock_path),
        str(lava_decision_path.relative_to(ROOT)): sha256(lava_decision_path),
        str(RULES_PATH.relative_to(ROOT)): rules_sha,
    }


def render(rows: list[dict[str, str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def build() -> dict:
    rows, inputs = build_rows()
    output = render(rows)
    provenance = {
        "schema_version": 1,
        "status": "DESCRIPTIVE_PARTIAL_CROSS_LAYER_INTEGRATION",
        "rows": len(rows),
        "evidence_categories": {value: sum(row["overall_evidence_category"] == value for row in rows)
                                for value in sorted({row["overall_evidence_category"] for row in rows})},
        "inputs": inputs,
        "builder_script_sha256": sha256(Path(__file__).resolve()),
        "output_path": str(OUTPUT_PATH.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(output).hexdigest(),
        "interpretation": "Descriptive five-pair integration only. It does not promote partial PLACO, failed-QC LAVA, or unavailable downstream analyses and performs no new hypothesis test.",
    }
    provenance_bytes = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8")
    immutable_write(OUTPUT_PATH, output)
    immutable_write(PROVENANCE_PATH, provenance_bytes)
    return provenance


def main() -> None:
    print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
