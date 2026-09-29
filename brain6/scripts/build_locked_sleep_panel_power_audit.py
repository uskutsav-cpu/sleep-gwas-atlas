#!/usr/bin/env python3
"""Audit all locked global-map sleep traits without inventing missing LAVA results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/power_optimized_sensitivity_v1"
OUT = BASE / "all_locked_sleep_trait_power_audit_v1.tsv"
PROV = BASE / "all_locked_sleep_trait_power_audit_v1.provenance.json"
ATLAS_TRAITS = ROOT / "results/atlas/traits.tsv"
FAMILY_CONFIG = ROOT / "brain6/config/brain6_locked_family.yaml"
CANONICAL = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
CANONICAL_CAUSES = ROOT / "brain6/results/lava/canonical_v3_not_run_cells.tsv"
CANONICAL_DECISION = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json"
DURATION_SCREEN = BASE / "lava_trait_screen_v1/sleep_duration_continuous_lava_univariate.tsv"
DURATION_PROV = BASE / "provenance.json"
DURATION_COMPLETION_PROV = BASE / "completed_trait_only_screen_comparison_v1.provenance.json"
CANONICAL_POWER_AUDIT_PROV = BASE / "canonical_family_trait_power_audit_v1.provenance.json"
REFERENCE_ROOT = ROOT / "ref/lava/ukb_v1.1"
REFERENCE_PROV = REFERENCE_ROOT / "reference.provenance.json"
EXTERNAL_DATA = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies")

AUDIT_MODULE_PATH = ROOT / "brain6/scripts/build_canonical_family_trait_power_audit.py"
sys.path.insert(0, str(AUDIT_MODULE_PATH.parent))
spec = importlib.util.spec_from_file_location("canonical_trait_power_audit", AUDIT_MODULE_PATH)
if spec is None or spec.loader is None:
    raise ImportError(f"Cannot load reference-overlap implementation from {AUDIT_MODULE_PATH}")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def locked_sleep_traits(path: Path = FAMILY_CONFIG) -> list[str]:
    match = re.search(r"^sleep_traits:\s*\[([^]]+)\]\s*$", path.read_text(encoding="utf-8"), re.M)
    if not match:
        raise ValueError("Locked Brain6 family config has no parseable sleep_traits list")
    traits = [value.strip() for value in match.group(1).split(",")]
    if len(traits) != len(set(traits)) or not traits:
        raise ValueError("Locked sleep trait list is empty or duplicated")
    return traits


def quantile(values: list[float], probability: float) -> float:
    values = sorted(values)
    if not values:
        return math.nan
    at = (len(values) - 1) * probability
    low, high = math.floor(at), math.ceil(at)
    if low == high:
        return values[low]
    return values[low] * (high - at) + values[high] * (at - low)


def effective_n(ncase: str, ncontrol: str, n_total: str) -> str:
    if ncase.isdigit() and ncontrol.isdigit() and int(ncase) > 0 and int(ncontrol) > 0:
        cases, controls = int(ncase), int(ncontrol)
        return f"{4 / (1 / cases + 1 / controls):.1f}"
    return n_total if n_total.isdigit() else "NOT_REPORTED"


def summarize_cells(rows: list[dict[str, str]]) -> dict[str, Any]:
    if len(rows) != 2495:
        raise ValueError(f"Expected 2,495 locus rows; observed {len(rows)}")
    counts = Counter(row["status"] for row in rows)
    reasons = Counter(row.get("reason", "") for row in rows if row["status"] == "NOT_RUN")
    h2 = [float(row["h2.obs"]) for row in rows if row["status"] == "TESTED" and row.get("h2.obs") not in {"", "NA"}]
    return {
        "tested": counts["TESTED"], "not_run": counts["NOT_RUN"], "failed": counts["FAILED"],
        "low_h2": reasons["LOW_LOCAL_H2_UNDERPOWERED"],
        "shared_reference_min_k": reasons["FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS"],
        "other_min_k": reasons["FEWER_THAN_MIN_K"],
        "h2_median": statistics.median(h2) if h2 else math.nan,
        "h2_q25": quantile(h2, .25), "h2_q75": quantile(h2, .75),
        "h2_observed_n": len(h2),
    }


def build() -> dict[str, Any]:
    if OUT.exists() or PROV.exists():
        raise FileExistsError("Refusing to overwrite the immutable all-sleep-panel audit")
    traits = locked_sleep_traits()
    atlas_by_trait = {row["trait_id"]: row for row in read_tsv(ATLAS_TRAITS) if row["domain"] == "sleep"}
    if set(traits) != set(atlas_by_trait):
        raise ValueError("Locked sleep panel and atlas trait metadata differ")
    reference_prov = json.loads(REFERENCE_PROV.read_text(encoding="utf-8"))
    if reference_prov.get("verification") != "SHA-256 verified after official HTTPS acquisition":
        raise ValueError("LAVA reference provenance is not checksum sealed")
    reference_hashes = {row["path"]: row["sha256"] for row in reference_prov.get("extracted_files", [])
                        if row["path"].endswith(".info")}

    canonical_rows = read_tsv(CANONICAL)
    canonical_by_trait: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in canonical_rows:
        canonical_by_trait[row["phen"]].append(row)
    cause_by_trait: dict[str, Counter[str]] = defaultdict(Counter)
    for row in read_tsv(CANONICAL_CAUSES):
        cause_by_trait[row["trait_id"]][row["reason"]] += 1

    screen_rows = read_tsv(DURATION_SCREEN)
    if {row["trait_id"] for row in screen_rows} != {"sleep_duration_continuous_dashti_2019"}:
        raise ValueError("Continuous duration screen has an unexpected trait identity")
    screen = summarize_cells([{
        "status": row["status"], "h2.obs": row["h2_obs"], "reason": row["reason"]
    } for row in screen_rows])

    inputs_sha256: dict[str, str] = {}
    per_trait_sources: dict[str, str] = {}
    all_reference_hashes: dict[str, str] = {}
    output_rows = []
    for trait in traits:
        meta = atlas_by_trait[trait]
        relative_harmonized_path = meta["harmonized_path"]
        source = EXTERNAL_DATA / relative_harmonized_path
        if not source.is_file():
            raise FileNotFoundError(f"Locked global-map sleep input is unavailable: {source}")
        if source.stat().st_size <= 0:
            raise ValueError(f"Locked sleep input is empty: {source}")
        source_hash = sha256(source)
        if source_hash != meta["harmonized_sha256"]:
            raise ValueError(f"Locked global-map source hash mismatch: {trait}")
        total, overlap, chromosome_hashes = audit.source_reference_overlap(
            source, source.stat().st_size, source_hash, reference_hashes)
        if total != int(meta["harmonized_rows"]):
            raise ValueError(f"Harmonized variant count differs from the locked atlas for {trait}")
        per_trait_sources[str(source)] = source_hash
        all_reference_hashes.update({f"chr{chrom}": digest for chrom, digest in chromosome_hashes.items()})

        local: dict[str, Any] | None = None
        analysis_status = "NOT_ANALYZED_IN_CANONICAL_V3_OR_POWER_SCREEN"
        if trait in {"insomnia", "longsleep"}:
            local = summarize_cells(canonical_by_trait[trait])
            analysis_status = "CANONICAL_LAVA_V3_FAILED_QC_NOT_PROMOTED"
            cause = cause_by_trait[trait]
            if local["low_h2"] != cause["LOW_LOCAL_H2_UNDERPOWERED"]:
                raise ValueError(f"Canonical local-h2 cause mismatch for {trait}")
            if local["shared_reference_min_k"] != cause["FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS"]:
                raise ValueError(f"Canonical shared-reference cause mismatch for {trait}")
        elif trait == "sleepdur":
            local = screen
            analysis_status = "POWER_OPTIMIZED_TRAIT_ONLY_SCREEN_NO_PROMOTION"

        case_count = meta["ncase"]
        control_count = meta["ncontrol"]
        total_n = meta["n_total"]
        output_rows.append({
            "sleep_trait": trait,
            "phenotype_label": meta["label"],
            "phenotype_definition": meta["phenotype_definition"],
            "study_version": meta["dataset_version"],
            "source_id": meta["source_id"], "PMID": meta["pmid"], "DOI": meta["doi"],
            "source_sample_size": total_n,
            "cases": case_count, "controls": control_count,
            "effective_N_case_control_approx_or_quantitative_N": effective_n(case_count, control_count, total_n),
            "ancestry": meta["ancestry"], "genome_build": meta["build"],
            "atlas_source_status": meta["source_status"],
            "atlas_harmonized_unique_variants": total,
            "atlas_source_SHA256": source_hash,
            "exact_SNP_ID_reference_overlap": overlap,
            "reference_overlap_pct": f"{overlap / total * 100:.6f}",
            "reference_overlap_interpretation": (
                "ID/chromosome overlap only; build differs from frozen GRCh37 LAVA reference; coordinate harmonization required"
                if meta["build"].lower() in {"hg38", "grch38"} else
                "Exact lower-case SNP-ID/chromosome match only; not allele/effect harmonization"),
            "SNP_h2": meta["h2"], "SNP_h2_SE": meta["h2_se"], "SNP_h2_scale": meta["h2_scale"],
            "local_analysis_status": analysis_status,
            "planned_loci": 2495 if local is not None else "NOT_ANALYZED",
            "local_h2_tested_loci": local["tested"] if local is not None else "NOT_ANALYZED",
            "local_h2_tested_pct": f"{local['tested'] / 2495 * 100:.4f}" if local is not None else "NOT_ANALYZED",
            "local_h2_low_support_NOT_RUN": local["low_h2"] if local is not None else "NOT_ANALYZED",
            "shared_reference_minK_NOT_RUN": local["shared_reference_min_k"] if local is not None else "NOT_ANALYZED",
            "other_minK_NOT_RUN": local["other_min_k"] if local is not None else "NOT_ANALYZED",
            "local_h2_failed": local["failed"] if local is not None else "NOT_ANALYZED",
            "local_h2_NOT_RUN_total": local["not_run"] if local is not None else "NOT_ANALYZED",
            "local_h2_median_TESTED": f"{local['h2_median']:.10g}" if local is not None else "NOT_ANALYZED",
            "local_h2_q25_TESTED": f"{local['h2_q25']:.10g}" if local is not None else "NOT_ANALYZED",
            "local_h2_q75_TESTED": f"{local['h2_q75']:.10g}" if local is not None else "NOT_ANALYZED",
            "local_h2_distribution_scope": "TESTED cells only; NOT_RUN cells have no estimate" if local is not None else "No local-h2 analysis/receipt exists",
        })
    inputs = [ATLAS_TRAITS, FAMILY_CONFIG, CANONICAL, CANONICAL_CAUSES, CANONICAL_DECISION,
              DURATION_SCREEN, DURATION_PROV, DURATION_COMPLETION_PROV,
              CANONICAL_POWER_AUDIT_PROV, REFERENCE_PROV]
    for path in inputs:
        inputs_sha256[str(path.relative_to(ROOT))] = sha256(path)
    source_manifest = ROOT / "brain6/manifests/locked_dense_input_audit.tsv"
    inputs_sha256[str(source_manifest.relative_to(ROOT))] = sha256(source_manifest)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0])
    with OUT.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    prov = {
        "schema_version": 1,
        "analysis_id": "brain6_all_locked_sleep_trait_power_audit_v1",
        "status": "PASS_ALL_12_LOCKED_SLEEP_SOURCES_HASHED_AND_REFERENCE_AUDITED",
        "scope": "All twelve sleep phenotypes from the locked Brain6 12x6 global map; local-h2 fields are populated only for the two canonical LAVA v3 sleep inputs and the separate continuous-duration trait-only screen.",
        "selection_basis": "Trait order and names are parsed from brain6/config/brain6_locked_family.yaml; no phenotype selected using association results.",
        "local_h2_method": "Canonical v3 h2.obs is summarized only for TESTED cells. Continuous sleep duration is from the separate outcome-blinded 2,495-locus univariate screen. Other locked sleep traits are explicitly NOT_ANALYZED.",
        "reference_overlap_method": "Exact lower-case SNP-ID and chromosome intersection against the checksum-sealed UKB v1.1 reference. It does not harmonize alleles/effects; hg38 source overlap is a triage metric and cannot imply GRCh37 LAVA compatibility.",
        "inputs_sha256": inputs_sha256,
        "locked_sleep_source_sha256": per_trait_sources,
        "reference_chromosome_sha256": all_reference_hashes,
        "canonical_run_id": "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b",
        "output": {"path": str(OUT.relative_to(ROOT)), "rows": len(output_rows), "sha256": sha256(OUT)},
        "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "Only insomnia and long sleep are in the canonical v3 seven-input LAVA family; ten other global-map sleep traits lack canonical v3 local-h2 receipts.",
            "Continuous sleep duration has a separate trait-only screen, which failed its predeclared promotion rule and supplies no bivariate local-rg results.",
            "SNP-ID overlap is not allele/effect harmonization or a sufficient LAVA QC check.",
            "The sleep-apnea source is hg38 while the frozen LAVA reference is GRCh37; its ID-overlap count is not evidence that it is LAVA ready.",
            "Tested-locus local-h2 summaries do not impute values for NOT_RUN loci.",
        ],
    }
    with PROV.open("x", encoding="utf-8") as stream:
        json.dump(prov, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return {"status": prov["status"], "rows": len(output_rows), "source_hash_verified": len(per_trait_sources), "output_sha256": sha256(OUT)}


def validate() -> dict[str, Any]:
    prov = json.loads(PROV.read_text(encoding="utf-8"))
    if prov.get("builder_sha256") != sha256(Path(__file__).resolve()):
        raise ValueError("Builder source hash differs from audit provenance")
    if sha256(OUT) != prov.get("output", {}).get("sha256"):
        raise ValueError("All-sleep audit table checksum mismatch")
    for name, expected in prov.get("inputs_sha256", {}).items():
        path = ROOT / name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Locked input hash mismatch: {name}")
    for name, expected in prov.get("locked_sleep_source_sha256", {}).items():
        path = Path(name)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"External sleep source hash mismatch: {name}")
    for name, expected in prov.get("reference_chromosome_sha256", {}).items():
        path = REFERENCE_ROOT / f"lava-ukb-v1.1_{name}.info"
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Reference chromosome hash mismatch: {path}")
    rows = read_tsv(OUT)
    if len(rows) != 12 or sum(row["local_analysis_status"] == "CANONICAL_LAVA_V3_FAILED_QC_NOT_PROMOTED" for row in rows) != 2:
        raise ValueError("Expected all twelve locked sleep traits and exactly two canonical v3 sleep inputs")
    if any(row["local_analysis_status"] == "NOT_ANALYZED_IN_CANONICAL_V3_OR_POWER_SCREEN" and
           row["local_h2_tested_loci"] != "NOT_ANALYZED" for row in rows):
        raise ValueError("An unrun trait was assigned fabricated local-h2 coverage")
    return {"status": "PASS", "rows": len(rows), "output_sha256": sha256(OUT)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build() if args.write else validate(), indent=2, sort_keys=True))
