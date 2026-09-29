#!/usr/bin/env python3
"""Build an outcome-blind, provenance-bound screen for higher-power sleep GWAS."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from statistics import median, quantiles

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "frailty_paper/analysis/power_optimized_sleep_sensitivity_v1"
RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
CANONICAL_RUN = ROOT / "work/lava-canonical-v3-production" / RUN_ID
INPUTS = {
    "sleep_panel": ROOT / "config/analysis_panel.tsv",
    "source_scan": ROOT / "frailty_paper/manifests/gwas_source_scan_sleep_panel.tsv",
    "sleep_h2": ROOT / "frailty_paper/results/frailty_v1/h2_sleep_panel.tsv",
    "resources": ROOT / "frailty_paper/manifests/all_acquired_resources.tsv",
    "canonical_family": CANONICAL_RUN / "canonical_family_manifest.tsv",
    "canonical_results": CANONICAL_RUN / "results/canonical_family_results.tsv",
    "canonical_not_run": ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv",
    "brain6_gwas": ROOT / "brain6/manifests/gwas_master.tsv",
    "execution_config": ROOT / "brain6/config/lava_execution_canonical_v3.json",
    "frailty_lava_lock": ROOT / "frailty_paper/config/lava_frailty_sensitivity_v1.yaml",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def case_control_neff(cases: int, controls: int) -> float:
    """GWAS effective N convention: 4 / (1/Ncase + 1/Ncontrol)."""
    if cases < 1 or controls < 1:
        raise ValueError("cases and controls must both be positive")
    return 4.0 / (1.0 / cases + 1.0 / controls)


def eligibility_change(baseline_eligible: int, candidate_eligible: int, denominator: int) -> dict[str, float | bool]:
    if denominator <= 0 or not 0 <= baseline_eligible <= denominator or not 0 <= candidate_eligible <= denominator:
        raise ValueError("eligible counts must be within the positive denominator")
    baseline_not_run = denominator - baseline_eligible
    candidate_not_run = denominator - candidate_eligible
    absolute_pp = 100.0 * (candidate_eligible - baseline_eligible) / denominator
    relative_not_run_reduction = (
        (baseline_not_run - candidate_not_run) / baseline_not_run if baseline_not_run else 0.0
    )
    return {
        "baseline_eligible_percent": 100.0 * baseline_eligible / denominator,
        "candidate_eligible_percent": 100.0 * candidate_eligible / denominator,
        "absolute_eligibility_change_percentage_points": absolute_pp,
        "relative_not_run_reduction": relative_not_run_reduction,
        "passes_promotion_gate": absolute_pp >= 10.0 or relative_not_run_reduction >= 0.25,
    }


def _quantile(values: list[float], q: int) -> float | None:
    if not values:
        return None
    if len(values) < 2:
        return values[0]
    return quantiles(values, n=4, method="inclusive")[q - 1]


def build() -> dict[str, object]:
    missing = [str(path) for path in INPUTS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing screening input(s): {missing}")
    panel = read_tsv(INPUTS["sleep_panel"])
    panel = [row for row in panel if row["domain"] == "sleep"]
    if len(panel) != 12 or len({row["trait_id"] for row in panel}) != 12:
        raise ValueError("locked source panel must contain exactly 12 unique sleep traits")

    source_scan = {row["resource_id"].removeprefix("sleep_panel_"): row
                   for row in read_tsv(INPUTS["source_scan"])}
    h2 = {row["trait"]: row for row in read_tsv(INPUTS["sleep_h2"])}
    resources = {row["resource_id"].removeprefix("sleep_panel_"): row
                 for row in read_tsv(INPUTS["resources"])
                 if row["resource_id"].startswith("sleep_panel_")}

    canonical_rows = read_tsv(INPUTS["canonical_results"])
    if len(canonical_rows) != 2495 * 7:
        raise ValueError(f"canonical v3 aggregate has {len(canonical_rows)} rows, expected 17,465")
    canonical: dict[str, list[dict[str, str]]] = {}
    for row in canonical_rows:
        canonical.setdefault(row["phen"], []).append(row)
    family = read_tsv(INPUTS["canonical_family"])
    family_traits = sorted({r["trait_id"] for r in family})
    if len(family_traits) != 7:
        raise ValueError("canonical v3 family must contain seven unique input GWAS")
    trait_domains = {row["trait_id"]: row["domain"] for row in read_tsv(ROOT / "config/analysis_panel.tsv")}

    trait_rows: list[dict[str, object]] = []
    for row in panel:
        tid = row["trait_id"]
        scan = source_scan.get(tid, {})
        resource = resources.get(tid, {})
        h = h2.get(tid, {})
        cells = canonical.get(tid, [])
        tested = [x for x in cells if x["status"] == "TESTED"]
        hs = [float(x["h2.obs"]) for x in tested if x.get("h2.obs") not in ("", "NA")]
        snps = [int(x["n_snps"]) for x in tested if x.get("n_snps") not in ("", "NA")]
        low_h2 = sum(x["reason"] == "LOW_LOCAL_H2_UNDERPOWERED" for x in cells)
        fewer_k = sum(x["reason"] == "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS" for x in cells)
        fewer_min = sum(x["reason"] == "FEWER_THAN_MIN_K" for x in cells)
        ncase, ncontrol = row["ncase"], row["ncontrol"]
        if ncase not in ("NA", "") and ncontrol not in ("NA", ""):
            neff: float | str = round(case_control_neff(int(ncase), int(ncontrol)), 1)
        else:
            neff = row["n_total"] if row["type"] == "continuous" else "NA"
        trait_rows.append({
            "trait_id": tid, "label": row["label"], "phenotype_definition": row["phenotype_definition"],
            "source_note": row["source_note"], "PMID": row["pmid"], "DOI": row["doi"],
            "source_id": row["source_id"], "sample_size": row["n_total"], "cases": ncase,
            "controls": ncontrol, "effective_N_4_over_reciprocal_sum": neff,
            "ancestry": row["ancestry"], "genome_build": row["build"],
            "source_variant_rows_not_unique_qc_count": scan.get("data_rows", "NOT_SCANNED"),
            "source_scan_status": scan.get("stream_status", "NOT_SCANNED"),
            "source_file": resource.get("file", row["raw_file"]),
            "source_bytes": resource.get("bytes", "UNKNOWN"),
            "source_sha256": resource.get("sha256", "UNKNOWN"),
            "global_h2": h.get("h2", "NOT_AVAILABLE"), "global_h2_scale": h.get("scale", "NOT_AVAILABLE"),
            "global_h2_se": h.get("se", "NOT_AVAILABLE"), "global_h2_z": h.get("z", "NOT_AVAILABLE"),
            "canonical_v3_scope": "IN_FAMILY" if tid in canonical else "NOT_IN_7_INPUT_FAMILY",
            "canonical_v3_tested": len(tested) if cells else "NA",
            "canonical_v3_not_run": sum(x["status"] == "NOT_RUN" for x in cells) if cells else "NA",
            "canonical_v3_failed": sum(x["status"] == "FAILED" for x in cells) if cells else "NA",
            "canonical_v3_low_h2_not_run": low_h2 if cells else "NA",
            "canonical_v3_reference_min_k_not_run": fewer_k if cells else "NA",
            "canonical_v3_min_k_not_run": fewer_min if cells else "NA",
            "canonical_v3_tested_local_h2_median": round(median(hs), 9) if hs else "NA",
            "canonical_v3_tested_local_h2_q25": round(_quantile(hs, 1), 9) if hs else "NA",
            "canonical_v3_tested_local_h2_q75": round(_quantile(hs, 3), 9) if hs else "NA",
            "canonical_v3_tested_locus_n_snps_median": round(median(snps), 1) if snps else "NA",
            "frailty_lava_local_h2": "RUNNING_PAIRWISE_SENSITIVITY_NOT_CANONICAL_UNIVARIATE",
            "note": "Source row counts are not unique/fully QC-passed variant counts; exact candidate LD-reference overlap remains unmeasured.",
        })

    v3_rows = []
    brain6 = {r["trait"]: r for r in read_tsv(INPUTS["brain6_gwas"])}
    not_run = {r["trait_id"]: r for r in read_tsv(INPUTS["canonical_not_run"])}
    for tid in family_traits:
        cells = canonical[tid]
        tested = [x for x in cells if x["status"] == "TESTED"]
        hvals = [float(x["h2.obs"]) for x in tested if x.get("h2.obs") not in ("", "NA")]
        source = brain6.get(tid, {})
        counts = not_run.get(tid, {})
        v3_rows.append({
            "trait_id": tid, "domain_in_locked_atlas": trait_domains.get(tid, "NOT_IN_PANEL"),
            "study": source.get("study", ""), "PMID": source.get("PMID", ""), "DOI": source.get("DOI", ""),
            "sample_size": source.get("sample_size", ""), "cases": source.get("cases", ""),
            "controls": source.get("controls", ""), "ancestry": source.get("ancestry", ""),
            "genome_build": source.get("genome_build", ""),
            "source_variant_count": source.get("variant_count", "NOT_VERIFIED_IN_CURRENT_MANIFEST"),
            "planned_loci": len(cells), "tested": len(tested), "not_run": len(cells)-len(tested),
            "failed": sum(x["status"] == "FAILED" for x in cells),
            "low_local_h2": counts.get("low_local_h2_underpowered", "NOT_RECORDED"),
            "fewer_than_min_k_shared_reference_variants": counts.get("fewer_than_min_k_shared_reference_variants", "NOT_RECORDED"),
            "fewer_than_min_k": counts.get("fewer_than_min_k", "NOT_RECORDED"),
            "tested_local_h2_median": round(median(hvals), 9) if hvals else "NA",
            "tested_local_h2_q25": round(_quantile(hvals, 1), 9) if hvals else "NA",
            "tested_local_h2_q75": round(_quantile(hvals, 3), 9) if hvals else "NA",
            "tested_locus_n_snps_median": round(median([int(x["n_snps"]) for x in tested if x.get("n_snps") not in ("", "NA")]), 1),
            "interpretation": "Domain follows locked panel; only insomnia and long sleep are sleep phenotypes in this 7-input Brain6 v3 family." if tid in {"insomnia", "longsleep"} else "Brain-disorder phenotype, not a sleep trait.",
        })

    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "current_sleep_trait_audit.tsv", trait_rows)
    write_tsv(OUT / "brain6_v3_seven_input_diagnostics.tsv", v3_rows)
    provenance = {
        "analysis_id": "power_optimized_sleep_sensitivity_v1_screen",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "scope": "12 locked sleep traits, plus an explicit scope reconciliation for the 7-input Brain6 canonical v3 family",
        "association_outcomes_accessed": False,
        "inputs_sha256": {name: sha256(path) for name, path in INPUTS.items()},
        "outputs_sha256": {
            "frailty_paper/analysis/power_optimized_sleep_sensitivity_v1/current_sleep_trait_audit.tsv": sha256(OUT / "current_sleep_trait_audit.tsv"),
            "frailty_paper/analysis/power_optimized_sleep_sensitivity_v1/brain6_v3_seven_input_diagnostics.tsv": sha256(OUT / "brain6_v3_seven_input_diagnostics.tsv"),
        },
        "limitations": [
            "No per-variant candidate-to-LAVA-reference overlap has yet been independently counted.",
            "Frailty local-LAVA pairwise receipts are ongoing and are not treated as canonical independent univariate h2 results.",
            "Canonical local-h2 quantiles summarize tested loci only; low-h2 loci are recorded by frozen reason and do not have an h2 estimate.",
            "Source row counts are not asserted to be unique or fully QC-passed variant counts.",
        ],
    }
    (OUT / "current_sleep_trait_audit.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return provenance


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
