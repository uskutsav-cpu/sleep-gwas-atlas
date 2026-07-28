#!/usr/bin/env python3
"""Audit the Phase 0 trait registry before raw-data acquisition or LDSC.

The audit emits one row per trait so unresolved source facts are visible rather
than becoming implicit assumptions in a shell script. ``--strict`` exits
non-zero when a trait marked CURATED still has a blocking field.
"""
import argparse
import os

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "trait_id", "label", "domain", "type", "source_note", "raw_file",
    "ncase", "ncontrol", "n_total", "pop_prev", "build", "status",
}


def numeric(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/traits.tsv")
    parser.add_argument("--out", default="results/tables/phase0_curation_audit.tsv")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = pd.read_csv(args.config, sep="\t", dtype=str)
    missing = REQUIRED_COLUMNS.difference(config.columns)
    if missing:
        raise SystemExit(f"ERROR: config is missing required columns: {sorted(missing)}")
    if config["trait_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate trait_id values in config")

    rows = []
    for _, trait in config.iterrows():
        issues = []
        if not str(trait["source_note"]).strip():
            issues.append("source_note_missing")
        if not str(trait["raw_file"]).strip():
            issues.append("raw_file_missing")
        if trait["build"] not in {"hg19", "GRCh37"}:
            issues.append("not_hg19_requires_separate_liftover_decision")
        if trait["type"] == "binary":
            ncase, ncontrol = numeric(trait["ncase"]), numeric(trait["ncontrol"])
            if not ncase or not ncontrol:
                issues.append("binary_ncase_or_ncontrol_unresolved")
            prevalence = numeric(trait["pop_prev"])
            if prevalence is None or not 0 < prevalence < 1:
                issues.append("population_prevalence_unresolved")
            if "pop_prev_citation" not in config.columns or not str(trait.get("pop_prev_citation", "")).strip() or str(trait.get("pop_prev_citation", "")).upper() == "UNRESOLVED":
                issues.append("population_prevalence_citation_missing")
        elif not numeric(trait["n_total"]):
            issues.append("continuous_n_total_unresolved")
        if "pmid" not in config.columns or not str(trait.get("pmid", "")).strip():
            issues.append("pmid_not_linked_in_pipeline_config")
        if "ancestry" not in config.columns or str(trait.get("ancestry", "")).upper() != "EUR":
            issues.append("eur_subset_not_linked_in_pipeline_config")

        status = str(trait["status"]).upper()
        rows.append({
            "trait_id": trait["trait_id"], "label": trait["label"], "domain": trait["domain"],
            "declared_status": status, "curation_ready": not issues,
            "issues": ";".join(issues) if issues else "none",
        })

    audit = pd.DataFrame(rows)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    audit.to_csv(args.out, sep="\t", index=False)
    curated_with_issues = audit[(audit["declared_status"] == "CURATED") & ~audit["curation_ready"]]
    print(f"Wrote Phase 0 audit: {args.out}")
    print(f"Curation-ready rows: {int(audit['curation_ready'].sum())} / {len(audit)}")
    print(f"CURATED rows with blocking issues: {len(curated_with_issues)}")
    if args.strict and len(curated_with_issues):
        raise SystemExit("ERROR: CURATED traits have unresolved blocking curation fields")


if __name__ == "__main__":
    main()
