#!/usr/bin/env python3
"""Freeze the existing, post-atlas Brain6 pair decisions and source ledger."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6"
LOCK = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
DECISIONS = ROOT / "extensions/brain6/work/overnight-v03/pair-decisions.reviewed.tsv"
PAIR_B = ROOT / "results/track_b/replication/pair_b_ldsc.tsv"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_new(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite {path}")
    path.write_text(text, encoding="utf-8")


def main() -> None:
    lock = json.loads(LOCK.read_text())
    decisions = pd.read_csv(DECISIONS, sep="\t", dtype=str).fillna("")
    global_map = pd.read_csv(OUT / "results/global/brain6_72_locked.tsv", sep="\t")
    source_registry = pd.read_csv(ROOT / "config/public_gwas_sources.tsv", sep="\t", dtype=str).fillna("")
    panel = pd.read_csv(ROOT / "config/analysis_panel.tsv", sep="\t", dtype=str).fillna("")
    replication = pd.read_csv(PAIR_B, sep="\t")
    if lock.get("selection_is_post_global_screen") is not True or lock.get("reviewed_by") != "Zunpeng":
        raise ValueError("Existing selection lock lacks post-atlas/reviewer provenance")
    if len(global_map) != 72 or len(decisions) != 6:
        raise ValueError("Unexpected locked global map or selection-family size")

    selected = []
    for _, d in decisions.iterrows():
        trait = d.sleep_trait
        disorder = d.disease_trait
        row = global_map.loc[
            global_map.brain_disorder.eq(disorder) & global_map.sleep_trait.eq(trait)
        ] if trait else global_map.iloc[0:0]
        if d.role == "PRIMARY" and (len(row) != 1 or not bool(row.iloc[0].significance_under_original_396_family)):
            raise ValueError(f"Primary pair is absent or not significant in original family: {disorder}/{trait}")
        selected.append({
            "disorder": disorder,
            "primary_sleep_trait": trait or "NO_ELIGIBLE_PAIR",
            "secondary_sleep_trait": "NA",
            "selection_reason": d.reason,
            "global_rg": row.iloc[0].rg if len(row) else "NA",
            "global_SE": row.iloc[0].se if len(row) else "NA",
            "global_P": row.iloc[0].p if len(row) else "NA",
            "original_q": row.iloc[0].original_BH_FDR_q if len(row) else "NA",
            "dense_GWAS_available": "NO_CURRENT_LOCAL_FILES; source/history only",
            "replication_available": (
                "DIRECTIONAL_REPLICATION; FinnGen R13 ADHD; cohort and phenotype caveats"
                if disorder == "adhd" else "NOT_ESTABLISHED_FOR_THIS_PAIR"
            ),
            "overlap_concern": (
                "UNKNOWN; individual-level overlap not verified"
                if disorder in {"adhd", "scz", "bipolar"} else
                "DOCUMENTED_POSSIBLE_UKB_OVERLAP; participant-level overlap unverified"
            ),
            "eligibility": "NO_ELIGIBLE_PAIR" if d.role == "NO_ELIGIBLE_PAIR" else "POST_ATLAS_PRIMARY",
            "notes": (
                "No pair passes original 396-family FDR; retain null in the six-disorder universe."
                if disorder == "alz" else
                "Post-atlas follow-up selection; original atlas FDR is inherited, not recalculated."
            ),
        })
    fields = list(selected[0])
    from io import StringIO
    buf = StringIO()
    writer = csv.DictWriter(buf, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected)
    target = OUT / "config/deep_tracks_v1.tsv"
    write_new(target, buf.getvalue())
    write_new(target.with_suffix(".tsv.sha256"), f"{digest(target)}  {target.name}\n")

    traits = panel.set_index("trait_id")
    registry = source_registry.set_index("source_id")
    ids = ["jansen_2019_insomnia_ukb", "dashti_2019_long_sleep", "demontis_2023_adhd_eur",
           "howard_2019_mdd_no23andme", "trubetskoy_2022_scz_eur", "mullins_2021_bipolar_eur",
           "nalls_2019_parkinson_public_proxy"]
    trait_for_source = panel.reset_index().set_index("source_id")["trait_id"].to_dict()
    ledger = []
    for source_id in ids:
        trait = trait_for_source[source_id]
        p, s = traits.loc[trait], registry.loc[source_id]
        ledger.append({
            "trait": trait, "study": p.source_note, "DOI": p.doi, "PMID": p.pmid,
            "GWAS_accession": "UNRESOLVED_IN_LOCKED_PANEL", "source_URL": s.download_url,
            "source_page_URL": s.source_page_url, "download_date": "UNKNOWN_FROM_CURRENT_LEDGER",
            "sample_size": p.n_total, "cases": p.ncase, "controls": p.ncontrol,
            "ancestry": p.ancestry, "age": "NOT_RECORDED", "sex_composition": "NOT_RECORDED",
            "cohorts": p.phenotype_definition, "UK_Biobank_overlap": "UNKNOWN; see cohort_overlap.tsv",
            "FinnGen_overlap": "UNKNOWN; see cohort_overlap.tsv", "23andMe_overlap": "UNKNOWN; see cohort_overlap.tsv",
            "PGC_overlap": "UNKNOWN; see cohort_overlap.tsv", "genome_build": p.build,
            "effect_allele_convention": "SOURCE_SPECIFIC; original dense file absent; see study source documentation",
            "effect_metric": "SOURCE_SPECIFIC; original dense file absent",
            "variant_count": "NOT_AVAILABLE_CURRENT_WORKTREE", "beta": "NOT_VERIFIED_CURRENT_FILE_MISSING",
            "SE": "NOT_VERIFIED_CURRENT_FILE_MISSING", "P": "NOT_VERIFIED_CURRENT_FILE_MISSING",
            "EAF": "NOT_VERIFIED_CURRENT_FILE_MISSING", "per_variant_N": "NOT_VERIFIED_CURRENT_FILE_MISSING",
            "INFO": "SOURCE_SPECIFIC_OR_NOT_RECORDED", "license": s.access,
            "expected_source_archive_bytes": s.archive_bytes, "expected_source_archive_SHA256": s.archive_sha256,
            "raw_file_path": "MISSING_FROM_CURRENT_WORKTREE", "raw_file_SHA256": "NOT_AVAILABLE",
            "current_readiness": "BLOCKED_MISSING_DENSE_FILE",
        })
    rep = replication.iloc[0]
    ledger.append({
        "trait": "adhd_replication_finngen_r13", "study": "FinnGen R13 F5 ADHD register phenotype",
        "DOI": "NA", "PMID": "NA", "GWAS_accession": "finngen_r13_F5_ADHD",
        "source_URL": "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_ADHD.gz?generation=1777989549884404",
        "source_page_URL": "https://www.finngen.fi/en/access_results", "download_date": "2026-09-01 (sealed receipt)",
        "sample_size": 495052, "cases": 5559, "controls": 489493, "ancestry": "FINNISH_EUR",
        "age": "NOT_RECORDED", "sex_composition": "NOT_RECORDED", "cohorts": "FinnGen R13",
        "UK_Biobank_overlap": "NONE_KNOWN; distinct FinnGen cohort; participant-level verification not supplied",
        "FinnGen_overlap": "DISCOVERY_OVERLAP_NOT_APPLICABLE_TO_EXTERNAL_REPLICATION",
        "23andMe_overlap": "UNKNOWN", "PGC_overlap": "UNKNOWN", "genome_build": "GRCh38",
        "effect_allele_convention": "FinnGen source release; harmonized to HM3 for LDSC",
        "effect_metric": "log odds / LDSC observed scale", "variant_count": 21326959,
        "beta": "YES_IN_SOURCE_RELEASE", "SE": "YES_IN_SOURCE_RELEASE", "P": "YES_IN_SOURCE_RELEASE",
        "EAF": "YES_IN_SOURCE_RELEASE", "per_variant_N": "CASE_CONTROL_N_IN_SOURCE",
        "INFO": ">0.6 source-level filter", "license": "PUBLIC",
        "expected_source_archive_bytes": 802195177, "expected_source_archive_SHA256": "05487573a3d954baebe5605042f119f2949f238705c865efbf872fd6e28e02a1",
        "raw_file_path": "NOT_RETAINED; replication ingestion/QC remains archived",
        "raw_file_SHA256": "05487573a3d954baebe5605042f119f2949f238705c865efbf872fd6e28e02a1",
        "current_readiness": "REPLICATION_RESULT_ARCHIVED; RAW_FILE_NOT_PRESENT",
    })
    columns = list(ledger[0])
    buf = StringIO()
    writer = csv.DictWriter(buf, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(ledger)
    gwas = OUT / "manifests/gwas_master.tsv"
    write_new(gwas, buf.getvalue())

    overlap = [
        ["insomnia__adhd", "UNKNOWN", "Participant-level discovery overlap was not verified; external FinnGen replication is a separate estimate."],
        ["insomnia__mdd", "DOCUMENTED", "Both discovery sources include UK Biobank; individual-level overlap and covariance are unverified."],
        ["longsleep__scz", "UNKNOWN", "No participant-level overlap evidence in the locked source ledger."],
        ["longsleep__bipolar", "UNKNOWN", "No participant-level overlap evidence in the locked source ledger."],
        ["longsleep__parkinson", "DOCUMENTED", "Sleep GWAS uses UK Biobank; PD GWAS includes UK Biobank proxy cases; individual overlap is unverified."],
    ]
    buf = StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(["pair", "overlap_classification", "evidence_and_limit"])
    writer.writerows(overlap)
    cohort = OUT / "manifests/cohort_overlap.tsv"
    write_new(cohort, buf.getvalue())
    print(json.dumps({"deep_tracks": str(target), "deep_tracks_sha256": digest(target),
                      "gwas_rows": len(ledger), "gwas_manifest": str(gwas),
                      "cohort_overlap": str(cohort), "pair_b_rg": float(rep.rg)}, indent=2))


if __name__ == "__main__":
    main()
