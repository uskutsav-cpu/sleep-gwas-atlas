#!/usr/bin/env python3
"""Build a source-level, non-replication context table for Xue et al. 2026."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "brain6/manifests/gwas_master.tsv"
OUTPUT = ROOT / "brain6/results/novelty/xue_2026_insomnia_psychiatric_context.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")
SOURCE = {
    "source_id": "xue_2026_insomnia_psychiatric",
    "citation": "Xue B, Niu M, Sun Y, et al. Sleep. 2026;49(1):zsaf317.",
    "pmid": "41065713",
    "doi": "10.1093/sleep/zsaf317",
    "url": "https://academic.oup.com/sleep/article-abstract/49/1/zsaf317/8279895",
    "pubmed_url": "https://pubmed.ncbi.nlm.nih.gov/41065713/",
    "online_publication_date": "2025-10-09",
    "issue": "2026-01",
    "reported_insomnia_neff": 314149,
    "method": "Bivariate MiXeR; conjFDR and ASSET cross-trait loci; MAGMA; SEISMIC",
    "reported_significant_insomnia_correlations": "ADHD; bipolar disorder; MDD; schizophrenia",
    "source_summary_limit": "The public abstract reports significance and per-disorder shared-locus counts, not pair-specific rg/SE/P or GWAS release identities.",
    "overlap_limit": "Source GWAS and participant identity with Brain6 are unverified; similar effective N is a compatibility clue, not proof of identical samples.",
}
DISORDERS = [
    ("adhd", 7),
    ("bipolar", 5),
    ("mdd", 15),
    ("scz", 19),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def input_metadata() -> dict[str, str | int]:
    with INPUT.open(encoding="utf-8", newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r.get("trait") == "insomnia"]
    if len(rows) != 1:
        raise ValueError(f"Expected one insomnia source row in {INPUT}, found {len(rows)}")
    row = rows[0]
    cases, controls = int(row["cases"]), int(row["controls"])
    # Standard effective sample-size approximation for a case-control GWAS.
    effective_n = round(4 / (1 / cases + 1 / controls))
    return {
        "brain6_gwas_master_sha256": sha256(INPUT),
        "brain6_insomnia_study": row["study"],
        "brain6_insomnia_pmid": row["PMID"],
        "brain6_insomnia_doi": row["DOI"],
        "brain6_insomnia_cases": cases,
        "brain6_insomnia_controls": controls,
        "brain6_insomnia_neff_from_case_control_counts": effective_n,
        "neff_formula": "round(4 / (1/cases + 1/controls))",
    }


def make_rows(meta: dict[str, str | int]) -> list[dict[str, str | int]]:
    rows = []
    for disorder, shared_loci in DISORDERS:
        rows.append({
            "source_id": SOURCE["source_id"],
            "citation": SOURCE["citation"],
            "pmid": SOURCE["pmid"],
            "doi": SOURCE["doi"],
            "online_publication_date": SOURCE["online_publication_date"],
            "issue": SOURCE["issue"],
            "sleep_trait": "insomnia",
            "brain6_disorder": disorder,
            "source_reports_significant_genetic_correlation": "true",
            "source_reported_shared_locus_count_total": 70,
            "source_reported_candidate_snp_count_total": 97,
            "source_reported_novel_association_count_by_disorder": shared_loci,
            "source_insomnia_neff": SOURCE["reported_insomnia_neff"],
            "brain6_insomnia_neff_from_locked_counts": meta["brain6_insomnia_neff_from_case_control_counts"],
            "neff_absolute_difference": abs(SOURCE["reported_insomnia_neff"] - int(meta["brain6_insomnia_neff_from_case_control_counts"])),
            "brain6_insomnia_cases": meta["brain6_insomnia_cases"],
            "brain6_insomnia_controls": meta["brain6_insomnia_controls"],
            "brain6_source_pmid": meta["brain6_insomnia_pmid"],
            "source_release_identity": "UNVERIFIED",
            "participant_overlap": "UNKNOWN; similar Neff alone does not establish identity",
            "pairwise_rg_se_p_available_in_reviewed_source": "not in abstract; full-text/supplement coefficients not source-verified here",
            "locus_level_match_to_brain6": "NOT_PERFORMED; source locus identifiers/results not independently extracted",
            "independent_replication": "NO; source overlap is unresolved and likely possible",
            "interpretation": "Prior published insomnia–psychiatric context only; no Brain6 estimate, test, locus, or novelty label is changed",
        })
    return rows


def tsv_bytes(rows: list[dict[str, str | int]]) -> bytes:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode()


def expected_provenance(rows: list[dict[str, str | int]], output_bytes: bytes, meta: dict[str, str | int], script: Path) -> dict:
    return {
        "schema_version": 1,
        "status": "PASS_SOURCE_LEVEL_PRIOR_CONTEXT_NOT_REPLICATION",
        "decision": "ADD_CONTEXT_ONLY; DO_NOT_CHANGE_LOCKED_BRAIN6_NOVELTY_OR_FDR",
        "scope": "Four insomnia–psychiatric pairs reported as significant in Xue et al. 2026; source-level comparison only.",
        "source": SOURCE,
        "source_accessed_date": "2026-09-26",
        "transcription_basis": "Publisher/PubMed abstract metadata; source-specific supplemental coefficients and release IDs were not extracted.",
        "brain6_input": meta,
        "comparison": {
            "rows": len(rows),
            "reported_significant_genetic_correlations": len(rows),
            "exact_gwas_release_matches": 0,
            "participant_overlap_verified": False,
            "pairwise_coefficients_compared": False,
            "locus_overlap_compared": False,
            "independent_replication": False,
            "novelty_classes_changed": False,
            "statistical_families_changed": False,
        },
        "script_sha256": sha256(script),
        "output_path": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(output_bytes).hexdigest(),
        "limitations": [
            "The article reports insomnia Neff=314,149; Brain6 case/control counts imply Neff=313,750. The close values do not prove source or participant identity.",
            "The abstract does not supply pair-specific rg, SE, P, or GWAS release identity; those values are not imputed.",
            "The reported shared-locus counts are reproduced as source summaries only and were not compared with Brain6 loci.",
            "The existing 2026-08-28 35-pair crosswalk remains the novelty classification source and still needs an exhaustive source-level refresh before submission.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    meta = input_metadata()
    rows = make_rows(meta)
    data = tsv_bytes(rows)
    prov = expected_provenance(rows, data, meta, Path(__file__).resolve())
    prov_data = (json.dumps(prov, indent=2, sort_keys=True) + "\n").encode()
    if args.write:
        if OUTPUT.exists() or PROVENANCE.exists():
            raise FileExistsError("Refusing to overwrite Xue et al. context outputs")
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(data)
        PROVENANCE.write_bytes(prov_data)
    else:
        if OUTPUT.read_bytes() != data:
            raise ValueError("Context TSV does not rebuild byte-identically")
        actual = json.loads(PROVENANCE.read_text(encoding="utf-8"))
        if actual != prov:
            raise ValueError("Context provenance differs from current input/source/code identities")
    print(json.dumps({"status": prov["status"], "rows": len(rows), "output_sha256": prov["output_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
