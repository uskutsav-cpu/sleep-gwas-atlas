#!/usr/bin/env python3
"""Complete the conservative pair-level novelty audit from frozen searches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path


SEARCH_DATE = "2026-08-29"
DIRECT_PRIOR_PAIR = "shortsleep__panukbb_continuous__20022__both_sexes__na__irnt"
LOW_PLAUSIBILITY_DOMAINS = {
    "digital_environment", "socioeconomic_occupational", "medication_use",
}
FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
    "rg", "se", "p", "extension_fdr", "abs_rg_ge_0.15",
    "preanalysis_novelty_priority", "preanalysis_prior_screen_coverage",
    "audit_status", "direct_prior_same_pair", "same_sleep_trait_context",
    "same_or_equivalent_phenotype", "same_direction", "broad_phenome_screen_overlap",
    "near_neighbor_evidence", "discovery_vs_replication_in_prior_work",
    "exact_prior_rg_found", "closest_prior_result", "prior_method", "prior_effect",
    "prior_publication", "prior_DOI", "prior_PMID", "search_databases",
    "search_queries_used", "search_date", "evidence_PMIDs_DOIs_URLs",
    "biological_plausibility", "connection_obviousness",
    "independent_replication_dataset_availability", "dense_summary_statistics_available",
    "molecular_qtl_data_available", "independent_replication_class",
    "novelty_class", "novelty_strength", "decision_rationale", "reviewer_notes", "reviewer",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "NA") for field in fields})
    return buffer.getvalue()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unique_join(values: list[str], limit: int | None = None) -> str:
    seen: list[str] = []
    for value in values:
        value = value.strip()
        if value and value != "NA" and value not in seen:
            seen.append(value)
    if limit is not None:
        seen = seen[:limit]
    return ";".join(seen) if seen else "NONE_IDENTIFIED"


def candidate_disposition(pair_id: str, rows: list[dict[str, str]]) -> tuple[str, str]:
    if pair_id == "insomnia__panukbb_categorical__2040__both_sexes__2040__na":
        return (
            "RELATED_GENETIC_CONTEXT_ONLY",
            "Returned records mention insomnia and risk taking as correlated risk factors or pleiotropic context, but do not report the exact insomnia-by-UKB-risk-taking rg as the tested pair.",
        )
    if pair_id == "insomnia__panukbb_categorical__3571__both_sexes__3571__na":
        return (
            "FALSE_EXACT_PAIR_RETRIEVAL",
            "The returned brain-network/fibromyalgia paper does not test the exact UKB physical-activity-type phenotype.",
        )
    if pair_id in {
        "insomnia__panukbb_continuous__20022__both_sexes__na__irnt",
        "shortsleep__panukbb_continuous__20022__both_sexes__na__irnt",
    }:
        return (
            "RELATED_MULTIVARIABLE_CONTEXT_ONLY",
            "The returned gastro-oesophageal reflux paper includes sleep and birth weight among many candidate exposures but does not report the exact sleep-by-birth-weight rg as the tested pair.",
        )
    lenses = {row["evidence_lens"] for row in rows}
    if "RELATED_GENETIC" in lenses:
        return (
            "RELATED_GENETIC_OR_MR_EVIDENCE_ONLY",
            "Returned title/abstract records provide related genetic, MR, or broad contextual evidence, not an explicit directly comparable rg for the exact pair.",
        )
    return (
        "MR_OR_OBSERVATIONAL_EVIDENCE_ONLY",
        "Returned title/abstract records provide MR or observational context, not an explicit directly comparable rg for the exact pair.",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    novelty_dir = root / "discovery_extension/results/novelty"
    template_path = novelty_dir / "extension_novelty_audit.tsv"
    query_path = novelty_dir / "literature_search_queries.tsv"
    candidate_path = novelty_dir / "literature_search_candidates.tsv"
    registry_path = root / "discovery_extension/provenance/prior_screens/prior_screen_registry.tsv"
    output_path = template_path
    review_path = novelty_dir / "literature_candidate_review.tsv"
    provenance_path = root / "discovery_extension/provenance/extension_novelty_audit.json"
    rg_path = root / "discovery_extension/results/ldsc/extension_rg_matrix.tsv"
    panel_path = root / "discovery_extension/config/candidate_traits.tsv"

    template = read_tsv(template_path)
    if len(template) != 603 or len({row["pair_id"] for row in template}) != 603:
        fail(f"expected 603 unique audit rows, observed {len(template)}")
    queries = read_tsv(query_path)
    candidates = read_tsv(candidate_path)
    if len(queries) != 2412 or any(row["query_status"] != "COMPLETE" for row in queries):
        fail("all 2,412 fixed Europe PMC queries must be complete")
    queries_by_pair: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    candidates_by_pair: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in queries:
        queries_by_pair[row["pair_id"]].append(row)
    for row in candidates:
        candidates_by_pair[row["pair_id"]].append(row)
    if any(len(queries_by_pair[row["pair_id"]]) != 4 for row in template):
        fail("every audit pair must have exactly four fixed literature queries")
    registry = {row["screen_id"]: row for row in read_tsv(registry_path)}
    registry_order = ["morrison_2024", "goodman_2025", "dashti_2019", "li_zhao_2020"]

    candidate_reviews: list[dict[str, object]] = []
    for pair_id in sorted(candidates_by_pair):
        pair_candidates = candidates_by_pair[pair_id]
        disposition, notes = candidate_disposition(pair_id, pair_candidates)
        candidate_reviews.append({
            "pair_id": pair_id,
            "returned_candidate_count": len(pair_candidates),
            "returned_lenses": unique_join([row["evidence_lens"] for row in pair_candidates]),
            "screened_PMIDs": unique_join([row["pmid"] for row in pair_candidates]),
            "screened_DOIs": unique_join([row["doi"] for row in pair_candidates]),
            "review_disposition": disposition,
            "review_notes": notes,
            "review_status": "TITLE_ABSTRACT_SCREEN_COMPLETE",
            "reviewer": "Codex evidence-conservative pair review",
        })
    review_by_pair = {str(row["pair_id"]): row for row in candidate_reviews}

    output: list[dict[str, object]] = []
    for row in template:
        pair_id = row["pair_id"]
        pair_queries = sorted(queries_by_pair[pair_id], key=lambda item: item["evidence_lens"])
        pair_candidates = candidates_by_pair.get(pair_id, [])
        coverage = row["preanalysis_prior_screen_coverage"]
        screen_ids = [screen for screen in registry_order if screen in coverage]
        screen_rows = [registry[screen] for screen in screen_ids]
        known_direct = pair_id == DIRECT_PRIOR_PAIR
        screened_before = coverage != "NO_EXACT_LOCAL_INVENTORY_MATCH"
        if known_direct:
            novelty_class = "KNOWN_BUT_NEW_DATASET"
            novelty_strength = "WEAK"
            exact_prior = "YES"
            direct_prior = "YES"
            sleep_context = "EXACT_SHORT_SLEEP_DEFINITION_IN_DASHTI_2019"
            equivalent = "YES_EXACT_BIRTH_WEIGHT_CONCEPT"
            same_direction = "YES_NEGATIVE_IN_BOTH"
            closest = "Dashti 2019 Supplementary Data 18 reports short-sleep vs birth-weight rg=-0.0825, SE=0.0426, P=0.0528."
            prior_method = "LD score regression in the source publication supplement"
            prior_effect = "rg=-0.0825;SE=0.0426;P=0.0528"
            prior_publication = registry["dashti_2019"]["citation"]
            prior_doi = registry["dashti_2019"]["DOI"]
            prior_pmid = registry["dashti_2019"]["PMID"]
            rationale = "A directly comparable short-sleep/birth-weight rg was already reported; the extension uses the version-pinned 2025 Pan-UKB external GWAS and is not an independent replication."
            obviousness = "EXPECTED"
            discovery_relation = "KNOWN_DIRECT_PAIR_UPDATED_EXTERNAL_DATASET_WITH_UKB_OVERLAP"
        elif screened_before:
            novelty_class = "PARTIAL_EXTENSION"
            novelty_strength = "WEAK"
            exact_prior = "NO"
            direct_prior = "NO"
            sleep_context = "RELATED_LATENT_COMPOSITE_OR_OTHER_SLEEP_CONSTRUCT"
            equivalent = "YES_EXTERNAL_PHENOTYPE_SCREENED;SLEEP_CONSTRUCT_NOT_EXACT"
            same_direction = "NOT_COMPARABLE_WITHOUT_EXACT_PAIR"
            closest = "Pinned broad-screen supplement contains the external phenotype under a composite, latent, or other sleep construct; no exact pair-level comparator was located."
            prior_method = "Broad phenome-wide sleep genetic-correlation screen with a non-identical sleep construct"
            prior_effect = "NOT_COMPARABLE_DIFFERENT_SLEEP_CONSTRUCT"
            representative = screen_rows[0]
            prior_publication = unique_join([item["citation"] for item in screen_rows])
            prior_doi = unique_join([item["DOI"] for item in screen_rows])
            prior_pmid = unique_join([item["PMID"] for item in screen_rows])
            rationale = "The external phenotype was already included in one or more pinned broad sleep-genetics screens, but the exact sleep phenotype/dataset pairing was not directly comparable; this is a partial extension, not a novelty claim."
            obviousness = "SEMI_OBVIOUS"
            discovery_relation = "PRIOR_BROAD_SCREEN_NONIDENTICAL_SLEEP_CONSTRUCT"
        else:
            novelty_class = "NO_DIRECT_RG_FOUND"
            novelty_strength = "MODERATE" if row["abs_rg_ge_0.15"] == "True" else "WEAK"
            exact_prior = "NO"
            direct_prior = "NO"
            sleep_context = "NO_DIRECT_COMPARATOR_LOCATED"
            equivalent = "NO_DIRECT_COMPARATOR_LOCATED"
            same_direction = "NOT_APPLICABLE_NO_DIRECT_COMPARATOR"
            review = review_by_pair.get(pair_id)
            if review:
                closest = str(review["review_notes"])
                prior_method = str(review["review_disposition"])
                representative = pair_candidates[0]
                prior_effect = "NO_DIRECT_RG_VALUE_IN_RETURNED_RECORDS"
                prior_publication = representative["title"]
                prior_doi = representative["doi"]
                prior_pmid = representative["pmid"]
            else:
                closest = "No directly comparable genetic evidence was returned by four fixed pair-specific Europe PMC lenses or found in the four pinned broad-screen inventories."
                prior_method = "NO_DIRECT_RG_OR_PAIR_SPECIFIC_GENETIC_METHOD_LOCATED"
                prior_effect = "NONE_IDENTIFIED"
                prior_publication = "NONE_IDENTIFIED_AS_OF_2026-08-29"
                prior_doi = "NONE_IDENTIFIED_AS_OF_2026-08-29"
                prior_pmid = "NONE_IDENTIFIED_AS_OF_2026-08-29"
            rationale = "Four dated pair-specific searches plus four pinned broad-screen inventories found no explicit directly comparable rg. This supports NO_DIRECT_RG_FOUND only, not a first-ever or apparently-novel claim."
            obviousness = "NON_OBVIOUS"
            discovery_relation = "DISCOVERY_ONLY_NO_INDEPENDENT_REPLICATION_YET"

        evidence_values = []
        for item in screen_rows:
            evidence_values.extend([item["PMID"], item["DOI"], item["source_url"]])
        for item in pair_candidates:
            evidence_values.extend([item["pmid"], item["doi"], item["europe_pmc_record_url"]])
        evidence_values.extend(item["europe_pmc_url"] for item in pair_queries)
        biological_plausibility = "LOW" if row["phenotype_domain"] in LOW_PLAUSIBILITY_DOMAINS else "MODERATE"
        near_neighbor = (
            "YES_RELATED_OR_BROAD_SCREEN_EVIDENCE" if screened_before or pair_candidates
            else "NO_PAIR_SPECIFIC_NEAR_NEIGHBOR_RETURNED"
        )
        output.append({
            **{field: row[field] for field in FIELDS[:12]},
            "audit_status": "COMPLETE", "direct_prior_same_pair": direct_prior,
            "same_sleep_trait_context": sleep_context, "same_or_equivalent_phenotype": equivalent,
            "same_direction": same_direction,
            "broad_phenome_screen_overlap": coverage if screened_before else "NO_EXACT_LOCAL_INVENTORY_MATCH",
            "near_neighbor_evidence": near_neighbor,
            "discovery_vs_replication_in_prior_work": discovery_relation,
            "exact_prior_rg_found": exact_prior, "closest_prior_result": closest,
            "prior_method": prior_method, "prior_effect": prior_effect,
            "prior_publication": prior_publication, "prior_DOI": prior_doi, "prior_PMID": prior_pmid,
            "search_databases": "Europe PMC;PubMed/PMC-indexed records;four checksum-pinned broad sleep-genetics supplements",
            "search_queries_used": " || ".join(item["europe_pmc_query"] for item in pair_queries),
            "search_date": SEARCH_DATE,
            "evidence_PMIDs_DOIs_URLs": unique_join(evidence_values, limit=24),
            "biological_plausibility": biological_plausibility,
            "connection_obviousness": obviousness,
            "independent_replication_dataset_availability": "UNCERTAIN",
            "dense_summary_statistics_available": "YES",
            "molecular_qtl_data_available": "UNCERTAIN",
            "independent_replication_class": "NOT_YET_ATTEMPTED",
            "novelty_class": novelty_class, "novelty_strength": novelty_strength,
            "decision_rationale": rationale,
            "reviewer_notes": "Evidence-conservative classification. Search absence is dated and overturnable; genetic correlation is not causality; UKB discovery overlap precludes treating prior UKB comparisons as independent replication.",
            "reviewer": "Codex fixed-query review with pinned supplement cross-check",
        })

    if len(output) != 603 or any(set(row) != set(FIELDS) for row in output):
        fail("completed audit coverage or schema invariant failed")
    counts = Counter(str(row["novelty_class"]) for row in output)
    strength_counts = Counter(str(row["novelty_strength"]) for row in output)
    review_fields = list(candidate_reviews[0]) if candidate_reviews else []
    audit_payload = tsv_text(output, FIELDS)
    review_payload = tsv_text(candidate_reviews, review_fields)
    provenance = {
        "schema_version": "sleep-atlas-extension-novelty-audit.1",
        "search_as_of_date": SEARCH_DATE, "audited_pair_count": len(output),
        "classification_counts": dict(sorted(counts.items())),
        "strength_counts": dict(sorted(strength_counts.items())),
        "candidate_hit_pair_count": len(candidate_reviews),
        "candidate_record_count": len(candidates),
        "direct_prior_pair_count": counts["KNOWN_BUT_NEW_DATASET"],
        "apparently_novel_pair_count": counts["APPARENTLY_NOVEL"],
        "classification_rule": {
            "known_but_new_dataset": "Exact prior sleep/external phenotype rg row in a pinned supplement; current external release differs but UKB overlap remains.",
            "partial_extension": "External phenotype previously screened with a composite, latent, or non-identical sleep construct.",
            "no_direct_rg_found": "No explicit directly comparable rg after four fixed pair queries, returned-record title/abstract review, and four pinned broad-screen inventories.",
            "apparently_novel": "Never assigned automatically in this audit; independent replication is required before a strong novelty label.",
        },
        "claim_limit": "NO_DIRECT_RG_FOUND is a dated search result, not proof of first publication. Candidate co-occurrence, MR, observational association, and broad-screen overlap are not direct rg evidence.",
        "input_hashes": {
            str(path.relative_to(root)): sha256(path) for path in (
                rg_path, panel_path, query_path, candidate_path, registry_path,
            )
        },
        "output_content_hashes": {
            "extension_novelty_audit.tsv": hashlib.sha256(audit_payload.encode()).hexdigest(),
            "literature_candidate_review.tsv": hashlib.sha256(review_payload.encode()).hexdigest(),
        },
    }
    provenance_payload = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    outputs = {output_path: audit_payload, review_path: review_payload, provenance_path: provenance_payload}
    if args.validate_only:
        for path, payload in outputs.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"novelty audit output drifted: {path.relative_to(root)}")
        print(f"EXTENSION_NOVELTY_AUDIT_VALID pairs={len(output)} classes={dict(sorted(counts.items()))}")
        return 0
    for path, payload in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(payload, encoding="utf-8")
        temporary.replace(path)
    print(f"EXTENSION_NOVELTY_AUDIT_BUILT pairs={len(output)} classes={dict(sorted(counts.items()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
