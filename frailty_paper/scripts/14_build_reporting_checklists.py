#!/usr/bin/env python3
"""Build item-level PRISMA 2020, PRISMA-S, and STREGA reporting tables."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

FIELDS = [
    "framework", "item_id", "topic", "reporting_requirement", "applicability",
    "status", "planned_manuscript_section", "manuscript_page_line",
    "protocol_or_data_evidence", "remaining_action", "guideline_source",
]

PRISMA_SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/8005924/"
PRISMAS_SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8270366/"
STREGA_SOURCE = "https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000022"


def row(framework, item_id, topic, requirement, section, evidence, action,
        status="PENDING", source=PRISMA_SOURCE, applicability="APPLICABLE"):
    return {
        "framework": framework, "item_id": item_id, "topic": topic,
        "reporting_requirement": requirement, "applicability": applicability,
        "status": status, "planned_manuscript_section": section,
        "manuscript_page_line": "NOT_DRAFTED",
        "protocol_or_data_evidence": evidence, "remaining_action": action,
        "guideline_source": source,
    }


def prisma_rows():
    f, s = "PRISMA-2020", PRISMA_SOURCE
    data = [
        ("1", "Title", "Identify the report as a systematic review.", "Title", "No manuscript drafted.", "Draft title after review scope and results are final.", "PENDING"),
        ("2", "Abstract", "Use the 12-item PRISMA 2020 for Abstracts checklist; see prisma_abstract.tsv.", "Abstract", "No manuscript drafted.", "Complete the separate abstract checklist after results are final.", "PENDING"),
        ("3", "Rationale", "Explain the review rationale against existing knowledge.", "Introduction", "Frozen questions in review/protocol.md.", "Write evidence-grounded rationale after screening and synthesis.", "PARTIAL_PROTOCOL_ONLY"),
        ("4", "Objectives", "State the review objectives or questions explicitly.", "Introduction", "Primary and secondary questions in review/protocol.md.", "Align manuscript objectives with the frozen protocol.", "PARTIAL_PROTOCOL_ONLY"),
        ("5", "Eligibility criteria", "Specify inclusion/exclusion criteria and how studies are grouped for synthesis.", "Methods: Eligibility", "Eligibility and evidence strata in review/protocol.md.", "Report exact criteria and synthesis groupings; record amendments.", "PARTIAL_PROTOCOL_ONLY"),
        ("6", "Information sources", "List every database, registry, website, organisation, reference list, and other source, with last-searched dates.", "Methods: Information sources", "PubMed acquisition logs; licensed database exports absent.", "Complete all authorized source searches and dates.", "PARTIAL"),
        ("7", "Search strategy", "Provide exact strategies for every database and website, including filters and limits.", "Methods / Supplement", "Exact PubMed searches in review/pubmed/search_log.tsv; other searches absent.", "Append verbatim strategies and query-level search logs for every source.", "PARTIAL"),
        ("8", "Selection process", "Describe reviewer counts, independent screening, inclusion decisions, and any automation.", "Methods: Study selection", "Dual-review process specified in review/protocol.md; screening not started.", "Report actual reviewer workflow and disagreements after screening.", "PARTIAL_PROTOCOL_ONLY"),
        ("9", "Data collection process", "Describe extraction reviewers, independence, investigator contacts, and automation.", "Methods: Data collection", "Extraction template exists; no included studies/extraction.", "Record actual extraction and verification procedures.", "PARTIAL_PROTOCOL_ONLY"),
        ("10a", "Data items: outcomes", "Define all sought outcomes and which compatible measures, time points, and analyses were collected.", "Methods: Data items", "Outcome hierarchy and extraction fields in review/protocol.md and review/extraction_template.tsv.", "Specify selection rules and report any amendments.", "PARTIAL_PROTOCOL_ONLY"),
        ("10b", "Data items: other variables", "Define other sought variables and assumptions for missing/unclear information.", "Methods: Data items", "Extraction template lists study and genetic fields; no extraction completed.", "State variable definitions and missing-data assumptions.", "PARTIAL_PROTOCOL_ONLY"),
        ("11", "Risk-of-bias assessment", "Name the tool, reviewer count/independence, and any automation used.", "Methods: Risk of bias", "JBI item-level appraisal and dual assessment specified in review/protocol.md.", "Report actual tool versions, judgments, evidence notes, and adjudication.", "PARTIAL_PROTOCOL_ONLY"),
        ("12", "Effect measures", "Define effect measures used for each outcome synthesis or presentation.", "Methods: Synthesis", "Narrative synthesis plan; no included effects.", "Specify measures by synthesis after evidence structure is known.", "PARTIAL_PROTOCOL_ONLY"),
        ("13a", "Synthesis methods: study eligibility", "Explain how studies were assigned to each planned synthesis.", "Methods: Synthesis", "Synthesis strata frozen in review/protocol.md.", "Report operational grouping and deviations.", "PARTIAL_PROTOCOL_ONLY"),
        ("13b", "Synthesis methods: data preparation", "Describe conversions and handling of missing summary statistics.", "Methods: Synthesis", "No study results collected.", "Document actual preparation decisions or state none were needed.", "PENDING"),
        ("13c", "Synthesis methods: tabulation", "Describe tables and visual displays used for individual and synthesized results.", "Methods: Synthesis", "No results tables generated.", "Specify actual display methods after extraction.", "PENDING"),
        ("13d", "Synthesis methods: synthesis", "Describe synthesis methods and rationale; if meta-analysis, models, heterogeneity methods, and software.", "Methods: Synthesis", "Protocol separates observational and genetic evidence and avoids unjustified pooling.", "State final synthesis approach and software; justify any quantitative synthesis.", "PARTIAL_PROTOCOL_ONLY"),
        ("13e", "Synthesis methods: heterogeneity", "Describe methods used to investigate sources of heterogeneity.", "Methods: Synthesis", "Planned subgroup/sensitivity comparisons in review/protocol.md.", "Report only analyses actually completed and their prespecification.", "PARTIAL_PROTOCOL_ONLY"),
        ("13f", "Synthesis methods: sensitivity", "Describe robustness/sensitivity analyses.", "Methods: Synthesis", "Planned genetic sensitivity rules in protocol/analysis plan.", "Report actual sensitivity analyses or state none were possible.", "PARTIAL_PROTOCOL_ONLY"),
        ("14", "Reporting-bias assessment methods", "Describe methods for assessing missing-result/reporting bias in each synthesis.", "Methods: Reporting bias", "No included studies or synthesis.", "Select and report a method only where appropriate and feasible.", "PENDING"),
        ("15", "Certainty methods", "Describe methods for assessing certainty/confidence for each outcome.", "Methods: Certainty", "No certainty assessment completed.", "Specify a suitable framework and apply it transparently.", "PENDING"),
        ("16a", "Study selection results", "Report counts from identification through inclusion, preferably with a flow diagram.", "Results: Study selection", "PubMed-only acquisition counts exist; screening is zero.", "Populate counts from immutable multi-source screening logs.", "PARTIAL_PUBMED_ACQUISITION_ONLY"),
        ("16b", "Excluded studies", "Cite near-eligible excluded studies and explain reasons.", "Results: Study selection", "Full-text exclusions queue is empty.", "Populate after full-text assessment; do not infer exclusions from abstracts.", "PENDING"),
        ("17", "Study characteristics", "Cite each included study and present its characteristics.", "Results: Study characteristics", "No studies included.", "Populate from checked extraction records.", "PENDING"),
        ("18", "Risk of bias in studies", "Present risk-of-bias judgments for each included study.", "Results: Risk of bias", "No studies appraised.", "Report item-level judgments and evidence locators after appraisal.", "PENDING"),
        ("19", "Individual study results", "For each outcome/study, present group summaries where relevant and effect estimate with precision.", "Results: Individual studies", "No studies extracted.", "Populate from verified reports; distinguish unavailable from null.", "PENDING"),
        ("20a", "Synthesis results: contributing studies", "Summarize characteristics and risk of bias of contributing studies.", "Results: Synthesis", "No synthesis or contributing studies.", "Summarize for each synthesis after extraction/appraisal.", "PENDING"),
        ("20b", "Synthesis results: estimates", "Present all synthesis results and precision/heterogeneity where relevant; state direction.", "Results: Synthesis", "No quantitative review synthesis; frozen atlas analyses are separate and reported in their own audit.", "Report only planned, eligible, frozen review results.", "PENDING"),
        ("20c", "Synthesis results: heterogeneity", "Present results of investigations into heterogeneity.", "Results: Synthesis", "No such analysis completed.", "Report all planned/completed investigations.", "PENDING"),
        ("20d", "Synthesis results: sensitivity", "Present all sensitivity-analysis results.", "Results: Synthesis", "No such analysis completed.", "Report all completed sensitivity analyses.", "PENDING"),
        ("21", "Reporting-bias results", "Present missing-result/reporting-bias assessments for each synthesis assessed.", "Results: Reporting bias", "No synthesis assessed.", "Report method-specific results or explain non-applicability.", "PENDING"),
        ("22", "Certainty results", "Present certainty/confidence assessments for each outcome assessed.", "Results: Certainty", "No certainty assessment.", "Present outcome-level judgments and rationale.", "PENDING"),
        ("23a", "Discussion: interpretation", "Interpret results in the context of other evidence.", "Discussion", "No completed evidence synthesis.", "Synthesize only after screening and analysis freeze.", "PENDING"),
        ("23b", "Discussion: evidence limitations", "Discuss limitations of the included evidence.", "Discussion", "No included studies.", "Address risk of bias, inconsistency, imprecision, and indirectness as supported.", "PENDING"),
        ("23c", "Discussion: review-process limitations", "Discuss limitations of the review process.", "Discussion", "Manual databases absent; screening not started.", "Disclose search, access, screening, and synthesis limitations.", "PARTIAL_KNOWN_GAPS"),
        ("23d", "Discussion: implications", "Discuss implications for practice, policy, and future research.", "Discussion", "No evidence synthesis.", "Keep implications proportionate to evidence certainty.", "PENDING"),
        ("24a", "Registration", "Give registry and number or state that the review was not registered.", "Other: Registration", "Protocol.md states no registration recorded.", "State registration status and chronology explicitly.", "PARTIAL_PROTOCOL_ONLY"),
        ("24b", "Protocol access", "State where the protocol can be accessed or that none was prepared.", "Other: Registration", "Frozen protocol at review/protocol.md.", "Provide stable access/location in final report.", "PARTIAL_PROTOCOL_ONLY"),
        ("24c", "Protocol amendments", "Describe and explain amendments to registration/protocol.", "Other: Registration", "Protocol records no amendments at v1 freeze.", "Update after any amendments and disclose timing/reason.", "PARTIAL_PROTOCOL_ONLY"),
        ("25", "Support", "Describe financial/nonfinancial support and funder/sponsor roles.", "Other: Support", "No manuscript or finalized funding statement.", "Collect author/funder disclosures.", "PENDING"),
        ("26", "Competing interests", "Declare review-author competing interests.", "Other: Competing interests", "No author declarations compiled.", "Obtain and report declarations.", "PENDING"),
        ("27", "Data/code/material availability", "State availability and location of forms, extracted data, analysis data, code, and other materials.", "Other: Availability", "Local templates/code exist; publication repository not selected.", "Freeze and deposit shareable materials; document licensed/restricted data.", "PARTIAL"),
    ]
    return [row(f, *x[:4], x[4], x[5], status=x[6], source=s) for x in data]


def abstract_rows():
    f, s = "PRISMA-2020-ABSTRACT", PRISMA_SOURCE
    items = [
        ("1", "Title", "Identify the report as a systematic review.", "Title", "No manuscript.", "Draft after completion."),
        ("2", "Objectives", "State the review objective/question.", "Abstract", "Frozen questions in protocol.", "Summarize objective after scope freeze."),
        ("3", "Eligibility", "Summarize inclusion and exclusion criteria.", "Abstract", "Criteria in protocol.", "Summarize final criteria."),
        ("4", "Information sources", "Name sources and last-search date.", "Abstract", "PubMed snapshot only; other sources pending.", "Update after all searches."),
        ("5", "Risk of bias", "State how study risk of bias was assessed.", "Abstract", "JBI method in protocol; assessments not done.", "Report actual method and status."),
        ("6", "Synthesis methods", "State methods used to present and synthesize results.", "Abstract", "Protocol synthesis plan exists.", "Describe completed synthesis."),
        ("7", "Included studies", "Give total included studies/participants and key features.", "Abstract", "No studies included.", "Populate from final extraction."),
        ("8", "Synthesis results", "Present main outcomes and study counts; summary estimate/precision if pooled.", "Abstract", "No review synthesis.", "Report only frozen synthesis results."),
        ("9", "Limitations of evidence", "Summarize limitations such as bias, inconsistency, and imprecision.", "Abstract", "Evidence appraisal not done.", "Summarize after appraisal."),
        ("10", "Interpretation", "Give overall interpretation and important implications.", "Abstract", "No completed synthesis.", "Keep conclusions proportionate."),
        ("11", "Funding", "State primary review funding source.", "Abstract", "Funding not compiled.", "Obtain author/funder statement."),
        ("12", "Registration", "Give registry and registration number.", "Abstract", "No registration recorded in protocol.", "State not registered, if unchanged."),
    ]
    return [row(f, *x[:4], x[4], x[5], status="PENDING", source=s) for x in items]


def prismas_rows():
    f, s = "PRISMA-S", PRISMAS_SOURCE
    items = [
        ("1", "Database name", "Name each database and its search platform.", "Methods: Information sources", "PubMed snapshot; licensed exports absent.", "Report exact database/platform inventory.", "PARTIAL"),
        ("2", "Multi-database searching", "If searched together on one platform, name platform and databases.", "Methods: Information sources", "No multi-database platform recorded.", "Report as not used or describe it if added.", "PENDING"),
        ("3", "Study registries", "List registries searched.", "Methods: Information sources", "No registry search logged.", "Search and record, or explicitly state none searched.", "PENDING"),
        ("4", "Online resources and browsing", "Describe purposefully browsed online/print sources and methods.", "Methods: Other sources", "Resource audit is distinct from review retrieval.", "Log any purposeful review-search browsing.", "PARTIAL"),
        ("5", "Citation searching", "State whether cited/citing references were searched and how.", "Methods: Other sources", "Citation chaining not logged.", "Perform and log if used, otherwise state not done.", "PENDING"),
        ("6", "Contacts", "State if studies/data were sought from authors, experts, or others.", "Methods: Other sources", "No contact logged.", "Record not done or document contacts/outcomes.", "PENDING"),
        ("7", "Other methods", "Describe other information sources/search methods.", "Methods: Other sources", "Supplementary retrieval methods not finalized.", "Log all methods used.", "PARTIAL"),
        ("8", "Full search strategies", "Provide exact strategies as run for every database/information source.", "Supplement: Search strategies", "Exact PubMed queries exist; licensed strategies absent.", "Append verbatim syntax and run logs for each source.", "PARTIAL"),
        ("9", "Limits and restrictions", "State no limits or describe limits and rationale.", "Methods: Search strategy", "PubMed query documented; review-wide limits not audited.", "Document per-source date/language/study/publication restrictions.", "PARTIAL"),
        ("10", "Search filters", "State whether published filters were used and cite them.", "Methods: Search strategy", "Filter use not assessed.", "Record no filter or cite/adapted filter.", "PENDING"),
        ("11", "Prior work", "State whether prior review strategies were adapted/reused and cite sources.", "Methods: Search strategy", "Adaptation record not finalized.", "Record adaptation history.", "PENDING"),
        ("12", "Updates", "Describe search update methods.", "Methods: Search updates", "NCBI index reconciliation done; no final update search.", "Run and log final update before review close.", "PARTIAL"),
        ("13", "Dates of searches", "Give last search date for each search strategy.", "Methods: Information sources", "PubMed snapshot dated 2026-09-22; manual sources absent.", "Record final date for every source.", "PARTIAL"),
        ("14", "Peer review", "Describe search-strategy peer review.", "Methods: Search strategy", "Independent peer review not documented.", "Arrange/record review or disclose not performed.", "PENDING"),
        ("15", "Total records", "Report records from each database and other source.", "Results: Study selection", "PubMed counts reconciled; other sources absent.", "Give per-source counts and reconcile flow.", "PARTIAL"),
        ("16", "Deduplication", "Describe deduplication process and software.", "Methods: Record management", "PMID/DOI/title audit exists for PubMed snapshot.", "Run combined-source dedup and preserve audit.", "COMPLETE_FOR_PUBMED_SNAPSHOT"),
    ]
    return [row(f, *x[:4], x[4], x[5], status=x[6], source=s) for x in items]


def strega_rows():
    f, s = "STREGA-STROBE", STREGA_SOURCE
    items = [
        ("1", "Title and abstract", "Identify the genetic association design in the title/abstract and give a balanced summary.", "Title / Abstract", "No manuscript.", "Draft after analysis freeze.", "PENDING"),
        ("2", "Background/rationale", "Explain scientific background and rationale.", "Introduction", "Frozen analysis plan and review protocol.", "Integrate after review and results are final.", "PARTIAL_PROTOCOL_ONLY"),
        ("3", "Objectives", "State specific objectives and prespecified hypotheses.", "Introduction", "Analysis plan v1 is frozen pre-new-results.", "Report hierarchy and any versioned amendments.", "PARTIAL_PROTOCOL_ONLY"),
        ("4", "Study design", "Present key design elements early.", "Methods: Design", "Existing analysis plan; frailty-specific study not run.", "Describe each source and analysis design.", "PENDING"),
        ("5", "Setting", "Describe settings, locations, and relevant dates, including recruitment/exposure/outcome periods.", "Methods: Setting", "Source metadata incomplete.", "Extract source-specific dates/settings from verified studies.", "PENDING"),
        ("6", "Participants", "Give eligibility, selection, and participant ascertainment; genetic extension: report how variants/genes were selected when relevant.", "Methods: Participants / Sources", "GWAS sources partly indexed; several unresolved.", "Report source inclusion, exclusions, and variant/locus selection.", "PARTIAL"),
        ("7", "Variables", "Define outcomes, exposures, predictors, confounders, and diagnostic criteria.", "Methods: Phenotypes", "FI is primary; phenotype definitions and alternatives require source verification.", "Give operational definitions and coding/direction.", "PARTIAL"),
        ("8", "Data sources and measurement", "Describe data sources and measurement comparability; genetic extension: genotyping, imputation, QC, and any haplotype inference.", "Methods: Data sources / QC", "Existing summary data; raw FI analysis inputs absent.", "Report verifiable source-specific assay/genotyping/QC; mark unavailable details.", "PARTIAL"),
        ("9", "Bias", "Describe efforts to address potential sources of bias.", "Methods: Bias control", "Overlap matrix documents cohort-level overlap; exact intersections unknown.", "Report overlap, selection, population-structure, and other bias controls.", "PARTIAL"),
        ("10", "Study size", "Explain how study size was determined.", "Methods: Participants", "Source sample counts are incomplete/variable.", "Report source-derived sample size and cases/controls.", "PENDING"),
        ("11", "Quantitative variables", "Explain handling of quantitative variables; genetic extension: state how quantitative-trait effect estimates are modeled/interpreted.", "Methods: Statistical analysis", "No new harmonization/QC.", "Specify scaling, transformations, effect units, and thresholds.", "PENDING"),
        ("12", "Statistical methods", "Describe statistical methods, confounding control, missing data, sensitivity/subgroups; genetic additions include stratification, relatedness, HWE, and haplotype methods where applicable.", "Methods: Statistical analysis", "Analysis plan locks methods; data-specific analyses not run.", "Report all models, QC, multiplicity, ancestry, relatedness, HWE/haplotype handling as applicable.", "PARTIAL_PROTOCOL_ONLY"),
        ("13", "Participants/results flow", "Report participant numbers at each stage and reasons for nonparticipation; genetic extension: numbers with genotyping attempted and successful.", "Results: Participant flow", "Frozen summary results only; original inputs unavailable.", "Report available sample flow and genotype counts from source papers.", "PENDING"),
        ("14", "Descriptive data", "Describe participant characteristics and missingness; genetic extension: ancestry/population descriptors and genotype/sample QC summaries.", "Results: Participants", "Source-level descriptions incomplete.", "Report available characteristics and explicitly mark unavailable data.", "PENDING"),
        ("15", "Outcome data", "Report outcome events/summary data; genetic extension: outcome-specific counts and quantitative estimates.", "Results: Outcomes", "No new frailty analyses.", "Report denominators and estimates from frozen tables.", "PENDING"),
        ("16", "Main results", "Give unadjusted/adjusted estimates with precision and clarify adjustment; genetic extension: identify discovery versus replication results.", "Results: Main results", "Frozen FI global-rg table exists; source inputs cannot be fully reprocessed.", "Report only frozen estimates, multiplicity family, and provenance caveat.", "PARTIAL_FROZEN_RESULT_ONLY"),
        ("17", "Other analyses", "Report other analyses, including subgroups and interactions.", "Results: Additional analyses", "No frailty-specific sensitivity or subgroup analyses.", "List all prespecified analyses and outcomes, including null/blocked.", "PENDING"),
        ("18", "Key results", "Summarize findings against objectives.", "Discussion: Key results", "No complete frailty paper results.", "Summarize without causal or replication inflation.", "PENDING"),
        ("19", "Limitations", "Discuss limitations, direction/possible magnitude of bias.", "Discussion: Limitations", "Reviewer audit records overlap and source-reprocessing limits.", "Integrate evidence-specific limitations.", "PARTIAL_AUDIT_ONLY"),
        ("20", "Interpretation", "Give cautious overall interpretation in light of objectives, limitations, multiplicity, and other evidence.", "Discussion: Interpretation", "Reviewer audit downgrades current result pattern to exploratory.", "Use bounded, noncausal language.", "PARTIAL_AUDIT_ONLY"),
        ("21", "Generalisability", "Discuss external validity/generalizability.", "Discussion: Generalisability", "Current FI/sleep sources are predominantly European and often UK Biobank.", "Restrict inference to represented ancestry/cohorts and discuss portability.", "PARTIAL_SOURCE_AUDIT_ONLY"),
        ("22", "Funding", "Report funding and funder role for the current and relevant underlying studies.", "Other: Funding", "Funding statements not compiled.", "Report project funding and distinguish source-study funding where relevant.", "PENDING"),
    ]
    return [row(f, *x[:4], x[4], x[5], status=x[6], source=s) for x in items]


def write(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def preserve_existing_locators(rows: list[dict[str, str]], path: Path) -> list[dict[str, str]]:
    """Keep reviewed locator cells when rebuilding checklist content."""
    if not path.is_file():
        return rows
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"framework", "item_id", "manuscript_page_line"}.issubset(reader.fieldnames):
            raise ValueError(f"Existing checklist has no locator-compatible schema: {path}")
        existing = {
            (item.get("framework", ""), item.get("item_id", "")): item.get("manuscript_page_line", "")
            for item in reader
        }
    for item in rows:
        locator = existing.get((item["framework"], item["item_id"]))
        if locator:
            item["manuscript_page_line"] = locator
    return rows


def main():
    parser = argparse.ArgumentParser()
    default_out = Path(__file__).resolve().parents[1] / "paper/reporting_checklists"
    parser.add_argument("--outdir", default=str(default_out))
    args = parser.parse_args()
    out = Path(args.outdir)
    datasets = {
        "prisma.tsv": prisma_rows(),
        "prisma_abstract.tsv": abstract_rows(),
        "prisma_s.tsv": prismas_rows(),
        "strega.tsv": strega_rows(),
    }
    for name, rows in datasets.items():
        path = out / name
        write(path, preserve_existing_locators(rows, path))
        print(f"Wrote {out / name}: {len(rows)} items")


if __name__ == "__main__":
    main()
