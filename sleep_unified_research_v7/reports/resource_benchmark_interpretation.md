# V7 resource benchmark: results and interpretation

**Build:** 2026-10-10, from frozen V4 result tables and the verified comparator snapshots named below.

**Population:** 1,637 result records (396 core, 1,200 extension, 41 historical outcome-side validation); the three original statistical families remain distinct.

**Scope:** descriptive source/coverage matching only. No GWAS or LDSC estimates, calibrated contrasts, independent replications, or biological findings were generated.

## Accessible row-level comparisons

| Snapshot | A exact source/estimand | B related, nonidentical | C incompatible | D absent in snapshot scope | E insufficient evidence | Apparent candidates | B/C/E among candidates | Eligible core/extension without A/B/E in this snapshot |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Human GWAS ATLAS release 3 (v20191115) | 0 | 210 | 0 | 1,344 | 83 | 240 | 100.00% | 1,321 |
| Morrison 2024 S10–S15 | 0 | 888 | 0 | 749 | 0 | 888 | 100.00% | 713 |
| Goodman 2025 SD24 | 0 | 492 | 0 | 815 | 330 | 822 | 100.00% | 765 |
| Fan 2026 SD2/SD4/SD11 | 0 | 40 | 0 | 1,549 | 48 | 88 | 100.00% | 1,453 |
| Union of four row-level snapshots | 0 | 1,183 | 0 | 366 | 88 | 1,247 | 100.00% | **361 qualified match-absent** |

The annotation table has 6,548 unique record/comparator rows: all 1,637 V4 records against each of four accessible row-level snapshots. The Human GWAS ATLAS release contains 1,393,615 unique pair keys; Morrison S10–S15 contains 8,418 extracted rows; Goodman SD24 contains 2,274 nonmissing score cells across **six** composite sleep-health scores; Fan SD11 contains 4,046 unique sleep/outcome rows. Input and generated-table hashes are in `../manifests/cross_resource_build_manifest.json`. Every candidate result row reference is retained; reference lists are not truncated.

### Why HGA has no class-A records

The initial matcher labeled 30 HGA rows exact from publication PMID, approximate sample size, ancestry text, and phenotype labels. Independent reviewers checked the actual result pairs and found opposite-signed rg values for 16 of those 30 candidates. More fundamentally, the metadata did not establish source-release identity, phenotype ascertainment/coding, cohort frame, effect direction, estimator/reference compatibility, or the exact estimand required by the preregistered A definition. The final builder conservatively classifies all 30 as E until those fields are adjudicated. V7 exports no HGA numeric result values; the sign comparison was used only to challenge the source-match classification.

All 30 candidates also have Morrison class-B coverage, so the union's precedence rule classifies them B. This is why HGA E increases by 30 while union B increases by 30 and union E remains 88. The corrected union has 1,247 apparent candidate records, all classified B/C/E under the exact-match rule. The 100% figure measures source/estimand matching difficulty among label/code candidates; it is not a statistical failure rate.

The 361 qualified match-absent records are core/extension rows that pass their existing h²-Z/intercept flags and have the required source, release, and definition metadata, with no A/B/E class across these four row-level snapshots. They comprise 85 core and 276 extension rows and are D in each of the four comparator columns. This is a release- and snapshot-bounded coverage statement. It does not establish global absence, novelty, a new association, or greater accuracy. The per-snapshot eligible counts in the table are not the all-snapshot qualified count.

Class D means no counterpart in a named complete table or dated release under the recorded matching rules; it is not evidence of novelty. Class E means comparability is unresolved. No C classification was assigned; uncertain candidates remain E rather than being forced into incompatibility.

## Comparator access and legal scope

Fan SD2/SD4/SD11 were read from their verified SSD copies. The annotation records their workbook row references and derived classes, not Fan rg, SE, P, or FDR values. The article's CC BY-NC-ND license means that the public-release status of derivative row references and classifications still requires accountable human rights review; this build does not determine permission to redistribute them. Publisher workbooks remain outside the repository.

SleepChart is carried forward only as context from the completed V4 audit, which found partial/related duration–disease coverage. Its inspected workbook is absent from the active checkout and was not re-downloaded, so SleepChart is neither in the four-snapshot union nor a negative search. CTG-VL has no direct export in this evidence set. Morrison is represented by its six published S10–S15 tables, not by a separate CTG-VL query.

## Scientific interpretation

The project is a source-resolved marginal-association resource, not a new method or biological-discovery study. Its defensible contribution is the curated 12-construct × 133-outcome-source configuration across two frozen testing families, explicit phenotype/source metadata, family-corrected results, and distinctions among estimated, underpowered, QC-ineligible, unavailable, and outcome-side-validation records. Existing resources already provide broad sleep-genetic screens, multivariate sleep models, source links, LDSC workflows, and downstream analysis tools. This benchmark does not establish unique associations, superior accuracy, or a first-of-kind resource.

The existing result package has relevant limits: the 41 external estimates reuse discovery sleep GWAS and are not independent two-trait replication; uncertainty intervals are marginal; discovery/validation sampling covariance is not calibrated; 62 sensitivity rows are diagnostics rather than tests of between-fit differences; and questionnaire phenotypes are not clinical diagnoses. Existing use cases (insomnia–frailty with its complaint-based definition, separately corrected insomnia–abdominal-pain, outcome-side validation labels, and a sensitivity estimate with a marginal interval crossing zero) are retrieval/interpretation examples, not new V7 findings.

The independent scientific reviewer found that a bounded resource paper is plausible, but that the integrated project does not establish a new biological result, inference method, independent two-trait replication, or measured user benefit. The V7 paired AI-answer benchmark was deferred before task sampling and remains untested; the V4 automated consumer review is not a baseline comparison or human usability study.

## Authoring decision

Investigators can begin drafting one carefully qualified resource manuscript using the existing V4 result figures/tables and this V7 coverage figure. Keep core, extension, and historical validation families separate. Do not frame the work as a new method, first discovery, independent replication study, global-absence claim, calibrated measurement comparison, or human-usability demonstration.

The evidence does **not currently support a Nature Communications submission-readiness claim**: no major new biological or methodological advance was demonstrated, and public data/code permissions and human approvals remain unresolved. Nature Communications describes published work as important advances for specialists and requires transparent data availability and, for central custom code, code availability; source-specific restrictions must be resolved with the journal and disclosed rather than assumed away. See the [Guide to Authors](https://www.nature.com/ncomms/submit/guide-to-authors), [reporting and data/code availability policy](https://www.nature.com/ncomms/editorial-policies/reporting-standards), and [policy guide](https://www.nature.com/ncomms/submit/policy).

No manuscript prose was written. Human rights, metadata/overlap review, ethics, authorship, funding, conflicts, AI-use disclosure, data/code statements, and final journal selection remain with the investigators.
