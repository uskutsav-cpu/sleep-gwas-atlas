# Independent finite assessment of the MVP confidence-interval issue

Reviewed 11 October 2026. Scope: Open Targets issues 3917/3949, their relevant comments, the small affected-study list, exact already preserved acquisition/metadata route, and the previously bounded source-schema observation. No GWAS body download, fit, preprocessing replacement, historical edit or broad audit was performed.

## Decision

**Retain the completed estimate only as the already qualified exploratory source/construct-transport result. Do not upgrade confidence-interval semantics to verified.** The new evidence establishes a relevant harmonised-export history for the exact accession, but does not demonstrate that the consumed original unharmonised body has that defect. It therefore strengthens provenance and uncertainty disclosure without establishing a reason to retract this conditional estimate or change its Z route after seeing the result. Confirmatory or source-schema-certified inference remains unsupported.

Material adverse evidence must be acted on if it concerns the exact consumed source and uncertainty construction. In particular, an original-export CI definition incompatible with the frozen 95% log-Wald formula, evidence of mixed/invalid original CIs, or a source warning against using the original uncertainty would require a new hold on inferential use; the old fit would remain preserved. None of those exact-source conclusions follows from the issue evidence inspected here.

## What the issue evidence actually says

Issue 3917 reports Boolean `True` replacing numeric CI endpoints in GWAS Catalog **harmonised** MVP output, illustrating GCST90480017. Its original unharmonised examples retain numeric endpoints. June 26 comments report 2,414 of 6,022 MVP studies affected and re-harmonisation after a stale processing-code problem; June 30 confirms propagation of the example fix. These are direct downstream/provider technical reports, not original publication statements defining confidence-interval coverage or construction. [Issue 3917 and comments](https://github.com/opentargets/issues/issues/3917).

Issue 3949 reports missing supplied SE, failure of a CI-based rescue against Open Targets' PZ QC, and the harmonised-CI defect. July 16 comments show corrected numeric CIs while SE remains missing, with beta/P-derived SE chosen. Later comments report successful re-harmonisation for most studies. The project used an operational PZ-QC criterion; it did not establish that every original MVP confidence interval was scientifically invalid. Numeric-CI restoration is distinct from validating interval semantics. [Issue 3949](https://github.com/opentargets/issues/issues/3949), [July 16 technical update](https://github.com/opentargets/issues/issues/3949#issuecomment-3079585398), [July 18 re-harmonisation update](https://github.com/opentargets/issues/issues/3949#issuecomment-3088969188).

The linked affected-study attachment contains 2,414 rows. **GCST90479148 is present**, with path:

`/nfs/ftp/public/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479148/harmonised/GCST90479148.h.tsv.gz`.

This is accession-specific evidence that its harmonised export was affected; it is not evidence that the sibling original file was affected. The 325,903-byte primary list was saved as `source_review/OpenTargets_CI_affected_studies.txt`. [Primary affected-study attachment](https://github.com/user-attachments/files/21101306/ci_affected_studies.txt).

## Relationship to our actual acquisition route

The existing acquisition's `member.url` identifies:

`https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479148/GCST90479148.tsv.gz`.

Its file is the original `GCST90479148.tsv.gz`, not `harmonised/GCST90479148.h.tsv.gz`. Exact source metadata states `is_harmonised: false`, GWAS-SSF v1.0, GRCh38/1-based, EUR N=388315 and MD5 `0e9231f2bc828189805436f9aa4fd3c8`, matching the already verified original source. [Exact original source metadata](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479148/GCST90479148.tsv.gz-meta.yaml).

Preserved provenance inspected read-only:

- `.../sleep_unified_research_v4/validation_raw_replay_v3/receipts/gwas_catalog_GCST90479148.json`, containing the exact original URL.
- `.../sleep_unified_research_v4/validation_pipeline_replay_v4/receipts_v4/gwas_catalog_GCST90479148.source.json`, identifying the original body's path/size/MD5/SHA.
- V8 `source_review/GCST90479148_meta_current.yaml`, the previously retrieved exact original metadata.
- V8 `source_review/MVP_bounded_CI_schema_diagnostic.json`, the existing first-1,000 eligible-row diagnostic with numeric, nearly log-symmetric original CIs and missing supplied SE.

Here `...` means `/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/`. No original body was reread or rehashed for this follow-up. Route identity matters scientifically because interpreting a defect in a transformation output as proof about the preceding source would be an unsupported attribution.

## Statistical implications and unresolved semantic contract

The existing diagnostic found 27/1,000 bounded eligible raw rows with >10% relative disagreement between supplied P and normal P from CI-derived Z. This is not a representative genome-wide error rate or proof of corrupted CIs. Score/SPA versus effect/Wald uncertainty and rounding can produce differences; exact attribution remains unestablished. The new issue reports do not settle those possibilities for this original accession.

Conversely, neither numeric endpoints nor near log symmetry certifies 95% Wald coverage. The consumed transformation assumes `SE_logOR=(log(ci_upper)−log(ci_lower))/(2×1.959963984540054)`. A file may be numerically well formed yet encode another interval construction. Its exact source-method contract remains necessary to call this verified effect uncertainty.

A beta/P-derived SE would produce matching normal P largely by construction. Open Targets' successful PZ QC after that route is therefore not an independent validation of the original CI route, and is not a basis to replace the frozen input after viewing a favorable correlation. It also does not independently validate P calibration, ascertainment or LDSC suitability. No P-versus-Wald equality gate is invented in this review.

The exact original accession metadata inspected here does **not** state CI confidence level or formula. The inspected GWAS-SSF specification defines numeric upper/lower CI fields but does not establish 95% coverage or log-Wald construction. Previously inspected primary MVP publication methods identify the association framework, but no accession-specific original-export CI contract was established in that review. This is a bounded statement about inspected evidence, not proof that no such documentation exists; the separate original-publication/supplement investigation may supply it. [GWAS-SSF field specification](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.0.2.tex), [MVP primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/).

## Claim and documentation consequences

The new result remains a conditional, selected exploratory estimate using a new FinnGen sleep source and reused MVP outcome. Its already quoted nominal CI/P must remain labelled conditional on CI95/log-Wald, effective-N, source-QC and LD assumptions. This technical evidence neither creates new biological support nor removes the prior limitation. It specifically prevents a statement that MVP CI semantics are resolved merely because a Catalog harmonisation bug was fixed.

Recommended additional source note: “GCST90479148 was listed among MVP studies with a historical GWAS Catalog harmonised-CI defect. This analysis used its original unharmonised numeric-CI source, not the affected harmonised export. The original CI confidence level/construction remains an explicit assumption; no P-versus-Wald equality criterion was used to certify or reject it.”

No new controller, large download, source retry, estimator rerun or historical correction is warranted by the evidence currently established. If primary documentation resolves the original CI semantics, update the source contract and corresponding interpretation prospectively. If it establishes incompatibility, withhold inferential use rather than preserving the stronger claim on the strength of clean software output.
