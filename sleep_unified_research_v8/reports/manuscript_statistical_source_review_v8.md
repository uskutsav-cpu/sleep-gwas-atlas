# Independent statistical and source review of the integrated manuscript draft

Reviewed 11 October 2026. Actual document: `/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/outputs/Cell_Genomics_V8_Manuscript_Review_2026-10-11/Integrated_Sleep_Atlas_Manuscript_Review_Draft.docx`. Builder was inspected as supporting authoring evidence but not executed or edited. Paragraph identifiers below are zero-based `python-docx` paragraph indices in this reviewed version. This is an authorized independent automated specialist review, not human investigator sign-off or human peer review.

## Verdict

**Retain as a qualified integrated scientific review draft, with three bounded corrections below.** The new numerical result and its reported uncertainty agree with the previously independently reviewed actual fit. No unqualified discovery, causal, mechanistic, disease-specific, source-difference, fully independent replication, empirical usability or journal-acceptance claim was found. The source-CI closure is represented honestly, and the three newly corrected source/dependence sentences were checked in the updated actual DOCX. No estimator rerun, new download, broad historical audit or numerical-review repetition was performed.

## Actionable corrections

### 1. Distinguish historical benchmark eligibility from current held interpretation

**Location:** paragraph 21, Results, comparator coverage.

**Exact wording:** “Among eligible core and extension records, 361 had no exact, related or unresolved counterpart in any of the four snapshots: 85 core and 276 extension records.”

**Evidence and implication:** The historical V7 `qualified_match_absent_across_four_snapshots.tsv` contains 361 records, including three that now have V8 source-direction holds: `core|sleep_apnea|ibd`, `core|sleep_apnea|crohn` and `core|sleep_apnea|uc`. Thus 361 is a historical benchmark-qualified coverage count; it is not the number currently unheld for biological interpretation. The subsequent mention of 144 held annotations is useful but does not identify this particular overlap. The 361 historical total is valid and need not be rewritten or rerun.

**Bounded replacement:** “Among records eligible under the historical benchmark criteria, 361 lacked an exact, related or unresolved counterpart in the four snapshots (85 core and 276 extension). Three of these core records now have V8 source-direction holds; this historical coverage count does not imply that all 361 are cleared for biological interpretation.”

If an additional current unheld descriptive coverage count is wanted, applying only these three holds gives 358 (82 core, 276 extension). Do not substitute this for the historical benchmark count or infer discovery novelty. This is a table-mask intersection, not a new statistical calculation or literature search.

### 2. Identify full-data totals and deletion uncertainty separately

**Location:** paragraph 29, Results, registry-insomnia experiment.

**Exact wording:** “Independent high-precision arithmetic reconstructed the reported ratio, SE, Z, P and interval from all 200 actual deletion ratios.”

**Evidence and implication:** The stock publication rg is the full-data ratio `gencov/sqrt(hsq1*hsq2)` from captured native totals. Its SE is independently reconstructed from the 200 centered deletion ratios; Z, P and marginal interval follow from that point estimate and SE. The bias-corrected jackknife point is a different captured quantity, approximately 0.398262. The quoted current sentence could make a reviewer believe the reported 0.399121 point is the deletion mean or the bias-corrected point. No numeric value is wrong.

**Bounded replacement:** “Independent high-precision arithmetic reconstructed the full-data ratio from the captured totals and its SE from all 200 actual deletion ratios, then recomputed Z, P and the marginal interval.”

Reference: V8 `source_review/statistical_scientific_actual_result_review_v8.md` and its companion JSON. Do not substitute the bias-corrected estimate for the stock reported ratio.

### 3. Repair scientific evidence navigation inside the manuscript bundle

**Location:** authoring supplement `supplementary_source_data/exploratory_source_transport.tsv`, fields `source_contract_ref` and `independent_actual_result_review_ref`.

**Actual values:** `../source_review/qualified_native_input_contract.json` and `../source_review/statistical_scientific_actual_result_review_v8.md`.

**Evidence and implication:** Both references resolve correctly in the original private V8 package, but do not resolve in the copied manuscript review bundle. The builder copied the evidence into `scientific_review_evidence/` while retaining the original TSV references. This is a bounded reproducibility/reviewer-access defect: the numeric candidate is not corrupted, but its machine-readable links to the source assumptions and independent result review are broken for the reviewer receiving this bundle.

**Bounded fix:** Either provide those two files at the bundle's existing expected `source_review/` relative paths, or change only the manuscript-copy references to `../scientific_review_evidence/...` and update its copied manifest accordingly. Preserve the original private V8 candidate and historical tables. No new reader/controller or statistical rerun is needed. Verify only these two relative targets in the resulting bundle.

## Corrections already confirmed in the updated DOCX

- Paragraph 31 now says “In the bounded inspected MVP rows, supplied SE was missing.” This correctly limits the observation instead of claiming a new whole-source missingness audit.
- The same paragraph now assigns the missing-INFO HapMap3 fallback explicitly to FinnGen, while MVP uses the preserved OR/CI/effective-N route. This prevents confusing MVP's supplied r2 filtering with FinnGen's no-INFO route.
- The dependence wording was corrected to attribute the identified statistical dependence to reuse of the same GWAS summary data; exact cross-enrollment/relatives remain unverified. Human overlap was not established by the source review.

## Numerical, family and source conclusions retained

The manuscript's new result, 815,347 allele-valid SNPs, rg 0.3991207522, SE 0.0506959384, nominal CI 0.2997585388–0.4984829656, P 3.46694×10⁻¹⁵, 200 finite positive deletion ratios spanning 0.390610–0.410274 and 4.8124% largest centered-variance contribution, agrees with the already reviewed actual evidence at quoted precision. This review did not rerun that arithmetic. The three displayed source estimates and intervals are rounded appropriately and distinguish the shared MVP outcome.

The 396 and 1,200 BH families remain separate. The 217 external-candidate classifications match the candidate's 23 threshold-positive, 18 additional directionally concordant, 17 QC-ineligible and 159 unavailable records; only 41 were estimated. Their outcome-side scope, not independent two-trait replication, is explicit. The one FinnGen experiment is outside those historical families and remains selected exploratory inference. The 36 source-held, 63 additional QC-failing and 297 displayed core counts agree with the current interpretation candidate, and the figure/method text says all 297 are plotted regardless of BH significance. The manuscript does not turn missing/ineligible/held categories into biological nulls or recalculate a reduced BH family.

Table 2 matches the V7 summary's four snapshot counts and union: 1,637 historical records × four gives 6,548 annotations; union related 1,183, absent 366, unresolved 88, exact zero. The manuscript expressly rejects exact-zero equivalence as a uniqueness/discovery claim, and does not claim a globally comprehensive literature benchmark. Historical positive counts are described as classifications rather than currently cleared discoveries. Processed-input reproduction and incomplete raw replay are properly separated; no uniform end-to-end source closure is claimed.

The bowel-source material does not claim original-author/HGA error, a proven all-row reversal or reinstated direction. The current signed-field hold remains mandatory, while retained P/q/classes are historical metadata. The raw uniform convention is unresolved, and a finite set of controls is not treated as a representative genomic error sample.

The new source paragraph accurately states that GCST90479148's documented Boolean-CI defect concerned a harmonised export and the consumed source was original/unharmonised. It does not treat the issue fix as certification of original 95% log-Wald intervals. The selected exploratory inference remains conditional on CI construction, constant effective N, INFO/N qualifications, approximate Finnish/EUR LD and unknown cross-enrollment. A source-specific contradictory interval definition would require withholding inferential use. P-versus-Wald equality is not used as an invented validity gate. The sampling interval is not presented as capturing all source/model uncertainty.

## Availability and claim boundary

The draft distinguishes the pushed authored-code reporting checkpoint from archival data deposition. It does not claim repository-wide CI is green, raw data are publicly redistributed, source rights are cleared, investigators have signed off or pilot responses exist. The literature/source prior art, unavailable molecular inference and absent human usability evidence are disclosed. The code-push statement is supported by the existing V8 handoff evidence; this review did not perform a new network/remote audit.

The paper is a scientific review draft requiring investigator decisions about central contribution, availability, authorship/ethics and the planned human feasibility study. Automated independent reviews do not replace those decisions. The three actionable fixes above improve scientific precision and reviewer evidence access without expanding the execution infrastructure or changing the completed science.

## Resolved addendum, 11 October 2026

The three identified corrections are resolved in the actual updated DOCX/bundle. A narrow follow-up checked only the changed benchmark/arithmetic wording and the two evidence targets: the text now labels 361 as historically benchmark-eligible, identifies three currently held records and gives the separate unheld count 358 (82 core/276 extension); it identifies the full-data ratio from captured totals and SE from all 200 deletion ratios before recomputing Z/P/interval; and both original relative evidence references resolve through the bundle's new `source_review/` files. The methods now express the log-interval transformation explicitly, and a six-query interpretation table has been added as authored demonstrations rather than measured usability.

**Final specialist verdict: retain the corrected qualified integrated scientific review draft; no unresolved numerical/source-claim correction from this review.** Human scientific approval, usable human-study evidence, source-specific public release clearance and submission readiness remain outside this automated review and are not established. No numerical rerun, figure/layout audit, historical recalculation or broad bundle audit was performed for closure.
