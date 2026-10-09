# Core processed-input native reproduction

Completed2026-10-09. Status: **COMPLETE_WITH_QUALIFICATIONS** for the processed-input native chain. Raw-to-harmonized-to-munged closure and scientific source/QC clearance are separate gates.

| Check | Verified result |
|---|---|
| Unchanged original commands |57, comprising45 standalone h² commands and12 correlation batches |
| Standalone h² estimates |45/45 `PRINTED_PRECISION_CONCORDANT` |
| Correlations |396/396 `PRINTED_PRECISION_CONCORDANT` |
| Original correction family |396, retained separately from the extension1200 and validation217 |
| BH q<0.05 |161 all rows:153 historical primary and8 sensitivity rows |
| Changed FDR decisions |0 |
| Independent delete-block arithmetic |45 h² and396 correlation estimates pass unchanged relative1e-12/absolute1e-15 tolerances; absolute P tolerance1e-300 |
| Historical unrounded identity |Unavailable; exact original full-precision reproduction is not claimed |

The original pinned LDSC commit is `6c673952cee74bd5c57aef1555a03b1c015399a0`, with Python3.9.23, NumPy1.21.5, pandas1.3.3 and SciPy1.7.3. The original runner SHA256 is `fa5ca8528caa1037cca55d5234e7e7721582abce479c95949bb81255e3428fa9`. Only its support/output PACKAGE is relocated to the new SSD namespace; original job identities, estimator code, inputs, options and safety floors remain fixed.

Every execution receipt verifies input and dependency hashes before/after, exact commands and requested scientific cardinality, and all produced-file hashes. The independent reviewer checked all57 command bindings,3960 literal summary fields,396 scalar P displays,358 h² NumPy scalar displays and each estimate's200 deletion values. Default blocks remain estimator-specific; this within-fit arithmetic does not establish cross-fit genomic alignment.

## Preserved precision adjudication

The initial root collator flagged the snoring–BMI tail probability while its point, SE, Z and component checks passed. Its initial implementation, receipt and failure are retained. Native P is1.02922734832538e-90. Reconstructing SE through algebraically equivalent pseudovalue variance amplified binary64 roundoff to a P relative error1.222817827813989e-12, just beyond the frozen1e-12 tolerance. Independent centered deletion variance gives P relative error7.650819321088045e-13; erfc of captured Z gives2.1181670324164024e-15. Both pass the unchanged tolerance. No estimator result, threshold, family or tolerance was changed. The separate adjudication is `reviews/independent_core_numerical_adjudication_v4.json`.

## Scientific QC remains qualified

| Outcome | Pairwise outcome intercept range | Historical all396 positives | Interpretation |
|---|---:|---:|---|
| LDL |1.221708709843142–1.2466023795867969 |4/12 |All12 exceed1.2; reproducing the failure does not clear it |
| HDL |1.4075843261724295–1.423266610795404 |8/12 |All12 exceed1.2 |
| Triglycerides |1.2045639965655714–1.2137836109143527 |8/12 |All12 exceed1.2 |

Across396 pairs,333 pass both pairwise h² Z≥4 and intercept≤1.2,48 pass Z while at least one intercept fails, and15 pass intercept while at least one Z fails. These diagnostic counts neither replace the historical primary/sensitivity classifications nor certify clinical or liability-scale validity.

MS has fitted observed total0.34811819110627373 (SE0.06958321424018578), with original P=0.005801 and K=0.002 producing reported liability h²5.982278329852824. Melanoma has fitted total0.23245358010751196 (SE0.06514611822503766), original P=0.010316 and K=0.02 producing3.730992444540936. Those supplied effective-N/sample-fraction conventions have not been certified as population variance fractions. Melanoma standalone Z=3.568187736137023 remains below4. Separately frozen total-N and balanced-ascertainment sensitivities remain to be executed; a physical-looking transformed point alone cannot clear the source.

## Evidence and resources

`tables/core_native_full_precision_rg_v4.tsv` and `tables/core_native_full_precision_h2_v4.tsv` contain the validated new precision, row-specific original log identities and preserved historical values. `statistical_validation/core_native_block_arithmetic_v4.tsv` preserves the initial flag; the independent adjudication and `reviews/independent_whole_core_receipt_v4_2.json` resolve it without overwriting. `tables/core_native_job_verification_v4.tsv` supplies the57 job checks.

The core stage took2265.251928 seconds. Maximum observed aggregate owned-worker RSS was1,094,139,904 bytes; sampling does not prove an instantaneous peak bound. Final native outputs were17,760,648 bytes. Final free internal storage was4,878,749,696 bytes. The monitor certifies that its owned process group was empty at completion.

The monitor receipt SHA256 is `910e629d9ba2c26615573a7b96b96eabd13dbf3e54dc5b9b08160fe77426fc5a`; the frozen native execution plan SHA256 is `f555dc441e8c93c528e153d2e8689edabb29ba17e4d8b0da68402ac4e4c98e88`. The independent adjudication report SHA256 is `7fed0ba8c0d0e808a3a02122deedc82247a1ff6ef1020fee5b279ccac83258c7`. All native captures and block vectors remain in the private SSD package; compact comparison tables and receipts are retained here.
