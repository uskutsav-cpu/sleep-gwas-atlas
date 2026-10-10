# Extension native reproduction — v4

Status: COMPLETE_VERIFIED for processed-input native reproduction. Raw-to-estimator preprocessing replay is PARTIAL and stopped with a preserved resource-guard failure; no estimator ran in that replay.

All 112 original jobs completed: 100 heritability jobs and 12 batches covering all 1,200 correlations. The completed continuation reused 107 independently verified successful receipts and ran only five unfinished original batches. The prior internal-space interruption and all preserved artifacts remain recorded as a failed attempt.

All 100 heritability estimates and 1,200 correlations agree with the historical stock logs at their printed precision. Historical full-precision estimates are unavailable, so this agreement is classified PRINTED_PRECISION_CONCORDANT. New validated full-precision values, SEs, confidence intervals, intercepts and both historical/native correction values are retained in the numerical tables.

Independent review verifies all 112 job identities, source/dependency/output receipts, 1,300 estimates and 3,800 finite 200-delete arrays. The complete 1,200-test BH family yields 603 positives in both native and frozen tables, with zero boundary changes. Maximum independently reconstructed P relative error is 6.91894e-13, below the unchanged 1e-12 tolerance. These checks do not establish genomic block alignment across different fits.

Evidence:

- tables/extension_native_full_precision_rg_v4.tsv
- tables/extension_native_full_precision_h2_v4.tsv
- tables/extension_native_job_verification_v4.tsv
- statistical_validation/extension_native_block_arithmetic_v4.tsv
- logs/extension_native_comparison_receipt_v4.json
- logs/extension_native_monitor_receipt_v4_3.json
- reviews/independent_whole_extension_review_v4.md
- reviews/independent_whole_extension_receipt_v4.json
- reviews/independent_whole_extension_v4.sha256

Independent nine-artifact seal: eb4051e7f26e76de71883c42ab0c4f4567f9c54749fc3ab9296e2ab23accaa22.

All 100 exact original extension raw bodies are recovered under the genuine v4_8 family terminal. The preprocessing/comparison pipeline has an independent actual-plan review and root admission (v7), covering 400 commands and zero estimator fits. It completed exact source-to-template preprocessing comparisons for 79 sources. Source 80 passed its immutable-source gate, then harmonization was stopped by the internal-space guard. The root receipt is `FAILED_PRESERVED_NO_AUTOMATIC_RETRY`; the pending marker is `PENDING_NOT_ADMISSIBLE`; no terminal seal exists. Twenty-one source chains lack completed outputs, including source 80's preserved partial temporary output. The existing runner explicitly refuses a retry in this attempt namespace. No new controller or retry was created.

Source recovery and unchanged preprocessing replay are separate from completed native reproduction. No source independence, primary novelty, covariance calibration or biological mechanism follows from numerical agreement.
