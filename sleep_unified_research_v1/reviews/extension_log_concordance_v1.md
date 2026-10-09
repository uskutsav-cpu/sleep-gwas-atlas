# Additive original-extension log concordance audit v1

Client audit date: 2026-10-08. This review incorporates newly recovered evidence after `statistical_validity_v1.md` and `numerical_reproduction_v1.md`; those earlier reports were not edited. Their descriptions of missing extension logs accurately describe the earlier audit state. **The original 100 extension h2 logs and 12 extension rg logs are now available and agree with the frozen 100/1,200 tables.**

## Recovery and independent method

The bounded `logs/archived_extension_recovery_receipt_v1.json` was read first. It records a completed sequential archive recovery from `/Volumes/Extreme SSD/Codex-Archive/2026-08-26.tar.gz`, with 360 retained artifacts, 227 receipt-pinned artifacts and no missing receipt-pinned artifacts. This reviewer did not reopen the large tar archive, read any GWAS source file or execute LDSC.

The newly recovered logs were read at `/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/logs/`. Their paths, archive members, sizes and hashes crosswalk to `tables/archived_extension_recovery.tsv`.

`reviews/extension_log_concordance_v1.py` uses only Python standard-library parsing, hashes, arithmetic and an independently written BH calculation. It imports no existing collator. Its JSON receipt hashes the frozen tables, bounded recovery receipt, recovery ledger and all 112 inspected logs. Its 1,300-row TSV records each table line, corresponding log path/line, SNP counts and P representations. It writes only this additive review's outputs.

All 112 rehashed logs match their recovery-ledger current hashes and sizes. Their ledger classification is `HISTORICAL_LOG_RECOVERED_CURRENT_HASH_ONLY`: none has a separately frozen expected historical log hash. The numerical agreement supports original-output concordance, but is not a claim of prior checkpoint SHA-256 identity for the log files themselves.

## Results

| Independent check | Result |
|---|---:|
| Unique extension h2 trait records/logs | 100/100 |
| H2 numerical fields checked | 995, all agree |
| Negative attenuation-ratio classifications | 5/5 agree |
| Recorded h2/intercept gates | 100/100 pass |
| Original sleep rg logs, each with 100 outcome rows | 12/12 |
| Exact Cartesian sleep–extension pairs | 1,200/1,200 |
| Authoritative serialized rg numerical fields | 12,000, all agree |
| Additional fixed-width P representations | 1,200 checked against scalar-P printing precision |
| Like-stage SNP-count fields | 3,600, all agree |
| BH-adjusted values from all original P values | 1,200/1,200 agree; 603 positives |
| Numeric mismatches | 0 |

The ten authoritative rg fields comprise rg, SE, Z, scalar P, outcome h2/SE/intercept/intercept SE and cross-trait intercept/intercept SE. The extra fixed-width P representation is retained as an eleventh printed representation. Its rounding is not substituted for the authoritative scalar P. All comparisons here are also exact as parsed binary64 values; the declared round-trip tolerance was relative 5×10⁻¹⁴ with zero absolute tolerance.

H2 checks include h2, SE, derived h2 Z, intercept/SE, input and regression SNP counts, lambda GC, mean chi-square and attenuation ratio when numeric. Five original logs instead state a negative ratio; their frozen `LT_ZERO`/`NEGATIVE_NOT_ESTIMATED` classifications are preserved, not replaced by an invented numeric value. Scale and gate labels also agree.

The fixed-width summary P differs from scalar P in 359 rows, including 341 summary zeros. These are display-precision differences; retained scalar P values reproduce every frozen P field. The summary representation agrees with scalar values at the documented four-decimal display precision. No threshold or family denominator was changed.

As an exact example, frozen `extension_trait_readiness.tsv:2` agrees with the recovered Vitamin D h2 log at line 26, including h2=0.0872, SE=0.0165 and input/regression SNP count=1,166,890. Frozen `extension_rg_matrix.tsv:2` agrees with `rg_insomnia.log:2636`: rg=−0.0216, SE=0.0253, scalar P=0.3932, extension read count=1,217,311, merged count=1,130,007 and valid-allele count=1,129,257.

All 12 first-pair sleep h2/intercept diagnostics in the original rg logs pass the recorded thresholds. These are pair-subset, observed-scale diagnostics. They do not replace the separately recovered 12 standalone sleep h2 records or establish source-level preprocessing validity.

## SNP-count stage distinction

Standalone extension h2 read counts range from 1,165,648 to 1,166,942; all 1,200 rg blocks initially report 1,217,311 external rows. The like-stage counts agree exactly with their respective frozen tables. Comparing these two read counts directly would compare different parser stages.

The independently inspected pinned `ldscore/sumstats.py:236–237,320–321` sends standalone h2 through a parser with `dropna=True`. At `:427–430`, external rg inputs use `dropna=False` and are filtered after merging. `ldscore/parse.py:99–115` implements the conditional missing-value removal. Thus identical munged files can legitimately produce these different initial log counts. The log count discrepancy does not establish a source mismatch. This reviewer did not independently recount source-file missing values; source-content checks remain a separate provenance/reproduction receipt.

## What this additive audit establishes

The historical extension now has complete independent printed-log concordance for all 100 h2 records and all 1,200 rg records, plus recovery-ledger current-hash agreement. The earlier missing-log evidence gap is closed at that level.

This is archival recovery and arithmetic/log validation. It is not a native replay, does not recover omitted unrounded estimator output or jackknife arrays, and does not reproduce raw-source filtering, genome-build mapping, allele harmonization or munging. Previously identified independence, heterogeneity-covariance, novelty and interpretation limitations remain in force. No MVP or other new validation outcome was examined by this review.
