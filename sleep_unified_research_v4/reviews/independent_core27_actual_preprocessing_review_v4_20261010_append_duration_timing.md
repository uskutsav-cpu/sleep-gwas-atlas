# Independent core actual-output append: duration and timing

2026-10-10. **QUALIFIED_TWO_SOURCE_APPEND_PASS** for only `accel_sleep_duration` and `sleep_timing`, under already admitted plan SHA-256 `8a839122ef7a9a0faefb6c3971f04abd8fa755ed630a6db1efe8f725ca54c6e3` and admission `11f00f5e7951da0a3d44e4587fff384d60e66b8f3122af0f0993bd92ab8cd136`. This new append leaves the previous partial review unchanged. It is not a family terminal certificate.

| Trait | Input rows | Harmonization rejects | Harmonized rows | HM3 template rows | Finite N/Z | Missing/nonfinite N/Z | Assigned constant N |
|---|---:|---:|---:|---:|---:|---:|---:|
| `accel_sleep_duration` | 11,997,351 | 5,514,837 | 6,482,514 | 1,217,311 | 1,176,749 | 40,562 | 85,449 |
| `sleep_timing` | 11,997,351 | 5,514,811 | 6,482,540 | 1,217,311 | 1,176,759 | 40,552 | 84,810 |

**No scientific content or QC discrepancy was found.** The executed full-stream comparator reports zero unequal rows in both harmonized and complete HM3 outputs, matching ordered decompressed hashes, exact SNP/allele/N/Z literals and missingness, and complete gzip CRC/EOF verification. No numerical tolerance was introduced. I independently compared the actual and archived compact QC metadata and all 15 ordered drop/remaining steps for each source; counts reconcile individually. Only execution paths and five additive `not supplied` absence fields differ. The 1,217,311-row HM3 left-join skeleton is not a finite-statistic count.

Stock-munge scientific counter lines also match the archived logs after excluding output paths. `accel_sleep_duration`: 5,305,528 outside the merge list and 237 allele mismatches; `sleep_timing`: 5,305,546 outside the merge list and 235 allele mismatches. Each reconciles to its harmonized rows and finite HM3 count. Harmonized and munged compressed SHA-256 values differ in both traits; exact ordered content does not imply archival byte identity, and the compression difference cause is not inferred.

Both traits preserve the Jones GCST007803 European UK Biobank source, hg19 coordinates with no liftover, ALLELE1 effect orientation, and their distinct rank-normalized `ACC_SLEEP_DUR_RAW_SIN` and `ACC_SLEEP_MIDP_RAW_SIN` column families. Timing retains the SPT-window midpoint phenotype definition. N remains the frozen assigned constant, not verified per-variant enrollment. The original INFO (0.6,1] gate rejects 122,113 duration rows and 122,142 timing rows. INFO is omitted from the final harmonized schema; a stock log of zero INFO≤0.9 removals cannot establish INFO>0.9. No revised source, N, INFO, phenotype or threshold contract is admitted.

All ten completed source/harmonize/munge/compare/cleanup worker records report return code 0, unchanged plan, no stop reason or metadata error, and owned cleanup. Subset peak observed RSS is 1,246,396,416 bytes. Source-gate and comparison-worker compact identities agree. This review consumes executed full-stream evidence and reads only compact metadata/QC/logs; it performs no large output/raw/reference scans, unchanged-library audit, fits, controller changes, admissions, SSD/BIO/code writes or active-worker inspection. Declared workflow agreement does not attest identical historical per-trait binaries.

This is preprocessing reproducibility only: no new association, discovery or independent replication. Selected UK Biobank device sampling, construct/clinical equivalence and shared-cohort/shared-trait covariance limitations persist. Family27 remains active under the parent’s status; later unconsumed outputs and terminal evidence are outside this append.

Compact evidence identities (SHA-256):

| Evidence | SHA-256 |
|---|---|
| `accel_sleep_duration` comparison | `be9814f5d5cd31fd19fb350ad296a4e25028ca5f55d90842fe8be3bf00e4b2ce` |
| `accel_sleep_duration` new QC | `45440eb058ace6df53a05c9ad1f2b332c74bf4e1cd30e47bfaff9da2f1ce40e4` |
| `accel_sleep_duration` archived QC | `0be686f4263155ad507d29f869e6d8c4d48210552412f2f9a5e7c7043426d928` |
| `sleep_timing` comparison | `bf5a6372a8c7fcb66238217233118de71c27b3b99e0ee1a7b0c67210b8777c63` |
| `sleep_timing` new QC | `4809b694950106653959680e03f3f13a4dbfb50c4a369b4134d2d2c4f28ccf15` |
| `sleep_timing` archived QC | `52a8eb9f36fb4738c49655726a614b324e4d0085e2e42c913e4392926fd58d1a` |
