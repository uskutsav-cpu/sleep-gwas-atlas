# Independent actual validation preprocessing review: terminal 13/13

Review date: 2026-10-10. **QUALIFIED_ALL13_TERMINAL_PASS for reproduction of the original filtering, serialization and QC.** All 13 frozen sources are complete under physical-reference plan SHA-256 `83004c6817d4231b44d6d0d73725b4f5e8350de2d90231cf46579aec4a3117e9`. This is a preprocessing result review with scientific qualifications; it establishes no new source admission, association or independent replication.

The prior seven-source Markdown and JSON are preserved verbatim as `independent_validation13_actual_preprocessing_review_v4_1.partial7.md` and `.partial7.json`. Their SHA-256 values are `34322ce1406452d44d526d8d5da7a439fef273a4ecdada851c3a57a1e62eee89` and `d9ad54e6d458e15aa88178f1918e90c738abe14429e7e6603507d0a6166c5fe3`. The first seven checked JSON result records remain unchanged. This update appends only plan members 8–13 and terminal evidence, using cached prior identities to close the final receipt without re-reading earlier source results. The pipeline and producer-profile reviews are reused. No GWAS/reference bodies, unchanged libraries, fits, controls or preprocessing are repeated; no SSD or protected BIO input is written.

## Retained content and source-level QC

| Source | Input rows | Rejected rows | Retained = finite Z / positive finite N rows | Assumed effective N, literal .12g | Maximum absolute Z | Historical h² Z / intercept | Historical necessary QC |
|---|---:|---:|---:|---:|---:|---:|---|
| MVP GCST90479148 abdominal pain | 19,703,221 | 18,877,192 | 826,029 | 274973.347782 | 6.01833649535 | 15.603559 / 1.132467 | PASS retained |
| FinnGen R13 K11_REFLUX | 21,326,506 | 20,169,711 | 1,156,795 | 138755.739075 | 5.35188784625 | 15.457688 / 1.068438 | PASS retained |
| FinnGen R13 K11_DIAHER | 21,326,252 | 20,169,420 | 1,156,832 | 65341.5764931 | 12.2186023782 | 12.223836 / 1.077791 | PASS retained |
| FinnGen R13 K11_OTHENTERCOL | 21,326,919 | 20,170,092 | 1,156,827 | 46410.9882519 | 5.71063584075 | 5.184337 / 1.018352 | PASS retained |
| MVP GCST90479330 drug-allergy history | 19,708,506 | 18,882,475 | 826,031 | 2137.44946378 | 4.84255698668 | 0.706311 / 0.999432 | FAIL: h² Z<4 retained |
| FinnGen R13 J10_BRONCHNAS | 21,325,654 | 20,168,784 | 1,156,870 | 6041.44054262 | 4.89775474091 | 1.546938 / 0.993464 | FAIL: h² Z<4 retained |
| FinnGen R13 SMOKING_DEPEND | 21,327,006 | 20,170,196 | 1,156,810 | 13933.2207752 | 6.7302802044 | 5.486192 / 1.031445 | PASS retained |
| FinnGen R13 J10_COPD | 21,326,060 | 20,169,181 | 1,156,879 | 95851.5846583 | 19.7119972182 | 14.543716 / 1.121515 | PASS retained |
| FinnGen R13 M13_ARTHROSIS_KNEE | 21,324,823 | 20,167,999 | 1,156,824 | 213523.320835 | 12.2344949368 | 22.640762 / 1.223616 | FAIL: intercept>1.2 retained |
| FinnGen R13 M13_ARTHROSIS | 21,325,719 | 20,168,896 | 1,156,823 | 317153.115034 | 11.4426725343 | 22.522232 / 1.219967 | FAIL: intercept>1.2 retained |
| FinnGen R13 M13_INTERVERTEB | 21,325,604 | 20,168,768 | 1,156,836 | 184776.518534 | 10.2385459997 | 18.823814 / 1.168145 | PASS retained |
| FinnGen R13 K11_CHOLELITH | 21,326,866 | 20,170,064 | 1,156,802 | 185876.661497 | 56.1942039552 | 5.741071 / 1.178645 | PASS retained |
| FinnGen R13 N14_CALCUKIDUR | 21,327,031 | 20,170,209 | 1,156,822 | 53057.1329907 | 11.2386142568 | 10.504049 / 1.061273 | PASS retained |

The full family reconciles exactly: **274,000,167 input rows = 259,622,987 terminal rejection dispositions + 14,377,180 retained rows**. All retained rows have reported finite Z and positive finite N, with zero literal missing Z/N. Each source reconciles individually and matches its frozen retained-row expectation. The newly consumed six contribute **127,956,103 input rows, 121,015,117 rejected rows and 6,940,986 retained rows**; the earlier seven totals are unchanged.

Rejection dispositions preserve the original non-HM3, allele/ambiguity, duplicate-rsID and MAF decisions. For the two MVP sources the combined invalid/OR/P/CI/r²≤0.9 category remains 349,827 rows for abdominal pain and 349,824 for drug-allergy history. It is not an isolated INFO-failure count. The companion JSON records every disposition count and the exact old/new QC metrics.

All 13 executed comparison receipts report full row-by-row ordered equality of the `SNP/A1/A2/Z/N` header and data, matching ordered decompressed SHA-256, unique SNP identifiers, compatible nonambiguous A/C/G/T alleles and full gzip CRC/EOF on each original and new derivative. These stream checks are consumed from the previously reviewed executed comparator; this reviewer independently checks the compact QC tables and arithmetic rather than opening or hashing large bodies again.

Every scientific QC row agrees by metric membership/order, value and notes. The compressed output SHA row correctly differs and is excluded from scientific-QC equality: MVP has nine total/eight scientific rows and FinnGen has eight total/seven scientific rows. Independently recomputing `4/(1/cases+1/controls)` for each new source reproduces the exact frozen `.12g` N string. This validates reproduction of the assigned constant, not observed per-variant N or enrollment.

The h²/intercept values in the table come from the already completed native validation table; none was fitted in this stage. **Nine sources retain historical necessary QC PASS; four remain QC-ineligible:** drug-allergy history and J10_BRONCHNAS have h² Z<4, while M13_ARTHROSIS_KNEE and M13_ARTHROSIS have intercept>1.2. Strong individual SNP Z, complete acquisition or successful preprocessing cannot reverse these exclusions. The PASS threshold is necessary historical QC, not evidence of variance-fraction interpretation, transport power or new source admission.

## Identity and terminal evidence

**All 13 have identical ordered decompressed content; none is compressed-byte-identical to its archived predecessor.** Old/new compressed SHA-256 values differ. The comparison establishes literal ordered native-estimator input reproduction and unchanged scientific QC, not exact archival byte restoration. The cause of the compressed difference is not inferred.

Completed producer/adapter evidence records full source size/MD5/SHA-256 and source gzip CRC/EOF verification, consistent with frozen member/queue identities. A source-gate-only receipt still truthfully marks CRC/EOF and raw replay incomplete; the later completed collector/adapter/comparison supplies the stronger evidence. Historical `full_resolution_local_retention:false` literals are preserved while adapters explicitly document exact local retention and a local route that does not claim newly queried HTTP headers.

The final execution receipt SHA-256 is `6bc9dc1700db14ef2a6a7f4bd1cb4454eaf41a9ab5aa5cd2726c54c8877e3c7f`; terminal seal SHA-256 is `a5a14dee54f4be03d2bdc249ccc3756834f08bea3681dd856dd61eebde961b89`. Both are recomputed from compact metadata and match the reported identities. The seal binds that exact receipt, the frozen plan and admission SHA-256 `fe413af06e6b9d7569554326944935886166aeff503372b17724d59861515bcc`. Frozen source membership and order match all 13 completed comparisons, and compact comparison/producer/QC/worker identities close against the master receipt. The earlier seven are checked through their preserved review identities; their result records are unchanged.

Both the exact pending path `validation_pipeline_pending_v4.json` and literal `PENDING` are absent, and no non-AppleDouble failure/abort addenda are present. Worker journal filenames containing `failure_journal` are routine journals, not evidence of a failed run. All 39 completed source/collector/compare worker records report return code 0, unchanged plan and verified owned cleanup; new workers have no stop reason or metadata errors. The master records owned cleanup and no termination requests. The parent reports session21520 exited0; this review does not query the completed process.

Observed peak owned RSS for this family is **640,827,392 bytes**, below the 2 GiB bound. This is worker-metadata evidence; no new resource census or controls were run. The inherited 190-fit gate receipts in the execution master were not re-audited.

## Scientific interpretation and stopping boundary

FinnGen’s rsID/reference-allele projection and assigned constant effective N are literally reproduced; absent per-variant INFO and N remain absent. No INFO=1 or INFO threshold is invented. The frozen queue’s EUR shorthand is preserved rather than treated as a newly certified population description: these are Finnish-ancestry sources, and adequacy of the external EUR LD reference remains unresolved. MVP’s original r² filter and CI-derived Z are reproduced, but its accession-specific CI/effect/P/test contract remains unresolved. Ordered equality validates replay of the conventions, not their phenotype, statistical or coordinate interpretation.

Clinical definition/ascertainment equivalence, build/effect semantics, INFO/N, Finnish LD adequacy, rights, cohort membership and power are not cleared by recovery or preprocessing. All four discovery-sleep/discovery-outcome/validation-sleep/validation-outcome participant intersections still need the prior human cohort review. Cohort names do not demonstrate individual or relative independence. Reusing discovery sleep leaves outcome-side external validation; the shared-sleep sampling covariance and winner’s-curse qualifications persist.

This stage produces **zero new associations, zero new h²/r_g fits and zero independent both-trait replications**. Original 217-candidate membership, 41 historical estimates and their qualifications remain unchanged. No sources remain unconsumed in this original-13 family. This review is complete for the current finite assignment; the newly launched core27 continuation and queued extension100 family are outside its certificate.
