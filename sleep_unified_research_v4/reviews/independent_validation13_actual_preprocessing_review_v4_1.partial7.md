# Independent actual validation preprocessing review: partial 7/13

Review date: 2026-10-10. **QUALIFIED_PARTIAL_PASS for the seven completed sources below; the 13-source campaign is not terminal.** This is an actual-result review of the completed comparison, collector, adapter and QC metadata under the frozen physical-reference plan SHA-256 `83004c6817d4231b44d6d0d73725b4f5e8350de2d90231cf46579aec4a3117e9`. It reuses the independent original-pipeline/producer-profile reviews; it does not repeat source-body hashes, unchanged-library checks, controls, fits or preprocessing. Only this review and its companion JSON are written; no SSD or protected input is modified.

The consumed subset is plan members 1–7, with completed non-AppleDouble `*.comparison.json` receipts in `validation_pipeline_replay_v4/receipts_v4`. The companion JSON preserves each compact receipt/QC identity, actual metrics, independent arithmetic checks, reviewed source IDs and the six-member continuation list. These seven source result records are to be reused unchanged in a later append; only newly completed sources should be reviewed after a root follow-up.

## Actual content and QC findings

| Source | Retained SNP rows = finite Z/positive finite N rows | Assumed effective N, literal 12g string | Maximum absolute Z | Historical h² Z / intercept | Historical necessary QC |
|---|---:|---:|---:|---:|---|
| MVP GCST90479148, abdominal pain | 826,029 | 274973.347782 | 6.01833649535 | 15.603559 / 1.132467 | PASS retained |
| FinnGen R13 K11_REFLUX | 1,156,795 | 138755.739075 | 5.35188784625 | 15.457688 / 1.068438 | PASS retained |
| FinnGen R13 K11_DIAHER | 1,156,832 | 65341.5764931 | 12.2186023782 | 12.223836 / 1.077791 | PASS retained |
| FinnGen R13 K11_OTHENTERCOL | 1,156,827 | 46410.9882519 | 5.71063584075 | 5.184337 / 1.018352 | PASS retained |
| MVP GCST90479330, drug-allergy history | 826,031 | 2137.44946378 | 4.84255698668 | 0.706311 / 0.999432 | FAIL: h² Z<4 retained |
| FinnGen R13 J10_BRONCHNAS | 1,156,870 | 6041.44054262 | 4.89775474091 | 1.546938 / 0.993464 | FAIL: h² Z<4 retained |
| FinnGen R13 SMOKING_DEPEND | 1,156,810 | 13933.2207752 | 6.7302802044 | 5.486192 / 1.031445 | PASS retained |

The historical h²/intercept values above are taken from the already completed native validation h² table. They are not new fits or newly measured QC. Five sources retain the old necessary signal/intercept PASS, while two remain QC-ineligible. A large maximum SNP Z or successful preprocessing cannot rehabilitate a low-h² source or make it a negative control. These thresholds remain necessary historical boundaries, not a guarantee of population-scale h², transport power, joint contrast coverage or source admission.

The seven completed outputs contain **7,436,194 retained rows**, all reported as finite Z and positive finite N with **zero literal missing Z/N rows**. The reviewed comparator checks unique SNP identifiers, nonambiguous compatible A/C/G/T alleles, the exact `SNP/A1/A2/Z/N` header, full row-by-row ordered equality and full gzip CRC/EOF on each original and new derivative. Both ordered decompressed SHA-256 values agree for every source. Each observed retained count equals its frozen plan count and collector `output_rows`; no retained-source selection was changed.

I independently compared the old and new compact QC tables by metric membership/order, value and notes. Every scientific filter-count and effective-N row matches exactly. The `munged_output_sha256` row appropriately records each version's own differing compressed output identity; it is excluded from the claim that QC values are identical. The comparator's `original_QC_metric_rows` counts include that identity row: MVP has nine total rows/eight scientific rows, and FinnGen has eight total rows/seven scientific rows.

The aggregate input count is **146,044,064**, with **138,607,870 terminal rejection dispositions** plus **7,436,194 retained rows**, reconciling exactly. Every source reconciles individually. Rejections preserve the original non-HM3, allele/ambiguity, duplicate-rsID and MAF decisions. MVP additionally preserves its combined invalid/OR/P/CI/r²≤0.9 rejection category: 349,827 rows for abdominal pain and 349,824 for drug-allergy history. That combined category is not a count of INFO failures alone and cannot certify the clinical outcome's CI/test semantics.

The case/control effective size recomputed independently as `4/(1/cases+1/controls)` yields the exact frozen `.12g` N string in each QC table. This establishes reproduction of the assigned constant, **not observed per-variant enrollment**, verified per-variant N, population size or prevalence. The comparator confirms that the serialized assigned N is positive and finite in the retained stream; it does not change the source's sample-size interpretation.

## Byte identity, producer evidence and scientific meaning

**None of these seven compressed munged files is byte-identical to its archived predecessor. All seven have exactly equal ordered decompressed content.** The old/new compressed SHA-256 values differ, and the QC file's output-identity row changes accordingly. The observed evidence does not determine the cause of compressed differences; do not attribute them to timestamp, filename or gzip metadata without inspecting it. A claim of exact archival byte restoration would be false. The appropriate claim is literal ordered native-estimator input reproduction and unchanged scientific QC.

Completed collector/adapter metadata documents full source size/MD5/SHA-256 and gzip CRC/EOF verification for each source, and source identity agrees with the frozen member and original queue identity. This reviewer consumes that already executed producer evidence rather than freshly reading/hashing large bodies. The source-gate-only receipt by itself still says CRC/EOF is not established; the subsequent completed collector and comparison establish the stronger full replay evidence. All 21 completed source/collector/compare worker metadata records report return code 0, unchanged plan, no stop reason or metadata error, and verified owned cleanup. Observed peak owned RSS in this subset is 504,070,144 bytes, below the 2 GiB bound; this is subset evidence, not a full-family/terminal resource claim.

The new producer receipt retains historical literals such as `full_resolution_local_retention:false`; the additive adapter explicitly records that the exact source body is now locally retained. The producer's local route does not claim newly queried HTTP headers. This documented distinction is inherited from the reviewed producer profile, rather than treated as a new acquisition or remote identity test.

FinnGen's original rsID/reference-allele projection and constant effective N are reproduced. Its missing per-variant INFO and N remain missing; no INFO=1 or assumed INFO threshold is created. MVP's literal r² filter and CI-derived Z are reproduced, but an exact source-specific Wald/log-CI/effect/P contract is still unresolved. Ordered equality proves the original conventions were replayed; it does not independently validate those conventions, GRCh38-to-reference coordinate interpretation, Finnish LD adequacy, clinical phenotype equivalence, ascertainment, overlap or power.

This stage adds **zero association estimates, zero h²/r_g fits and zero independent replications**. It preserves the original 217-candidate membership, 41 historical estimates and their qualified classifications. Reusing the same sleep GWAS still means outcome-side external validation; neither a larger recovered file nor a successful pipeline supplies independent sleep enrollment, winner's-curse-free effect comparisons or a shared-trait sampling covariance. All prior phenotype/source/INFO/N/LD/clinical/rights qualifications remain.

## Remaining sources and stopping boundary

The following six frozen sources were not consumed by this partial review: `finngen_r13_J10_COPD`, `finngen_r13_M13_ARTHROSIS_KNEE`, `finngen_r13_M13_ARTHROSIS`, `finngen_r13_M13_INTERVERTEB`, `finngen_r13_K11_CHOLELITH`, and `finngen_r13_N14_CALCUKIDUR`. Their frozen retained-row expectation totals **6,940,986**, making the complete original family expectation 14,377,180. At the reviewed snapshot, PENDING exists and the terminal seal/master receipt do not. A source-only receipt or in-progress collector for the next member is not a completed comparison result.

No polling/waiting is performed after this snapshot. The next authorized update should append only newly completed source records and then, at terminal, assess the all-13 aggregate/terminal evidence once. Until that occurs this review must remain `QUALIFIED_PARTIAL_PASS`, not a terminal all-source certificate or new scientific admission.
