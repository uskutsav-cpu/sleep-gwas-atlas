# Independent statistical review of V7 benchmark results

**Scope:** read-only review of the serialized comparator results and the frozen V4 records; it did not rerun GWAS/LDSC or change any results. This is an automated independent-agent review, not human peer review or empirical replication.

## Data review and correction

The first V7 result table had 6,548 rows (1,637 records × four comparators) and 361 qualified rows, but 30 Human GWAS ATLAS rows were class A using publication PMID, approximate sample size, ancestry label, and phenotype label. Review of the actual HGA pair table found opposite rg signs for 16 of those 30 candidates, including sleep-duration/Crohn's disease and insomnia/IBD examples. An aligned insomnia/MDD example showed that the comparison was discriminating. PMID/N/population matching alone did not establish source release, phenotype coding or ascertainment, effect direction, or estimator compatibility. The final builder therefore classifies all 30 as E, insufficient evidence.

## Final count and integrity checks

| Comparator | A | B | C | D | E | Candidates | B/C/E candidates |
|---|---:|---:|---:|---:|---:|---:|---:|
| Human GWAS ATLAS | 0 | 210 | 0 | 1,344 | 83 | 240 | 240 |
| Morrison | 0 | 888 | 0 | 749 | 0 | 888 | 888 |
| Goodman | 0 | 492 | 0 | 815 | 330 | 822 | 822 |
| Fan | 0 | 40 | 0 | 1,549 | 48 | 88 | 88 |
| Union | 0 | 1,183 | 0 | 366 | 88 | 1,247 | 1,247 |

Every candidate-failure proportion is 100%. The 30 reclassified HGA rows are Morrison B, so the union precedence retains them as B. The 361 qualified match-absent IDs (85 core, 276 extension) exactly match an independent recomputation of the source/QC predicate; each is D in all four comparator columns. Candidate reference counts equal `candidate_result_count` for all comparators, with no duplicate references. The Goodman SD24 parser covers six scores (SHS-ADD and PC1–PC5). The four manifest-listed generated TSV hashes match the serialized tables.

The 100% candidate-failure proportion measures source/estimand matching difficulty, not a statistical failure rate. Class D remains specific to a named snapshot. No new biological result follows from this comparison.

## Material limits

The cross-resource result does not verify a global search, participant independence, statistical effect-size concordance, or public redistribution rights. Fan row-reference derivatives require human rights review. SleepChart is not in the four-snapshot union; CTG-VL has no direct export here. The original V4 validation set shares discovery sleep sources and is not fully independent two-trait replication.
