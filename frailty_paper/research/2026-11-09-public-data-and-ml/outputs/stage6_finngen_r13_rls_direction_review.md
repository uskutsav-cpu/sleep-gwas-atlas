# FinnGen R13 RLS: two-variant follow-up

We queried the public [FinnGen R13 RLS summary file](https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_G6_RLS.gz) using its tabix index for two variants selected after a cross-study direction discrepancy. FinnGen's [data documentation](https://finngen.gitbook.io/documentation/data-description) identifies `alt` as the effect allele and `beta` as its log-odds effect. The [RLS endpoint description](https://risteys.finngen.fi/endpoints/G6_RLS) uses register-based codes; ascertainment need not match the other RLS studies.

| Variant and effect allele | FinnGen R13 beta (SE) | p | Effect allele frequency |
|---|---:|---:|---:|
| rs113851554 T | +0.691794 (0.0304694) | 4.04 × 10⁻¹¹⁴ | 0.0698 |
| rs78206187 G | +0.456858 (0.0503174) | 1.09 × 10⁻¹⁹ | 0.0277 |

Both matched the stated GRCh38 coordinate, ref/alt pair and rsID. The positive rs78206187 G direction agrees with the EU-RLS-GENE and INTERVAL summary results and the conditional current All of Us browser result, but differs from the negative Akçimen meta-analysis. This strengthens the observation of a **source discrepancy**, not its explanation. The Finnish estimate and other estimates have different phenotype definitions, selection, allele frequencies and possible participant/relative overlap; EU-RLS-GENE includes Finnish cohorts. Do not pool them or count these findings as independent replications without an overlap audit. The candidate set and endpoints were postselected. Neither RLS summary results nor these variants validate the project scores or establish a shared frailty mechanism.

[Lookup](stage6_finngen_r13_rls_two.csv), [source and hash receipt](stage6_finngen_r13_rls_two_receipt.json); reproduce with `python3 work/rls2024/finngen_r13_two.py`. The query retrieves indexed ranges only, so the complete 803 MB compressed file's published checksum was not verified. No individual-level data were accessed. Finnish reference allele frequencies differ from the other public studies, so transportability and LD should not be assumed.
