# Finnish-reference LD check at the TenK10K–frailty locus

The publicly documented [FinnGen LD service](https://docs.finngen.fi/working-outside-the-sandbox/ld-browser) can query SISu 4.2, a Finnish imputation-reference panel, for pairwise r². I queried the two SuSiE-coloc signal-index variants from the expanded TenK10K–FI comparison and the molecular peak's lead variant, using a 1 Mb window and an r² threshold of zero. These comparisons use GRCh38 variant allele strings, not rsID alone.

| Variant pair | SISu 4.2 r² | 503-person European reference r² |
|---|---:|---:|
| QTL fitted hit `5:136187682:T:G` vs FI fitted hit `5:136235970:T:C` | 0.3592 | 0.3170 |
| QTL fitted hit vs molecular peak lead `5:136131620:C:G` | 0.9942 | 0.9957 |
| FI fitted hit vs molecular peak lead | 0.3558 | 0.3152 |

The two fitted hits are moderately correlated in both references, whereas the molecular lead and fitted QTL hit are in very strong LD. These pairwise results are consistent with, but do not independently verify, the prior model's preference for separate signals. The FI study and TenK10K are not Finnish studies, so SISu is a complementary population check rather than demonstrably the better reference for that analysis. Neither r² nor a nearest gene assignment identifies a shared causal variant.

The service returned 1,991–2,058 of the 3,754 SNPs in the existing FI model per anchor; the requested window spans only part of that larger model region, and other availability or filtering differences may also contribute. Among returned exact matches, median absolute r² differences versus the EUR panel were 0.004–0.009; the 95th percentile was 0.102–0.106. These describe selected anchor pairs only and cannot establish agreement for the entire LD matrix. The service returns r² and D-prime, not the signed correlations needed for SuSiE-coloc; it does not fix the earlier model's elevated mismatch diagnostics, small external panel, rank deficiency, or uncertain sample-size assumptions. The source's sample count and relationship to study participants have not been independently verified.

The [machine-readable receipt](stage6_tenk10k_sisu_ld_comparison.json) records exact request URLs, response hashes, row counts and matched-set summaries. Individual responses and comparisons are `stage6_tenk10k_sisu_anchor_[1-3].json` and `stage6_tenk10k_sisu_anchor_[1-3]_comparison.csv.gz`; reproduction script: `work/tenk10k_sisu_ld_comparison.py`. The anchor itself is omitted by the service and has a null self-entry in the receipt, not a zero correlation. This is a public reference-data sensitivity check, not an independent participant-level validation of the score.
