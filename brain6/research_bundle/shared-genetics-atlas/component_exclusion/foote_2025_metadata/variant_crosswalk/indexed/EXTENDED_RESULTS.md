# Published frailty-factor lookup: 87 selected variants

We queried the seven published Foote factor files at 87 post-selected project variants after uniquely lifting their GRCh37 coordinates to GRCh38. These include the FI/FFS diagnostic crossings and five previously selected anchors. The table reports exact rsID, position and allele matches in the harmonized release; all remaining rows were absent, not invalid or negative associations. This is a post-selected descriptive crosswalk, not an independent replication cohort or a newly corrected association scan.

| Foote factor | Exact matches / 87 | Missing / 87 | Matched markers below the paper's 7.14×10⁻⁹ seven-factor threshold |
|---|---:|---:|---:|
| General frailty | 32 | 55 | 25 |
| Limited social support | 76 | 11 | 0 |
| Unhealthy lifestyle | 76 | 11 | 0 |
| Multimorbidity | 69 | 18 | 0 |
| Metabolic problems | 39 | 48 | 0 |
| Poorer cognition | 76 | 11 | 0 |
| Disability | 64 | 23 | 0 |

All 25 general-factor markers lie on chromosome 5: 23 are FFS diagnostic markers and two are the prior anchors rs2431108 and rs7703746. They are linked regional markers, not 25 independent loci. The exact-marker general-factor association at rs1363103 is p=9.16×10⁻²¹; at rs13164188 it is p=1.07×10⁻²⁰. Both also appear in the post-selected FI–insomnia and FFS follow-up. Their appearance in several analyses does not supply independent confirmation: the factor includes insomnia, and frailty and sleep input data overlap. The factor model is not a matched FI-minus-insomnia or FFS-minus-exhaustion score.

The factor files were filtered to remove variants with Q heterogeneity p<7.14×10⁻⁹, among other possible input and processing exclusions. The 55 general-factor absences cannot be assigned to a particular cause without unfiltered files; no association strength can be inferred for them. The apparent zero marker counts for other factors are limited to this selected set and threshold, not evidence of no regional or genome-wide relationships. The significance cutoff is the paper's stated factor-GWAS criterion and is not a correction for our post-hoc selection, additional sleep phenotypes or claims of shared causality.

All source fields and statuses for 609 marker–factor combinations are in `all_extended.tsv`; `extended_summary.json` records counts. The five-anchor results from a separate indexed query were checked to have identical statuses and source values. The lookups rely on remote indexed GRCh38 data and **do not validate whole-file checksums**; the separately running complete-stream checks concern different, original GRCh37 source files. The [protocol](EXTENDED_PROTOCOL.md), factor receipts and saved metadata document the liftover and queries. The [earlier factor review](../../REVIEW.md) explains how these latent constructs differ from component-excluded analyses.

A [neighborhood inspection](ABSENCE_RESULTS.md) finds that 50 of 55 absent general-factor markers have no released records within ±250 kb. The paper describes pruning around heterogeneous signals, but the cause of each gap remains unknown.

**General-factor verification update:** a later [complete original-file read](../VERIFIED_GENERAL_FACTOR.md) passed the publisher MD5 and reproduced all 32 matches and 55 absences for this factor; all matched effect, frequency, SE, p and Q values agree after allele orientation. The harmonized whole file and the six other original factor files remain unverified. An absence from both filtered files is still not proof of absence from the unfiltered input or of no association.
