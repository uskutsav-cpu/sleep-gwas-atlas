# Two selected variants: Estonian insomnia and Finnish clinical ADHD

The two Brain6 candidates appearing in the selected Estonian insomnia table both have exact coordinate, rsID and allele-pair matches in the public FinnGen R13 clinical ADHD release. Neither supports a nominal association with **both** phenotypes in this post-selected comparison:

| Variant | Estonian insomnia p | FinnGen ADHD p | Larger p |
|---|---:|---:|---:|
| rs1093463, chromosome 3 | 0.539 | 0.0384 | 0.539 |
| rs1989903, chromosome 7 | 0.465 | 0.156 | 0.465 |

The rs1093463 ADHD result is nominal only and was selected after prior Brain6 and sleep-table examination; it is not a corrected or independently verified replication claim. Its T-allele Finnish ADHD beta is +0.0504 (SE 0.0244). The rs1989903 G-allele beta is +0.0281 (SE 0.0198). Both Estonian insomnia estimates have wide confidence intervals that include zero. The larger p-value is an intersection–union statistic only if both marginal p-values are valid; it is not a colocalization test and the table is not a prespecified family.

The Estonian sleep result comes from an unreviewed preprint's **selected-lead table**, not a complete regional GWAS. Its registry phenotype differs from Brain6's original questionnaire measure. FinnGen's endpoint is clinical ADHD and the source descriptions report 5,559 cases and 489,493 controls. Participant and family overlap across all discovery and follow-up cohorts has not been independently audited. No two-trait replication, shared causal variant or gene mechanism follows, and a weak result does not prove absence of an association. The earlier six-lead assessment and failed LAVA gates remain unchanged.

The [protocol](PROTOCOL.md), [results](results.tsv), [source receipt](receipt.json), [lookup](lookup.py) and [independent indexed-row check](verification.json) document this comparison. The lookup uses a version-pinned public source and hash-matched tabix index; a previous [whole-source byte stream](../../finngen_adhd_raw_lookup/whole_source_hash_receipt.json) matched the recorded MD5 but did not fully decompress the gzip. The independent verification refetched both source rows; it does not validate the source model or ascertainment.
