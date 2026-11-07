# MVP insomnia and FinnGen clinical ADHD: six-lead comparison

Among the six original Brain6 regional leads, five have exact matching records in both the previously checked MVP insomnia source and the newly queried FinnGen R13 clinical ADHD source. **None crosses the descriptive six-lead 0.05 benchmark** for the two associations jointly. The sixth lacks an exact indexed FinnGen row and remains in the denominator.

| Brain6 lead | MVP insomnia p | FinnGen ADHD p | max-p × 6, capped at 1 |
|---|---:|---:|---:|
| rs2582894 | 0.000933 | 0.0653 | 0.392 |
| rs7160630 | 0.000791 | 0.0133 | 0.0798 |
| rs62264767 | 0.0000499 | 0.0261 | 0.157 |
| rs2431108 | 0.0895 | 0.000108 | 0.537 |
| rs6452785 | 0.0000310 | Missing | Not estimable |
| rs1476535 | 0.000000120 | 0.367 | 1.000 |

The chromosome 3 lead rs62264767 has strong MVP insomnia evidence and nominal FinnGen ADHD evidence, but the six-lead numerical benchmark does not pass. Its closely tagged, separately selected neighbor **rs62264764** has MVP insomnia p = 0.00004918 and FinnGen ADHD p = 0.0279; this is an exploratory pair of marginal observations, not an additional independent locus or a prespecified seventh test. The chromosome 5 lead rs2431108 has stronger ADHD-side evidence but weaker MVP insomnia evidence. A non-passing p-value does not demonstrate absence of association.

These calculations were made after the results and markers had already been examined. Max-p is an intersection–union test of both marginal associations if the marginal p-values are valid; it does not require the two statistics to be independent. Multiplication by six is **not** a confirmed familywise correction for the broader discovery and subsequent choices. The source comparisons do not show that both traits share a causal variant. The sleep and ADHD effects have different units and case definitions; variant imputation quality and ancestry-specific LD need attention.

MVP and FinnGen are separately named sources, but individual and family overlap with the discovery studies and between all sources has not been independently audited. The FinnGen index was fully hashed, and a [separate whole-source byte stream](../finngen_adhd_raw_lookup/whole_source_hash_receipt.json) matched the historical source length and MD5; the compressed source was not retained or fully decompressed. The variant query itself used indexed ranges. See the [MVP source analysis](../mvp_six_leads/RESULTS.md), [FinnGen clinical ADHD lookup](../finngen_adhd_raw_lookup/RESULTS.md), and [rs62264764 MVP analysis](../rs62264764_mvp/RESULTS.md) for phenotype, matching and source limits. These results do not validate AlphaGenome, establish a gene mechanism or repair failed LAVA quality gates.

[Protocol](PROTOCOL.md), [code](analyze.py), [joined results](results.tsv), and [input-hash receipt](receipt.json) record this post-hoc comparison. Both constituent tables trace to their own receipts; matching a file hash does not verify cohort independence or fix selection bias.
