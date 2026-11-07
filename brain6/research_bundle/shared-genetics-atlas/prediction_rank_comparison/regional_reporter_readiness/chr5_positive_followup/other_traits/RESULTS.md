# rs12659965 in additional sleep and hospital-frailty sources

The exact rs12659965 G/A contrast appears in both additional sources. The reported **MVP insomnia p=0.01019** and **FinnGen hospital-frailty p=0.02391** are nominal, post-selected associations; they are not corrected for the wider search.

| Source | Source effect and p | A-allele direction |
| --- | --- | --- |
| MVP clinical-record insomnia | G-allele OR 1.014; p=0.01019 | Approximate log-odds −0.0139, from the rounded OR |
| FinnGen R12 HFRS | A-allele beta −0.02840, SE 0.01257; p=0.02391 | Negative |
| FinnGen R13 registry insomnia, previously checked | A-allele beta −0.02803; p=0.00002635 | Negative |
| Childhood ADHD symptoms, previously checked | A-allele beta −0.00657; p=0.01575 | Negative |

The signs agree in these four displayed estimates, but the phenotypes, populations and effect scales differ. MVP reports no SE in this row, and its OR is rounded; its p-value is taken directly from the source. FinnGen R12 and R13 have overlapping participants, and exact cohort and relative overlap across all sources has not been audited. HFRS is hospital-ascertained and is not the UKB Fried score. The separately observed UKB Fried score association and uncertainty about its effect-allele convention are documented in the [parent review](../RESULTS.md).

The maximum of the MVP, HFRS and symptom marginal p-values is **0.02391**. That is only a descriptive intersection–union statistic if all three input p-values are valid. It is not adjusted for the many candidates, outcomes and prior searches that led to this row; it does not establish a shared causal variant, independent three-trait replication, a gene, or a disease mechanism. Weak or nominal results do not establish either presence or absence of a true effect. The selected SK-N-SH reporter observation and moderate reference LD provide experimental and regional context, not evidence that the same variant mediates all traits. The failed local-QC gates remain unchanged.

The [pre-query protocol](PROTOCOL.md), [full-file scan](query.py), [all matching source rows](receipt.json), [verification](verification.json), and [run log](run.log) record the lookup. Both source files were read completely and matched their previously saved hashes (19,703,815 MVP rows and 21,294,560 HFRS rows); exact GRCh38 G/A matches were unique. These checks do not validate the source GWAS models, cohort independence or all possible multiallelic representations. No individual data or AlphaGenome predictions were used.
