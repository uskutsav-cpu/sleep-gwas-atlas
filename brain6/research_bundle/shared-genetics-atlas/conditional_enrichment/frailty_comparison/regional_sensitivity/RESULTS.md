# Regional sensitivity of the three-frailty conditional comparison

This post-hoc analysis removes whole chromosomes or a broad APOE-containing window from the previously matched marker set and recalculates the descriptive insomnia p≤0.01 fractions among frailty p≤0.001 markers. It reuses physical-block counts; the [protocol](PROTOCOL.md) defines the exact windows. These are SNP counts, not independent loci or tests of a causal gene.

| Excluded from matched support | FI-conditioned fraction | FFS-conditioned fraction | HFRS-conditioned fraction |
|---|---:|---:|---:|
| Nothing | 6.13% | 5.64% | 8.60% |
| chr19:45,000,001–50,000,000 GRCh37 | 6.13% | 5.65% | 8.62% |
| All chromosome 19 | 6.17% | 5.67% | 8.65% |
| All chromosome 5 | 5.11% | 4.90% | 8.35% |
| All chromosome 11 | 6.02% | 5.62% | 7.93% |

The broad chr19 window contains only 2 of the 1,459 HFRS-conditioned markers also crossing the insomnia threshold; removing all chr19 removes 12. Thus the descriptive HFRS fraction in this analysis is not primarily driven by that one window. Chromosome 5 contributes 129 of 553 FI joint markers and 134 of 684 FFS joint markers; chromosome 11 contributes 217 of 1,459 HFRS joint markers. Linked SNPs may count the same signal many times, and this scan does not establish that any chromosome or gene causes the observed differences. The exclusion choices were inspected after the original analysis and are not independent replication or a multiplicity-corrected test.

`exclusion_results.tsv` contains all leave-one-chromosome results and the broad chr19-window check, including denominators and removed counts. `hfrs_block_counts_ranked.tsv` retains the full 562-block table sorted by observed HFRS joint counts, while `receipt.json` and `run.py` record the source hash and computation. These checks do not remove FI sleep items or FFS exhaustion, address FinnGen sample overlap, perform conjFDR, or validate a shared causal locus. See the [common-support comparison](../RESULTS.md) for its broader limitations.
