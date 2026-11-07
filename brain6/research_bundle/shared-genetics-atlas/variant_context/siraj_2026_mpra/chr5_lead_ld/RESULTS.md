# Chromosome 5 reporter candidate and shared-trait leads: reference LD

The post-selected reporter-assay candidate **rs12658032** is linked to the previously studied chromosome 5 insomnia–ADHD leads in the 1000 Genomes Phase 3 European-labeled reference sample:

| Pair | r² | 95% bootstrap interval for r² | Counted alternate alleles |
|---|---:|---:|---|
| rs12658032 – rs2431108 | 0.776 | 0.725–0.829 | G and C; r = −0.881 |
| rs12658032 – rs77960 | 0.773 | 0.718–0.826 | G and A; r = −0.879 |
| rs2431108 – rs77960 | 0.995 | 0.986–1.000 | C and A; r = +0.998 |

All three pairwise estimates use 503 reference individuals. Phase-aware haplotype correlations use the listed alternate alleles; 2,000 bootstrap resamples are over individuals. The resulting LD supports investigating rs12658032 in this locus rather than treating it as an unrelated regional marker. It does **not** prove it is causal, establish perfect tagging in discovery or Finnish cohorts, identify a gene, or validate a shared disease mechanism. Bootstrap intervals reflect this reference panel only, and the remote whole VCF was not independently hashed.

The association inputs show the G allele of rs12658032 has negative insomnia and ADHD beta estimates, whereas the C allele of rs2431108 has positive estimates; the negative signed reference LD is consistent with these directions being on opposite haplotypes. This is descriptive consistency, not independent replication. In the selected single-signal Brain6 approximations, rs12658032 has insomnia PIP≈0.086 and ADHD PIP≈0.0046; the numbers are model-conditional, not evidence that a variant causes either trait. Its adult-frontal 501-bp AlphaGenome DNase ALT/REF prediction is −0.0201 log2, a prediction that is not a measured reporter effect or a risk estimate.

The [assay follow-up](../RESULTS.md) reports a skew FDR≈0.017 for the OL33 SK-N-SH library but ≈0.756 for OL32; the OL33 row did not meet the study's specified per-row activity threshold. We do not classify it as the authors' emVar or claim cross-library replication. Differences in cell context, sample and sequence library, prior selection, multiple testing and source activity rules remain material. The rs62264764 finding is on chromosome 3 and is distinct from this analysis; the original LAVA gate remains failed.

[Pairwise table](pairwise_ld.tsv), [source receipt](receipt.json) and [code](analyze.py) preserve reference allele mappings, counts, bootstrap seed and panel hash. Source-association and prior inputs are recorded in the [rank comparison](../../../prediction_rank_comparison/variant_ranks.tsv) and [dense chromosome 5 QTL match](../../chr5_dense_qtl/Testis_ENSG00000251574_insomnia_matched.tsv).

An additional [developmental caQTL and expression lookup](../chr5_developmental_lookup/RESULTS.md) found no nominal association among 14 selected developmental variant–peak rows, while the available expression association is in testis. Neither finding identifies a brain effector gene.
