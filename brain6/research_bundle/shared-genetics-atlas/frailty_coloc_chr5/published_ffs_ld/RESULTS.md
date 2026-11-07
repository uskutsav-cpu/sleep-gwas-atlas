# LD of the published chromosome 5 FFS lead

The previously published FFS lead rs7703746 and the chr5 candidates were queried in the phased 1000 Genomes Phase 3 GRCh37 European reference, with 503 EUR-labeled individuals. All four rsIDs, positions and allele pairs matched unique reference records. The analysis was post-hoc; [PROTOCOL.md](PROTOCOL.md) specifies the comparison.

| Reference pair | Counted alleles | r | r² | Individual-bootstrap 95% interval for r² |
|---|---|---:|---:|---:|
| rs7703746–rs1363103 | G–T | +0.727 | 0.529 | 0.477–0.580 |
| rs7703746–rs13164188 | G–T | +0.733 | 0.537 | 0.486–0.590 |
| rs7703746–rs2431108 | G–T | −0.734 | 0.538 | 0.489–0.590 |

These are reference-panel LD estimates, not GWAS effect correlations. The bootstrap resamples the 503 reference individuals; it does not account for differences from the UKB, TwinGene or Finnish study samples. The FIN subgroup is not separately assessed here. Moderate LD can explain why multiple SNP associations may coexist in the region, but cannot establish that insomnia, FI, FFS and ADHD share one causal variant. The Ye discovery table identifies rs7703746 as an FFS lead; a later [adapter audit](../../replication/ye_2023/allele_adapter/RESULTS.md) supports the source effect-allele convention, while precise interpretation of its phenotype scale remains important. Nor does this calculation establish independent two-trait replication or a target gene.

`pairwise_ld.tsv` contains frequencies, sample counts and estimates; `reference_receipt.json` records the reference records, remote URLs and panel hash, and `calculate.py` reproduces the calculations with the documented external panel. The [source and replication review](../../replication/ye_2023/REVIEW.md) and [FI–insomnia coloc](../RESULTS.md) provide the distinct association evidence.
