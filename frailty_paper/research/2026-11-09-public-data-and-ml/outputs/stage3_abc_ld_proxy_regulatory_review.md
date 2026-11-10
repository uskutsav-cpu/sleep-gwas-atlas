# LD proxy overlap with ABC regulatory predictions

The [complete overlap table](stage3_abc_ld_proxy_regulatory_links.csv) contains 1,531 source prediction records for proxies of 13 candidate markers. They represent 147 candidate–target-label pairs across 128 source cell/tissue labels; records may be repeated for correlated proxies, intervals and tissues. The input was 726 nonself candidate–proxy records with European reference r² ≥ 0.8 within 250 kb. A source record is included when a proxy position overlaps an ABC interval; see [the receipt](stage3_abc_ld_proxy_overlap_receipt.json) and [LD methods](stage3_candidate_reference_ld_proxies.md).

| Candidate | ABC records | Distinct source target labels | Distinct reference proxies | Source cell/tissue labels |
| --- | ---: | ---: | ---: | ---: |
| rs12524473 | 107 | 6 | 10 | 41 |
| rs12536335 | 7 | 1 | 1 | 7 |
| rs12671330 | 6 | 3 | 2 | 1 |
| rs1909650 | 1 | 1 | 1 | 1 |
| rs2155292 | 83 | 28 | 3 | 19 |
| rs2604551 | 535 | 22 | 3 | 104 |
| rs2815764 | 171 | 12 | 4 | 64 |
| rs2838788 | 25 | 6 | 3 | 7 |
| rs28528997 | 18 | 10 | 2 | 4 |
| rs34290943 | 168 | 19 | 2 | 56 |
| rs6058750 | 141 | 31 | 10 | 28 |
| rs7734045 | 268 | 7 | 28 | 93 |
| rs817509 | 1 | 1 | 1 | 1 |

For example, a proxy for rs2815764 overlaps NEGR1-labeled promoter predictions in multiple source contexts, and rs7734045 proxies overlap predictions naming several target labels. These are starting hypotheses, not proof that the index marker or proxy regulates the target in sleep or frailty. ABC scores are model scores, not probabilities of disease causality. The source includes promoter, genic and intergenic classes, and some source cell labels are fetal, cultured, non-neural or otherwise not established as disease-relevant. For 110 of 147 candidate–target pairs, the source symbol, chromosome and TSS match one gene among the tested Ensembl v92 protein-coding genes; 37 pairs remain unresolved in this limited mapping. See [the identifier mapping](stage3_abc_ld_proxy_gene_id_mapping.csv). This is not an exhaustive annotation or functional validation. The thresholded source does not cover all tissues, and absence of a record is not negative evidence.

The reference is 503 European individuals; LD may differ across ancestry or the GWAS samples. Multiallelic representation, the GRCh37 coordinate assumptions and source 0-based interval convention remain important limitations. This procedure did not intersect GWAS credible sets with predictions or show allele-specific effects. The underlying ABC source and file hash are documented in the [exact-marker ABC review](stage3_abc_regulatory_links_review.md). Each new candidate–target pair in the variant–gene table is explicitly labeled as a model prediction, and the network keeps LD and regulatory links as distinct edges.
