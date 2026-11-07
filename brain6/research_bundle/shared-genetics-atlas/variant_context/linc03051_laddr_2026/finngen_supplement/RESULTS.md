# LINC03051 RNA features and FinnGen depression: published context

We inspected Table S3 from [Munro et al. 2026](https://doi.org/10.1016/j.ajhg.2026.06.005), which reports selected colocalizations between GTEx molecular features and FinnGen R12 traits. Twelve rows involve LINC03051, all for the paper's **depression/dysthymia** endpoint. Two rows match the tissue and residual RNA-feature identifiers of previously selected insomnia xTWAS results:

| Matched RNA feature | FinnGen depression–RNA CLPP | Prior insomnia xTWAS p | Prior insomnia–RNA coloc PP4 |
|---|---:|---:|---:|
| Hippocampus, residual PC2 | 0.0554 | 9.39×10⁻⁹ | 0.943 |
| Anterior cingulate, residual PC2 | 0.0157 | 3.13×10⁻⁸ | 0.761 |

These are **two different pairwise analyses**, not a three-trait colocalization or an ADHD result. The paper's FinnGen statistic is eCAVIAR-style CLPP, the sum of the products of the two credible sets' variant probabilities, with a reporting threshold of 0.01. It is not equivalent to the insomnia analysis's coloc PP4 and should not be compared numerically as though it were the same posterior. Matching a tissue and RNA feature identifier is useful context but does not independently verify all weights, processing or causal assumptions. The two analyses reuse GTEx molecular data, and overlap within each source and between tissues is not independently audited.

The FinnGen top shared variants for the two rows are GRCh38 chr3:117897197 C>T and chr3:117907302 C>T, respectively; neither is rs62264764. Their linkage to the insomnia or ADHD signals was not evaluated here. These findings suggest a testable RNA–depression hypothesis, but do not demonstrate that the same variant or LINC03051 mediates depression, insomnia or ADHD, and are not independent replication of both Brain6 traits. The prior nucleus-accumbens ADHD follow-up used residual PC14; the nucleus-accumbens depression row in this supplement is **PC2**, a different feature.

A [later reference-panel LD comparison](reference_ld/RESULTS.md) gives r²≈0.397 and 0.349 between these respective depression–RNA lead variants and rs62264764 in 503 European reference participants. This describes tagging in that panel; it is not a causal colocalization or study-matched LD estimate.

The published supplement contains 63,741 rows for 692 reported endpoint labels; none of those labels contains `ADHD` or `HYPERKIN`. The paper screened 1,302 FinnGen traits, and the supplement is restricted to results above its CLPP threshold. Therefore this search cannot establish whether ADHD was tested, whether another endpoint code represents it, or whether a result below threshold exists. No absence-of-effect claim follows.

The [source article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13361010/) describes the FinnGen R12 method. [All twelve gene rows](gene_rows.tsv), [the two identifier matches](matched_sleep_models.tsv), [extraction code](extract.py), [source receipt](receipt.json) and [verification](verification.json) record the evidence. The full Europe PMC supplemental archive and workbook are stored on the SSD; SHA-256 hashes document the exact downloaded bytes, but no independent publisher checksum was available. The first extraction attempt stopped because it used the wrong local model-set label; the corrected script produced the table. No AlphaGenome predictions were used to train a model.
