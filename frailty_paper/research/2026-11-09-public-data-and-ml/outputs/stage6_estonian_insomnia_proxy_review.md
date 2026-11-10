# Estonian insomnia: exploratory nearby proxies

The [Estonian follow-up table](https://pmc-oa-opendata.s3.amazonaws.com/PMC12668095.1/media-2.xlsx), sheet ST5, reports an insomnia association for rs1093463, 6,099 bp from 22-marker candidate rs817509 on GRCh37. In 503 1000 Genomes phase 3 European reference individuals, their **alternate-allele dosages are identical in this sample** (r² = 1.00). The reference VCF has candidate A/G and proxy C/T; the project's counted G and the table's effect T are positively aligned in that panel. An Ensembl GRCh38-to-GRCh37 positive-strand mapping was used for the proxy position.

For proxy T, the Estonian table reports beta = 0.01149, SE = 0.01871 and p = 0.539 (approximate OR 1.012, 95% CI 0.975–1.049; N = 205,165). Its direction is positive, but the interval includes a range of effects on either side of zero: it is **not a statistically clear replication**. This is an association at a correlated marker, not a direct rs817509 test; r² = 1 in a small external panel does not guarantee identical LD in the actual Estonian target, and exact allele and multiallelic representation in both studies need confirmation. It provides no FI/FFS result and cannot validate a 22-marker score.

Two other near pairs were assessed for perspective: rs17312861–rs1989903 (12,430 bp apart) had EUR r² = 0.091, and rs2155292–rs745631 (35,109 bp apart) had r² = 0.397. Neither is a strong proxy in that reference, even though both are nearby. These three pairs were selected post hoc based on proximity to the 33 selected insomnia rows, not through a genome-wide or prespecified proxy search. The selected Estonian results arise from a preprint, and sample and relative disjointness with every project's source has not been established. No causal locus, shared sleep–frailty mechanism, or score-level performance follows from these comparisons.

The [reproducible calculation](stage6_estonian_nearby_proxy_ld.py) and [receipt](stage6_estonian_nearby_proxy_ld.json) record the positions, alleles and correlations. The remote VCFs and locally downloaded indices were not independently checksum-verified. Mapping and allele records should be checked again before any confirmatory use. The [direct-marker review](stage6_estonian_sleep_replication_selected_variants.md) describes the table's selection and rs113851554.

## Expanded proximity screen

The [expanded screen](stage6_estonian_all_nearby_ld.py) evaluated every nonidentical selected Estonian insomnia lead within 1 Mb of either project list: 15 candidate–lead pairs, using 503-person European reference LD. In addition to rs817509–rs1093463, two pairs had r² > 0.8:

| Project marker | Estonian lead | EUR r² | Lead insomnia beta (SE) | p |
|---|---|---:|---:|---:|
| rs12536335 (FOXP2 region) | rs1989903 | 0.904 | 0.00960 (0.01314) | 0.465 |
| rs34290943 | rs73082337 | 0.924 | −0.03299 (0.01849) | 0.0743 |

These lead effects both have confidence intervals spanning zero; neither is a statistically clear association in the Estonian follow-up. Under the reference panel's allele coding, both directions correspond to the project's counted alleles, but this proxy orientation needs verification in the target study and does not establish the candidate's effect. In particular, nearby genes are positional labels, and these comparisons do not assign an effect to FOXP2. All screens are postselected and restricted to leads; the remaining pairs had r² ≤ 0.397. Correlation estimates use the external panel, and perfect or high correlation there does not establish target-sample LD, gene causality, full-score performance, or shared frailty risk. The full [pairwise receipt](stage6_estonian_all_nearby_ld.json) includes signed correlations and alleles.
