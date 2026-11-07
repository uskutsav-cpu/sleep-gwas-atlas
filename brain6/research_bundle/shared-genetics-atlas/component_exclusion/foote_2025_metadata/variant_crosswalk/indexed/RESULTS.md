# Foote frailty-factor lookup for five project anchors

An indexed lookup in the published, harmonized GRCh38 genomic-SEM files found unique coordinate, rsID and allele matches for 33 of 35 anchor–factor combinations. rs1476535 was absent at the queried coordinate in the general and multimorbidity files. Absence cannot distinguish lack of source coverage, heterogeneity filtering or other exclusions and is not evidence of no association.

| Previously selected anchor | Published general-factor p | Interpretation |
|---|---:|---|
| rs2431108, Brain6 chr5 | 1.84×10⁻²² | Below the paper's 7.14×10⁻⁹ seven-factor GWAS threshold. |
| rs7703746, previously reported FFS lead | 1.36×10⁻¹⁷ | Below that same threshold. |
| rs62264764, Brain6 chr3 | 7.86×10⁻⁵ | Does not cross that threshold at this exact marker. |
| rs817509, chr3 frailty candidate | 7.47×10⁻⁵ | Does not cross that threshold at this exact marker. |
| rs1476535, Brain6 chr7 | Not returned | No conclusion about association. |

The source's factor-specific Q heterogeneity p for rs2431108 is 9.54×10⁻⁹, just above the release's removal threshold of 7.14×10⁻⁹. Its retention therefore should not be interpreted as evidence of negligible heterogeneity. The other six factors and all source fields are in `all_anchors.tsv`; treating a single non-significant marker as proof that its region has no factor association would also be incorrect. This table is a **post-selected descriptive lookup**, not a new significance family or independent validation.

The general factor includes insomnia among its 30 inputs; several contributing studies overlap the project sources. Its association at rs2431108 does not establish sleep-independent frailty, cross-trait pleiotropy or a common causal allele across insomnia, ADHD and FI/FFS. The latent factors differ from FI and FFS; their files exclude some significant Q-heterogeneity signals. The lookups used remotely indexed data and verified coordinate/allele matches but **did not verify the complete harmonized-file checksum**. The harmonized build-38 release is distinct from the original build-37 files. The original GRCh37 general-factor file was subsequently fully read and passed its publisher checksum; see the verification update below. The [lookup protocol](PROTOCOL.md), seven receipts and saved harmonized metadata preserve the precise limitations and reproduction details. Source: [Foote et al. 2025](https://www.nature.com/articles/s41588-025-02269-0).

An [original-release row crosscheck](ORIGINAL_RANGE_RESULTS.md) confirms four of these general-factor records after accounting for effect-allele flips. It does not verify either complete file or make the two releases independent.

The [complete original general-factor read](../VERIFIED_GENERAL_FACTOR.md) passed its publisher MD5 and matched all selected indexed records after allele alignment; the full harmonized file and other six factor originals have not been checksum-verified.
