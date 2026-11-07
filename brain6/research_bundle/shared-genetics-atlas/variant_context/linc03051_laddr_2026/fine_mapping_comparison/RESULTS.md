# Molecular leads versus Brain6 chromosome 3 fine mapping

All six selected variants matched the archived Brain6 insomnia and ADHD regional inputs by rsID, GRCh37 position and allele pair. In the **archived single-signal model**, the four published depression/RNA lead markers have much smaller conditional disease PIPs than the two highlighted sleep–ADHD markers:

| Marker | Role | Insomnia conditional PIP | ADHD conditional PIP |
|---|---|---:|---:|
| rs1585084 | Depression–RNA top variant, hippocampus | 0.00000037 | 0.0000376 |
| rs9822630 | Depression–RNA top variant, anterior cingulate | 0.00000081 | 0.0000112 |
| rs13073833 | RNA lead, hippocampus feature | 0.00000081 | 0.0000112 |
| rs4597702 | RNA lead, anterior-cingulate feature | 0.00000086 | 0.0000123 |
| rs1456193 | Published insomnia xTWAS GWAS lead | 0.0615 | 0.00229 |
| rs62264764 | AlphaGenome-prioritized candidate | 0.0402 | 0.0312 |

These PIPs are conditional on each trait's archived model and selected regional variant support, not direct probabilities that a molecular feature or gene is causal. The four molecular/depression lead variants can tag other variants, and the model allows one signal per trait; it cannot exclude shared secondary signals. The earlier small European [reference LD check](../finngen_supplement/reference_ld/RESULTS.md) found moderate LD between depression leads and rs62264764. Neither positional differences nor small PIPs prove distinct disease and molecular causal variants. Conversely, a high molecular–insomnia coloc posterior from a different RNA phenotype and data source does not establish ADHD mediation.

The [full comparison](comparison.tsv) preserves effect-allele beta, SE and approximate two-sided Wald p-values; values are not independent replication or multiple-testing-adjusted results. [Code](compare.py) verifies all source inputs and the stored PIP file against the original pilot hashes and matches alleles; [receipt](receipt.json) records them. This is a post-selected descriptive comparison, not a new fine-mapping analysis. The original LAVA quality gates remain failed.
