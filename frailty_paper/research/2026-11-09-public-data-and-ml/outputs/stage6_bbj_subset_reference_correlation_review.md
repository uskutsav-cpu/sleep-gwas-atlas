# Reference coverage of the BBJ insomnia marker subsets

The earlier BioBank Japan (BBJ) insomnia follow-up used overlapping 14-, 15- and 17-marker subsets of the project's 22-marker equal-weight count. To assess how much the missing alleles vary in a reference panel, we recalculated each person's full 22-marker count and the three subset counts in 1000 Genomes Phase 3. These are reference genotypes, **not BBJ participant scores or outcome data**.

| Reference sample | n | Full vs 14-marker correlation | Full vs 15-marker correlation | Full vs 17-marker correlation |
|---|---:|---:|---:|---:|
| Japanese (JPT) | 104 | 0.970 | 0.996 | 0.999 |
| East Asian pooled (EAS) | 504 | 0.974 | 0.995 | 0.999 |
| European pooled (EUR) | 503 | 0.876 | 0.905 | 0.944 |

For the 14-marker subset the omitted count has SD 0.565 in JPT versus 1.467 in EUR; the full count has SD 2.319 and 3.037 respectively. The 15-marker subset adds rs4886140 to the 14-marker set, and the 17-marker subset adds rs2815764 and rs817509. In these small reference samples, omission affects the count's variation differently across groups. JPT is a subset of EAS, not an independent replication; neither reference group establishes BBJ genotype frequencies or sample LD. These populations are heterogeneous and not representative of everyone with those labels.

Correlation is invariant to constant shifts, so a high correlation does not mean that the partial and full raw scores have the same scale or value. It cannot establish that the scores have the same outcome association, discrimination, calibration or clinical usefulness. The marker subsets were chosen after inspecting data, not as prespecified alternative scores. For the 17-marker BBJ association sensitivity, rs2815764 and rs817509 also have low estimated case minor-allele dosage sums; reference-count correlation cannot fix the instability or lack of a validated score-level test. The [BBJ coverage analysis](stage6_bbj_insomnia_marker_coverage_sensitivity.md) and [case dosage audit](stage6_bbj_insomnia_candidate_review.md) describe those limitations.

Reproduce with `python3 work/bbj_subset_reference_correlation.py`. The [receipt](stage6_bbj_subset_reference_correlation.json) records input hashes, subsets and omitted markers. Each variant had one exact biallelic reference position, identifier and allele-pair match, and complete called diploid reference genotypes. The remotely indexed VCF was not verified by a full-file checksum. This check uses the equal-weight count only, not the separately oriented weighted 22-marker score or 15-marker prototypes.
