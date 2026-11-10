# Sensitivity of the 22-marker count to removing a region

I removed each of 18 descriptive marker groups in turn from the existing unweighted 22-marker shared-direction count, then compared the reduced and full counts in 503 EUR-labeled 1000 Genomes reference individuals. These are **correlations between two genotype sums**, not genetic correlations between traits or evidence of accuracy, risk, or clinical validity. The full and reduced scores overlap by construction, so high correlations are expected.

| Removed marker group (positional gene labels) | Markers removed | Pearson correlation of reduced with full count | Spearman rank correlation |
|---|---:|---:|---:|
| rs12536335, rs12671330, rs17312861 (FOXP2) | 3 | 0.898 | 0.893 |
| rs7731706, rs7734045 (TRPC7, AC009014.3) | 2 | 0.924 | 0.919 |
| rs9536347, rs1909650 (OLFM4, PCDH8P1) | 2 | 0.954 | 0.951 |
| rs6452792 (CTC-498M16.4) | 1 | 0.971 | 0.969 |
| rs2155292 (NCAM1) | 1 | 0.972 | 0.971 |

The FOXP2-labeled set had the greatest influence on the full count among these 18 groups in this reference panel. Its removal lowers the count's sample SD from 3.037 to 2.778 allele-count units; the group itself has SD 1.335. Because the removed markers covary with remaining markers, the change in variance cannot be attributed to this group by simply squaring its SD. These comparisons provide a way to describe redundancy and score sensitivity, not a reason to discard that region or a measure of its contribution to disease. Nearest-gene annotations do not establish causal genes.

## Method and boundaries

Groups are connected components of markers on the same chromosome no more than 1 Mb apart, a descriptive distance rule. This produced one three-marker group, two two-marker groups and 15 single-marker groups. It is not a definition of independent loci, and different choices can yield different groupings. The same fixed count was used for all deletions; no weights were refit. Pearson correlations use raw allele counts, and Spearman correlations use midranks to account for ties. The [complete group results](shared22_leave_region_reference.csv) contain all 18 comparisons.

The analysis matched the original score's assumed GRCh37 positions, rsIDs and allele pairs to biallelic 1000 Genomes phase 3 records; all 503 EUR-labeled samples had complete calls. It reproduced the previously reported full-count SD. Pooled reference subjects may be related and are not representative of every target population; results can change with ancestry, sample and LD. Some rsIDs have multiallelic representations whose compatibility with the original GWAS needs further checking. Individual genotypes and participant scores were not saved. The score file SHA-256 is `c419c1617dc7bf496f5213a4b2a8fae2de7a36b44984ba6d0ea0564dc3e40b16`; executable script: `work/shared22_leave_region_reference.py`.

This analysis cannot determine whether deletion improves or harms association or prediction. A prospective comparison of the fixed full count and prespecified reduced counts requires an independent dataset with measured outcomes; the user has specified free public data, and such an eligible participant-level dataset has not been established. Interpret alongside the [variance decomposition](shared22_count_variance_review.md) and [score guide](shared22_score_research_guide.md).
