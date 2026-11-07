# Cohort allele-orientation sensitivity

The public Atkins cohort file and its README provide conflicting evidence about how to interpret the TwinGene effect at rs2396766. The README describes EA as the effect allele for both cohorts. With that interpretation, and aligning both results to the published meta-analysis effect allele A, the UK Biobank beta is +0.02157 and the TwinGene beta is +0.00527. The published article instead reports UK Biobank/TwinGene directions +−.

| Interpretation | TwinGene beta on A | Approximate pooled beta | Approximate heterogeneity p | Difference in cohort EAFs on the assumed common allele |
| --- | ---: | ---: | ---: | ---: |
| Both cohort effects use the stated EA | +0.00527 | +0.02070 | 0.268 | 0.0552 |
| Hypothetically flip TwinGene effect and frequency | −0.00527 | +0.02013 | 0.0680 | 0.00035 |
| Published meta-analysis table | Direction − | +0.020 (source extract +0.0201) | 0.07 | Not available in the table |

The hypothetical flip is consistent with the published directions, pooled beta and rounded heterogeneity p, and makes the reported frequencies closer. It also gives a closer approximate pooled estimate for rs1363103. These observations support investigating the orientation but **do not establish it**: the README says otherwise, frequencies can vary across cohorts, and the paper used METAL with genomic control. The simple inverse-variance and two-study calculations here are not an exact reproduction. Author clarification and original harmonization code are needed before choosing an adapter. The TwinGene estimate for rs2396766 is imprecise (SE 0.01431; source p 0.713). The published p=0.07 for heterogeneity is above 0.05 even though the paper's text calls heterogeneity significant.

A complete, gzip-validated scan of the 16,557,800-row cohort source confirmed that each of the four selected rsIDs appears exactly once and matches the saved extract. The selected rows, both model results, source checksums, verification receipt, and reproducible code are in this directory. An [unsent author question](AUTHOR_QUESTION_DRAFT.md) records the specific ambiguity. The data concern the same FI cohorts as the original study, include sleep-related components, and do not establish independent replication or a causal shared sleep–frailty variant. [Paper and table](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/).
