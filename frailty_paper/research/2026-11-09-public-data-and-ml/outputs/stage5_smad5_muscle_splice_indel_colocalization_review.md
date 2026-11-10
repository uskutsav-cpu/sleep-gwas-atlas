# SMAD5 skeletal-muscle splice event: indel coverage sensitivity

For the selected eQTL Catalogue QTD000285 LeafCutter event `5:136132962:136153592:clu_31306_+`, the frozen molecular fine-mapping reports **two signals**. We added eligible exact-allele credible-set indels to the earlier single-signal analysis, while leaving its other filtering unchanged. These additions do not make a single-signal model appropriate for a two-signal locus.

| Trait | Indels added | L1 variants before → after | L2 variants before → after | L1 covered PIP sum before → after | L2 covered PIP sum before → after |
|---|---:|---:|---:|---:|---:|
| Mak FI | 8 | 37 → 42 of 52 | 14 → 17 of 18 | 0.529 → 0.663 | 0.755 → 0.877 |
| Mak physical frailty | 8 | 37 → 42 of 52 | 14 → 17 of 18 | 0.529 → 0.663 | 0.755 → 0.877 |
| FinnGen insomnia | 5 | 38 → 43 of 52 | 0 → 0 of 18 | 0.531 → 0.665 | 0 → 0 |

L1 and L2 total PIP sums are 0.950 and 0.956, respectively. The three L2 indels have FinnGen alternate-allele frequencies outside the earlier 1%–99% filter and were excluded from the baseline; this is not evidence that they have no effect. The [coverage table](stage5_smad5_splice_indel_coloc_coverage.csv) and [exclusions](stage5_smad5_splice_indel_exclusions.csv) preserve exact counts and reasons. Only eight selected molecular-set indels were added; other excluded variants and indels were not recovered. QTL imputation r² was unavailable for these variants and not filtered, consistent with the original analysis.

At the original baseline priors (p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, splice-prior SD 0.15 and original trait-specific scales), H3/H4 are 0.99655/0.00078 for Mak FI, 0.02028/0.01207 for Mak physical frailty, and 0.04336/0.08982 for FinnGen insomnia. H1 is 0.00267, 0.96765 and 0.86681 respectively; H0/H2 supply the remainder. The numbers change little from the earlier analysis, but **the single-signal H3/H4 values are not reliable conclusions about sharing** because the molecular phenotype has multiple signals and FinnGen covers none of L2. The [sensitivity table](stage5_smad5_splice_indel_coloc_sensitivity.csv) varies priors; this is not a multi-signal colocalization.

Original baseline values were recomputed before adding indels, and the updated values agree with the `coloc` v5.2.3 formulas to less than 4 × 10⁻¹⁶ ([check](stage5_smad5_splice_indel_coloc_formula_check.csv)). This verifies arithmetic only. The event was selected after observing QTL evidence; source build and allele conventions are inherited assumptions; GWAS and model-training sample overlap, LD, ancestry, effect scales and missing variants constrain interpretation. A multi-signal method using appropriate ancestry-matched LD and well-covered variants is the next analysis, not an automatic conclusion from these posteriors. Full AlphaGenome splice predictions and independent allele-specific validation remain missing.

A later [public-reference LD check](stage5_smad5_indel_reference_ld_review.md) found seven of the eight credible-set indels in the broader 1000 Genomes Phase 3 VCF. This is useful preparation, but the eighth exact allele and a complete, ancestry-appropriate locus LD matrix remain unresolved; the single-signal limitation above remains.

A provisional [European-reference locus LD matrix](stage5_smad5_muscle_eur_reference_ld_readme.md) is now available for 2,345 SNPs and indels. Its molecular-set coverage remains incomplete, and its 503-person reference is rank-deficient; no multi-signal colocalization has yet been performed.
