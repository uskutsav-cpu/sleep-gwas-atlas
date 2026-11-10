# Two candidate variants: effects, uncertainty and gene context

This is a detailed follow-up of the two strongest rows in the exploratory joint screen of MVP insomnia and Finnish hospital frailty (HFRS). The counted allele in every row is **G**, compared with A. Effects for continuous frailty traits are on their source-specific scales; do not compare their magnitudes directly with one another or with insomnia odds ratios. Confidence intervals below are approximate unadjusted Wald intervals calculated from reported beta and SE, not selection-adjusted intervals.

| Variant and source | Effect on G | Approximate 95% CI | Reported p |
|---|---:|---:|---:|
| **rs12536335**: FI meta-analysis | −0.0158 | −0.0225 to −0.0091 | 2.63 × 10⁻⁶ |
| FFS, allele convention inferred | −0.00579 | −0.00944 to −0.00214 | 0.0017 |
| Finnish HFRS | −0.0728 | −0.0975 to −0.0481 | 7.85 × 10⁻⁹ |
| UK Biobank FI component | −0.0189 | −0.0257 to −0.0121 | 4.6 × 10⁻⁸ |
| TwinGene FI component, inferred adapter | +0.0392 | +0.0107 to +0.0677 | 0.0070 |
| **rs1909650**: FI meta-analysis | +0.0169 | +0.0102 to +0.0236 | 6.32 × 10⁻⁷ |
| FFS, allele convention inferred | +0.00311 | −0.00058 to +0.00681 | 0.097 |
| Finnish HFRS | +0.0452 | +0.0202 to +0.0701 | 0.000386 |
| UK Biobank FI component | +0.0181 | +0.0113 to +0.0250 | 2.1 × 10⁻⁷ |
| TwinGene FI component, inferred adapter | −0.00285 | −0.0308 to +0.0251 | 0.842 |

The two component estimates form the FI meta-analysis and must not be counted as independent confirmations of it. For rs12536335 the component estimates differ in sign and each excludes zero under the inferred TwinGene adapter. That is a material warning, but the adapter is **not author-confirmed**; it could reflect coding or other source differences, and we cannot assert biological heterogeneity until coding is resolved. For rs1909650 the TwinGene interval is broad and includes either sign; it is not evidence of an opposing effect. Its FFS interval also includes zero, which neither establishes an FFS association nor proves no effect.

For clinical insomnia in MVP, the source reports G-versus-A odds ratios of approximately **0.969** for rs12536335 and **1.024** for rs1909650 (p=4.57 × 10⁻⁸ and 1.58 × 10⁻⁵ respectively). The source rounded the original odds ratios to three decimals and reports no SE for these rows, so intervals cannot be recovered reliably and are left blank. The first odds ratio was obtained by inverting the source A-versus-G value; it inherits source rounding. These are genotype associations with an insomnia case definition, not absolute probabilities or effects of changing sleep. UKB and FinnGen insomnia beta estimates and Wald intervals are preserved in `top_two_candidate_effects.csv`; their scales and definitions are source-specific.

## Nearby genes

- **rs12536335** overlaps the FOXP2 gene body. The allele-specific reference annotation includes intronic and noncoding transcript terms; overlap does not identify a causal variant or prove FOXP2 mediates either phenotype. The HFRS source also identifies this variant as intronic in FOXP2 and reports consistent HFRS direction in UK Biobank, but those results reuse cohorts and the dementia-excluded score does not exclude sleep or exhaustion.
- **rs1909650** lies 15,348 bp from the annotated pseudogene PCDH8P1. An Ensembl GRCh37 gene-body search found LINC01065 65,952 bp away and the closest protein-coding gene, OLFM4, 165,703 bp away; PCDH8 is 369,120 bp away. These are distances to gene bodies within a ±1 Mb search, not variant-to-gene functional links. Its reference rsID has multiple alternate alleles, so the exact G/A contrast must be verified in each study.

No fine-mapping, colocalization, single-cell effect or causal gene was established for either marker. Nearby genes can be informative leads but may not be targets of the associated variant. The HFRS includes broad codes capable of reflecting drowsiness or vital exhaustion; the original FI includes insomnia/tiredness. The candidate set was selected from prior cross-trait findings and the present inspection was post hoc. Participants may overlap across source GWAS, and ancestry and phenotype differences limit direct comparisons. Neither marker, alone or as part of the prototype, gives a validated individual risk score.

Source values and all calculations: `top_two_candidate_effects.csv` (SHA-256 `8eef39a9431135eb35deae332f715d27a304d8dcb9955e63030c80f69bbcfdf8`), drawn from the project’s harmonized candidate evidence, cohort lookup, raw MVP lookup and checksum-matched FinnGen HFRS file. Gene inventory: [Ensembl GRCh37 overlap query](https://grch37.rest.ensembl.org/overlap/region/human/13:52791895-54791895?feature=gene;content-type=application/json), downloaded response SHA-256 `ccfaa504a845e514e0c2399a4ed23b7a6c18d05c30b14e18241e618f31c5e293`. More detailed provenance is in `integrated_15_candidate_review.md`, `candidate_twingene_sensitivity.md` and `rs12536335_finngen_hfrs_followup.md`.
