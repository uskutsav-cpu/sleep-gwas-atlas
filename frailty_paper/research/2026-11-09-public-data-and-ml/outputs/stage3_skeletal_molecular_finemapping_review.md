# Muscle QTL credible sets: implications for the first colocalizations

The eQTL Catalogue release-7 SuSiE outputs report **one** skeletal-muscle credible set for SMIM32 expression and **two** for the selected SMAD5 splice event. These are molecular-QTL fine-mapping results, not GWAS credible sets or colocalization. The two splice credible sets indicate that the event's molecular architecture may contain multiple signals, making our earlier single-signal splicing colocalizations unsuitable for a definitive shared-versus-distinct conclusion.

| Molecular phenotype | Credible sets | Unique variants per set | Minimum within-set r² reported by source |
|---|---:|---:|---:|
| SMIM32 expression | 1 | 98 | 0.808 |
| SMAD5 selected splice junction | 2 | 52 and 18 | 0.897 and 0.565 |

The exact rs7734045 `chr5_136202653_A_C` variant was not listed in any of these three credible sets. That is not proof that its posterior inclusion probability is zero, that it has no molecular association, or that another correlated causal variant is unrelated to the GWAS. Its GTEx and Catalogue marginal QTL associations remain associations. SuSiE credible-set membership, lead variants, priors and purity are model- and dataset-dependent; gene causality is not established.

## Coverage of the previous single-signal comparisons

The molecular credible sets were computed over approximately 2 Mb. All of their listed variants lie within the prior 1 Mb candidate-centered query, but that does **not** mean the prior query covers all variants tested in each molecular model. Comparing exact GRCh38 chromosome/position/ref/alt IDs with the earlier colocalization inputs yields:

| Molecular credible set | Mak FI/FP input | FinnGen insomnia input |
|---|---:|---:|
| SMIM32 expression | 75 / 98 | 74 / 98 |
| SMAD5 splice signal L1 | 37 / 52 | 38 / 52 |
| SMAD5 splice signal L2 | 14 / 18 | 0 / 18 |

The FinnGen zero is **zero exact matches in the previously filtered colocalization input**, not proof of no coverage in the original FinnGen release or zero LD with the set. Earlier filters removed rare and palindromic SNPs, indels, duplicate IDs and variants missing from the SMIM32 expression intersection. The specific reasons have not yet been decomposed. Missingness of a complete molecular credible set, and multiple splice signals, are material weaknesses of those tests. In particular, the earlier H3/H4 values **must remain exploratory and must not be treated as evidence that the molecular and GWAS signals are distinct or shared**. More complete variant harmonization and a multi-signal approach with appropriate LD and all credible-set coverage are needed. [Exact overlap table](stage3_qtl_credible_set_coloc_support.csv).

## Source integrity and duplicate records

The [summary](stage3_skeletal_qtl_credible_set_summary.csv) and [source rows](stage3_skeletal_qtl_credible_set_variants.csv) use the Catalogue's QTD000281 expression and QTD000285 leafcutter credible-set releases. The Catalogue can duplicate a variant across rows when it has multiple rsIDs. Here the raw row counts exceeded reported credible-set sizes by 4, 4 and 1; collapsing exact variant IDs recovered the published sizes of 98, 52 and 18 and posterior sums near 0.95. Duplicated rsIDs must not be counted as independent variants or added twice. Source compressed files passed complete gzip reading and have local SHA-256 hashes recorded in [the receipt](stage3_skeletal_qtl_credible_set_receipt.json); provider checksums were not available in the inspected directory.

Sources and methods: [eQTL Catalogue release notes](https://www.ebi.ac.uk/eqtl/Release_notes/) and [data-access guidance](https://www.ebi.ac.uk/eqtl/Data_access/). The previous [expression](stage3_smim32_skeletal_finngen_coloc.md) and [splicing](stage3_smad5_splice_coloc.md) colocalization reports describe their initial model and priors. This fine-mapping check changes the interpretation of those reports, not the raw posterior values.

A subsequent [coverage audit](stage3_smad5_L2_finngen_coverage_audit.csv) found all 18 L2 variants at exact chromosome, position, reference and alternate in the full FinnGen source, with minor-allele frequencies approximately 0.91%–0.95%; all fell below the original 1% filter. Three are indels and one is palindromic, with overlapping exclusion flags. Presence in the file does not establish sufficient imputation quality, power, or eligibility for colocalization. A separately labeled [lower-frequency sensitivity](stage3_smad5_low_frequency_sensitivity.md) includes 14 of the 18 variants; its single-signal conclusions remain limited.
