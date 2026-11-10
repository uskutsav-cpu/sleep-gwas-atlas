# Exploratory SMAD5 skeletal-muscle splicing colocalization

**Subsequent QC finding:** The [molecular credible-set audit](stage3_skeletal_molecular_finemapping_review.md) found incomplete credible-set coverage, including zero of 18 SMAD5 L2 variants in the filtered FinnGen input, and two SMAD5 molecular signals. The single-signal H3/H4 estimates below are not a reliable conclusion of shared or distinct signals.

I tested one SMAD5 skeletal-muscle splice junction in the rs7734045 region against Mak FI, Mak modified Fried physical frailty, and FinnGen insomnia. The source-selected event is `5:136132962:136153592:clu_31306_+` (GRCh38) from eQTL Catalogue release-7 GTEx muscle leafcutter dataset QTD000285. It is a specific splicing phenotype, not total SMAD5 expression. At the candidate marker, the Catalogue reports p = 8.92 × 10⁻⁵⁹ for this event. The GTEx portal v10 identifies the same junction boundaries under a different cluster identifier; exact cross-release cluster membership has not been established.

| Pair | H1: splicing association only | H3: distinct signals | H4: one shared signal | H3 + H4 |
|---|---:|---:|---:|---:|
| Splicing vs Mak FI | 0.0027 | 0.9966 | 0.00072 | 0.9973 |
| Splicing vs Mak physical phenotype | 0.9674 | 0.0202 | 0.0124 | 0.0326 |
| Splicing vs FinnGen insomnia | 0.8674 | 0.0431 | 0.0895 | 0.1326 |

These are exploratory, model-dependent probabilities under illustrative p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵, with a splicing effect prior SD of 0.15. For FI the single-signal model favors distinct signals; for physical frailty and insomnia it mostly favors a splicing-only association. The H3 + H4 values for the latter two are small, so conditional H4 alone would be misleading. This does not prove that no signal is shared, that this splice event is causally linked to disease, or that results apply to all SMAD5 transcripts or cell types.

The [sensitivity results](stage3_smad5_splice_coloc_abf_sensitivity.csv) vary p12 from 10⁻⁶ to 5 × 10⁻⁵ and splicing-prior SD from 0.075 to 0.30. H4 ranges from 0.000072–0.00371 for FI, 0.00123–0.0592 for physical frailty, and 0.00941–0.331 for insomnia. The event was selected after seeing prior QTL evidence, so these are not confirmatory or multiple-testing-adjusted analyses.

## Evidence and limitations

The `.cc` QTL file retains selected molecular traits corresponding to fine-mapped signals rather than all splice events; per the [Catalogue release notes](https://www.ebi.ac.uk/eqtl/Release_notes/), this filtering is designed to retain summary statistics for selected events, and the file is not comprehensive. The [extraction receipt](stage3_catalogue_smad5_splicing_extract_receipt.json) records the event, source and hashes. The Catalogue's [methods](https://www.ebi.ac.uk/eqtl/Methods/) describe inverse-normalization of LeafCutter junction usage. A standard-normal prior scale is a working assumption, not verified participant-level phenotype variance.

To avoid another scan of each whole GWAS, the splicing variants were intersected with the existing SMIM32 eQTL–GWAS matched SNP sets. This retains 2,404 variants for each Mak trait and 2,212 for FinnGen, including rs7734045, but further restricts coverage and can omit sQTL variants that were absent from the expression input. Alleles, exact positions and rsIDs were checked; assumed build conversion and complex alleles still need verification. The 1 Mb window is candidate-centered rather than the complete gene cis window. QTL and GWAS donor overlap is unknown.

FI and physical-frailty prior scales use the data-derived estimates explained in the [SMIM32 FI/FP analysis](stage3_smim32_skeletal_mak_coloc.md); insomnia uses log-odds prior SD 0.20. This custom single-signal ABF implementation was algebraically checked on a small example and subsequently checked against the unmodified coloc 5.2.3 `approx.bf.estimates` and `combine.abf` functions on the same inputs; the full package workflow was not run. Multiple causal signals, LD/population mismatch, missing variants, low GWAS power and selection can change conclusions. The [machine-readable receipt](stage3_smad5_splice_coloc_receipt.json) and local script `work/stage3_smad5_splice_coloc.py` record inputs and exact model settings. This analysis is not fine-mapping, multi-signal colocalization, independent replication or an individual risk estimate.

A later [MAF sensitivity analysis](stage3_smad5_low_frequency_sensitivity.md) showed that all 18 L2 variants are present in the full FinnGen source at MAF below 1%; lowering the filter includes 14. The sensitivity does not resolve the multiple-signal problem.

[Reference-formula validation](stage3_coloc_reference_formula_review.md) checks arithmetic only; credible-set coverage and single-signal model limitations still govern interpretation.
