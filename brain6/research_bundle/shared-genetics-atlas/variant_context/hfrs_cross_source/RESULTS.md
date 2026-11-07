# Hospital frailty at selected Brain6 variants: UK Biobank and FinnGen

Among 425 previously selected sleep–ADHD variants, 420 matched both the original UK Biobank hospital frailty risk score (HFRS) lookup and the Finnish HFRS lookup using the recorded GRCh38 position and alleles. **46 SNPs passed a 425-slot post-hoc max-p benchmark in both HFRS sources:** 19 in the chromosome 3 region and 27 in chromosome 7. These are correlated variant records in two selected regions, not 46 independent discoveries. Their aligned effects had the same sign across sources; across all 420 matched variants, 405 had the same sign. This descriptive count is dominated by correlated markers and is not a formal test of directional replication.

| Selected variant | UKB HFRS p | Finnish HFRS p | max-p × 425 | ALT direction across HFRS sources |
| --- | ---: | ---: | ---: | --- |
| rs62264764, chr3 | 0.000126925 | 4.5221×10⁻⁷ | 0.05394 | Both negative |
| rs62264767, chr3 | 0.000116018 | 4.62594×10⁻⁷ | 0.04931 | Both negative |
| rs1476535, chr7 | 3.93814×10⁻⁹ | 1.45603×10⁻⁸ | 0.00000619 | Both positive |
| rs2431108, chr5 | 0.00441979 | 0.00110375 | 1 | Both positive |
| rs77960, chr5 | 0.00476056 | 0.000967275 | 1 | Both positive |
| rs4916723, chr5 | 0.12052 | 0.0091939 | 1 | Opposite estimated signs |

The chromosome 3 candidates rs62264764 and rs62264767 sit on opposite sides of a threshold here, but their close values and [strong reference LD](../../replication/rs62264764_mvp/RESULTS.md) do not justify concluding that their underlying effects differ. The chromosome 7 result is stronger in this selected comparison, with rs2894699 having the smallest max-p (4.98×10⁻⁹). The opposite signs for rs4916723 are noisy estimates, not a formal test of effect heterogeneity. This comparison assesses HFRS associations only; it does **not** establish a replicated three-trait sleep–ADHD–frailty locus, a common causal variant, gene, or an AlphaGenome effect. The variant family and this contrast were selected after earlier results, so the numerical threshold is not a project-wide error guarantee. No threshold crossing is evidence of absence.

The max-p test concerns association in both sources if the source marginal p-values are valid; it does not require the source tests to be independent and does not test whether the effects have the same direction. UKB and FinnGen differ in ancestry, hospital ascertainment, phenotype distributions and model details; sample independence and relatedness have not been established. Raw effect sizes are not pooled or tested for equality. The UKB source reports N=407,643 in its rows while its paper and release description state 407,463. The previously mapped builds and allele conventions, source-specific GWAS calibration and variant frequencies have not been independently authenticated. The original local LAVA QC failures remain unresolved.

## Reproducibility

The [protocol](PROTOCOL.md), [script](compare.py), [all matched rows](results.tsv) and [summary and input hashes](summary.json) document this secondary join of the [UKB candidate lookup](../brain6_ukb_hfrs/RESULTS.md) and [Finnish candidate lookup](../brain6_mvp_finngen_hfrs/RESULTS.md). All 425 candidate slots remain in the adjustment, including five unavailable Finnish records. The match requires identical recorded rsID, region and GRCh38 position/allele string; this checks consistency between prior crosswalks, not independent verification of them. An additional calculation independently recomputed each max-p, each adjustment, all direction flags and the 46-marker count from the output. These arithmetic checks do not validate cohort models or causal interpretation.
