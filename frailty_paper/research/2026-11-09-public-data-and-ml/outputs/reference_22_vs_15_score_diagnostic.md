# Two alternative sparse genetic-score prototypes: technical comparison

I constructed three **separate** 22-variant research sums using the already aligned FI, FFS and UK Biobank insomnia marginal coefficients in `expanded_22_candidate_risk_gene_table.csv`. The coefficients, counted alleles and positional labels are in `unvalidated_22_variant_score_prototype.csv`. These are post-selected, uncalibrated illustrations, not clinical scores or probabilities. Do not add the FI, FFS and insomnia sums: they represent different source outcomes and scales.

In 503 European-labeled 1000 Genomes reference individuals, I computed those three sums and the previously prepared 15-marker sums. The lists overlap at **4 variants**; they use different selection criteria, so the 22-marker set is not an extension that contains all 15 original markers. Correlations between the two alternative sums for the same trait were:

| Trait-specific sums | Pearson r: original 15 vs new 22 |
|---|---:|
| FI | 0.346 |
| FFS | 0.198 |
| UKB insomnia | 0.293 |

The substantially different numerical patterns demonstrate that these **two particular post-selected candidate lists and weights** produce scores with limited linear agreement in this reference sample. They do not show that one score predicts better, that personal disease risk changed, or that the lists exhaust relevant variants. Correlations are influenced by variant selection, weights, allele frequencies, LD and population structure. No FI, FFS, insomnia or future outcome data were available for these individuals.

For completeness, the new 22-marker sums correlate with one another at **0.783** (FI–FFS), **0.955** (FI–insomnia), and **0.797** (FFS–insomnia). These are correlations of constructed sums, not genome-wide genetic correlations or clinical phenotype correlations. They arise partly because the candidates were selected for cross-trait evidence and most aligned marginal effects have the same estimated sign. They must not be presented as independent confirmation of shared biology.

## Technical checks and limitations

All **33 unique variants** in the union of the two lists matched the named rsID, assumed GRCh37 position and the specified biallelic contrast in the public 1000 Genomes phase 3 reference VCF, with complete calls for all 503 individuals. For variants whose two lists used opposite counted alleles, dosages were reversed correctly when computing each separate sum. The previous 15-marker cross-score correlation matrix was reproduced to numerical precision. The recorded reference matching does not resolve the original Song table's unverified build, multiallelic representation in original GWAS, the inferred FFS effect convention, or the TwinGene FI coefficient issue. Fourteen of the 22 rsIDs have more than two listed alleles in broader Ensembl reference data. The same 503 samples include several EUR-labeled subgroups and may include relatives; no ancestry-PC adjustment or calibrated uncertainty was applied. Subgroup descriptions in the JSON are not evidence of differences between populations.

Both sets were selected after looking at association results, and their GWAS can share participants. FI contains sleep/tiredness-related items, FFS includes exhaustion, and the follow-up HFRS can contain related broad diagnosis codes. Marginal weights were not jointly optimized for LD, externally tuned, assessed for missing genotypes, or tested against observed outcomes. Neither set demonstrates discriminative ability, calibration, incremental value over other predictors, clinical usefulness or causal association. The next meaningful evaluation requires an endpoint, fixed model, separate development/tuning and genuinely independent genotyped participants with measured outcomes, ideally including FI with sleep and fatigue removed as a sensitivity.

`reference_22_vs_15_score_diagnostic.json` contains the complete aggregate correlations, subgroup summaries, variant match checks and score-file hashes. The new score CSV SHA-256 is `9db36af13ee01fd721ffba239bce70f164b65cefce586ac1b2a855b606942224`. Script `work/expanded22_score_check.py` fetches the indexed reference VCF; no person-level genotypes or scores were saved. The earlier 15-marker prototype is unchanged.
