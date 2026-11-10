# Clinical insomnia follow-up for all 22 candidates

All 22 selected candidate variants had a unique match by rsID, mapped GRCh38 position and allele pair in the complete FinnGen R13 clinical insomnia GWAS. **Nine** have nominal FinnGen insomnia p < 0.05, and all 22 have the same estimated effect direction as MVP insomnia when aligned to the same allele. Weak estimates can change sign; the 22 matching signs alone are not proof of an effect or a common causal mechanism. The 22 candidates were selected based on prior cross-trait, MVP insomnia and FinnGen hospital-frailty results, so the screen is post-selected.

Six of the 22 have matching point-estimate directions and unadjusted p < 0.05 in **all six** reviewed sources: FI, FFS, UKB insomnia, MVP insomnia, FinnGen insomnia, and FinnGen hospital frailty (HFRS). This extends the [earlier six-candidate follow-up](six_candidate_cross_trait_followup.md), in which two of the six had that pattern. The four additional variants below were in the broader 22, not in that earlier six-marker subset.

| Variant | Counted allele | Nearest gene-body label | FinnGen insomnia beta (SE) | FinnGen insomnia p | FI p | FFS p | HFRS p |
|---|---|---|---:|---:|---:|---:|---:|
| rs12536335 | G | FOXP2 | −0.0183 (0.00670) | 0.00625 | 2.63×10⁻⁶ | 0.0017 | 7.85×10⁻⁹ |
| rs817509 | G | RP11-384F7.1 | +0.0249 (0.00832) | 0.00272 | 0.0223 | 0.0099 | 1.21×10⁻⁷ |
| rs2838788 | C | ADARB1 | −0.0246 (0.00687) | 0.000345 | 0.000791 | 0.0040 | 0.000275 |
| rs34290943 | C | LAMB2 | −0.0296 (0.00908) | 0.00110 | 2.60×10⁻⁵ | 0.00011 | 0.00401 |
| rs12671330 | C | FOXP2 | +0.0237 (0.00668) | 0.000384 | 2.03×10⁻⁹ | 0.0017 | 0.00718 |
| rs7734045 | C | AC009014.3 | +0.0153 (0.00726) | 0.0352 | 0.00545 | 0.015 | 0.000325 |

FinnGen insomnia beta is a log-odds coefficient per counted allele. See [the full table](22_finngen_insomnia_followup.csv) for allele pairs, SEs, approximate intervals and source p-values; [the prior 22-variant table](expanded_22_candidate_risk_gene_table.csv) contains all six-source inputs except the FinnGen insomnia column. All p-values here are unadjusted, selection and multiple analysis rounds matter, and not every source is independent: FI and FFS use UK Biobank, while FinnGen insomnia R13 and FinnGen HFRS R12 may share people. Neither nominal significance in six analyses nor direction agreement establishes one shared causal variant, a causal gene, effect differences or an independently validated risk score.

The two FOXP2-labeled markers are 290,884 bp apart and show reference-panel r² approximately 0.106 for the specifically tested contrasts in 503 EUR-labeled 1000 Genomes individuals. They should not be counted as independent loci solely from their distinct rsIDs; nor does the observed LD identify a shared causal variant. Gene names are nearest gene-body labels, not proof of the effector gene; for rs34290943 the source paper discussed DALRD3 in a broader analysis while the closest-body mapping here is LAMB2. The FI UKB and TwinGene component betas for rs12536335 differ in sign under the inferred TwinGene adapter. Original source build, multiallelic coding and phenotype definitions continue to matter. FI and FFS contain sleep/fatigue-related items; HFRS includes broad codes that can reflect drowsiness or exhaustion.

## Verification

We scanned all 21,327,010 records of the checksum-verified FinnGen R13 insomnia file, matching coordinates, rsID and allele pair and aligning ALT beta to each candidate's counted allele. The six previously checked candidates reproduce the prior aligned betas and p-values after allele conversion. [Receipt and counts](22_finngen_insomnia_followup.json); script `work/inspect_22_finngen_insomnia.py`. No new model, independent outcome-based prediction, or genotype data were generated.
