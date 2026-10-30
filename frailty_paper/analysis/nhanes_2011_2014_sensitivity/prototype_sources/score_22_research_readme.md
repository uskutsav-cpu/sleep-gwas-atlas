# Research-only 22-marker scoring tool

[score_22_research.py](score_22_research.py) computes three separate raw additive sums from [the 22-variant prototype weights](unvalidated_22_variant_score_prototype.csv): one using FI weights, one FFS weights and one UK Biobank insomnia weights. It is an implementation of an **unvalidated, post-selected research prototype**, not a clinical or personal risk calculator. Its outputs are not probabilities, have no validated threshold or calibration, and cannot tell anyone whether they will become frail or develop insomnia. Never add the three scores together: the phenotypes and scales differ. No participant outcomes were used to train or validate it.

## Input and use

A long-format UTF-8 CSV must have one row for each sample and each of the 22 variants:

```csv
sample_id,variant,dosage_allele,other_allele,dosage
example,rs12536335,G,A,1.1
```

The example is only one row and is deliberately insufficient to score. `dosage` is the expected **diploid copy count** of `dosage_allele` between 0 and 2, including fractional imputed dosages. Allele strings must match the same two-allele contrast and strand as the prototype; when the input dosage allele is the model's other allele, the tool uses 2 minus dosage. This reversal is valid only for a verified biallelic diploid contrast. A multiallelic record, strand ambiguity, mismatched genome build, sex-chromosome ploidy or genotype-to-variant mismatch cannot be solved by matching rsID and two letters. The tool does not inspect VCF files, infer assembly or handle missing genotypes; resolve these using appropriate source-specific QC **before** scoring. Fourteen of the 22 rsIDs have multiple alleles in broader reference records, and the Song coordinate build remains an assumption.

The two-label input cannot show whether a third allele is present, so the script cannot enforce the biallelic condition. At an A/G/T site an A/T genotype has A dosage 1 and G dosage 0, but `2 − A dosage` incorrectly gives 1. Resolve the per-allele dosage and GWAS effect contrast before scoring; see the [shared-marker audit](score_15_22_shared_marker_crosswalk.md).

As of 29 October 2026, reversed dosage input fails by default. Prefer dosage of the model's counted allele. To restore the prior complement behavior, pass `--allow-verified-reversal` only after independently verifying every reversed row is a compatible diploid biallelic contrast. The flag records an analyst assertion; it cannot verify the record, source GWAS contrast, or additional alleles. This changes earlier command behavior and does not validate the scores.

In the [22-marker FFS alignment audit](score22_ffs_alignment_audit.csv), all 22 cached FFS rows match the prototype's reported positions and allele pairs, but **all 22** have source A1/A2 reversed relative to the prototype's counted/other alleles. All stored FFS weights equal the negated source BETA, assuming A1 is the effect allele. Fourteen of these rsIDs have more than two alleles in the project's reference annotation. This is an arithmetic check, not proof that a third allele occurs in the actual source or target; verify each source contrast, target genotype representation, build, and ploidy before using a reversed weight. [Prior source-table checks](ffs_effect_allele_source_consistency_audit.md) support the A1 convention at published FFS leads, but do not resolve candidate-specific multiallelic coding. Reproduce with `python3 work/audit_22_ffs_allele_swaps.py` from the project root.

Run using a local Python 3 installation:

```bash
python3 score_22_research.py --dosages research_dosages.csv --output research_scores.csv
```

By default the script reads the weight file next to itself; `--weights` specifies an alternative. The output must not already exist. Every sample must have exactly one valid row for each model variant. Duplicate or empty CSV column names, malformed rows, duplicate variants, unmatched alleles, unknown variants, missing rows, nonfinite or out-of-range dosages fail the whole run; no partial score file is written for these validation failures. The output has one row per sample and raw sums. No normalization is applied: means and SDs would have to be defined from an appropriate training/reference population, never inferred as absolute disease risk. Use deidentified sample IDs and follow applicable data access and security rules.

## Validation status

The script's allele reversal, dosage validation and arithmetic were checked on synthetic data. A separate genotype-only reference check ran all 503 European-labeled 1000 Genomes samples in memory, reversed the input allele convention for 11 of 22 variants, and reproduced the earlier three-score correlation matrix to within 6 × 10⁻¹⁶. The [receipt](score22_reference_implementation_check.json) records the check; no individual genotypes or scores were retained. These checks verify implementation, not clinical utility. The chosen weights are marginal, post-selected and not LD-adjusted; FI includes sleep-related items, the FFS effect-allele convention is supported by a source-table audit rather than author confirmation for every record, and source overlap and coding require review. The 22-marker scores have not been calibrated or evaluated against independent FI, FFS, insomnia or future disease outcomes. See [the risk assessment](mentor_risk_assessment_updated.md), [FFS allele audit](ffs_effect_allele_source_consistency_audit.md), and [validation design](independent_outcome_validation_cohort_feasibility.md).

A concrete three-allele example and the 22-marker reference allele inventory are in the [multiallelic representation audit](shared22_multiallelic_representation_audit.md). The software accepts preprocessed allele strings and cannot verify whether a source or sample is genuinely biallelic.
