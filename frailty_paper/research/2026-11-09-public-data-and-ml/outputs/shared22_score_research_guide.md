# A 22-marker shared-direction research score

**Reversal safeguard, 29 October 2026:** The count scorer now rejects reversed-allele dosage input by default. Supply counted-allele dosages, or use `--allow-verified-reversal` only after independently verifying all reversed rows are compatible diploid biallelic contrasts. The flag cannot establish this from two allele labels, and it does not validate the score.

This prototype supplies a single research sum for the user's requested shared-risk question. It is an **unweighted allele count**, not a validated genetic risk score or absolute probability. The [22-marker definition](unvalidated_22_shared_direction_score.csv) gives the specific counted and other alleles, assumed GRCh37 positions and positional gene labels. The accompanying [scoring script](shared_score_22_research.py) calculates the count using the input format and strict checks in [the scoring README](score_22_research_readme.md).

## Definition

The 22 variants are the post-selected candidates that had nominal p < 0.05 in both MVP clinical insomnia and FinnGen hospital-frailty (HFRS) GWAS. For each, the counted allele is the allele with a higher estimated MVP insomnia odds ratio; the estimated HFRS beta for the same allele is also positive. This agreement is a characteristic of this selected set and the observed estimates, not proof that the variant or allele causes both conditions. Nearby-gene labels are positional and do not establish effector genes.

For 22 verified biallelic diploid variants, the score is:

**Shared-direction allele count = sum of the dosages of the 22 listed counted alleles.**

Each allele has weight 1, so the theoretical range is 0–44 for diploid dosages. Imputed dosages can make the sum fractional. This equal-weight definition is a transparent demonstration, not a fitted model. It ignores differences in effect magnitude and LD, and it should not be interpreted as a calibrated measure of shared risk. The selection used HFRS, which is distinct from both FI and FFS. All 22 also have positive FI point estimates after orientation to the count allele; 21 have positive FFS point estimates, but many estimates are imprecise and the exception does not establish protection. The estimated directions are not evidence of a shared causal variant or an intervention effect.

The [eight-source evidence matrix](shared22_eight_source_evidence_matrix.csv) aligns each candidate to this score’s counted allele; its [review](shared22_eight_source_evidence_review.md) explains source overlap, nominal selection and gene-label limitations.

The earlier [22-candidate table](expanded_22_candidate_risk_gene_review.md) uses the earlier weighted-prototype allele orientation. Eight alleles are reversed relative to this count; use the [crosswalk](shared22_allele_orientation_crosswalk.csv) and the correspondingly aligned effect signs when comparing them. Matching alleles across these project files does not independently verify the original GWAS coding.

The [regional LD review](shared22_reference_ld_interpretation.md) shows five same-chromosome pairs within 1 Mb and aligns their phased correlations to the score's counted alleles. In particular, three markers carry FOXP2-region labels. Their presence does not establish three independent signals or a causal gene.

A [genotype-only variance decomposition](shared22_count_variance_review.md) finds that covariance of those five local pairs contributes 13.0% of this particular count's sample variance in 503 European-labeled reference individuals. It does not measure variance explained in either trait or a change in prediction.

## Technical check and interpretation

Using the named alleles in 503 EUR-labeled 1000 Genomes reference individuals, the unweighted count correlates with the existing weighted FI, FFS and insomnia sums at **0.977, 0.781 and 0.972**, respectively. These are correlations among deliberately related constructed sums, not correlations with actual traits, independent evidence, prediction accuracy or genome-wide genetic correlations. Candidate selection and allele orientation make high agreement unsurprising. No outcomes were available. [Reference-check receipt](shared22_reference_diagnostic.json).

The implementation was checked using synthetic direct and reversed allele dosages and verified against the expected allele count. A valid score requires every candidate in each sample, a matching verified **biallelic diploid** contrast and dose in [0,2]; the existing validation rejects missing rows, duplicate rows, mismatched alleles and nonfinite or out-of-range dosages. It does not ascertain the genome build, imputation quality, palindromic strand, exact source coding or population suitability. Fourteen rsIDs have additional alleles in broader reference data; the same rsID is insufficient for defining the source contrast. The original Song build is assumed, FI cohort effect conventions require confirmation, and the FFS A1 convention was checked at published lead loci but not independently confirmed for every candidate. Existing FI/FFS GWAS share UK Biobank data, and FinnGen HFRS may overlap the FinnGen insomnia follow-up. These limitations remain even when the software runs successfully.

For a research evaluation, freeze the selection and count, define one outcome and time horizon, ascertain discovery/selection–evaluation independence and genotype QC, then assess association, discrimination, calibration and incremental value over clinical predictors and existing genome-wide scores in an independent cohort. Because one count may perform differently for FI, FFS and insomnia, test these as separate prespecified outcomes; do not assume a joint phenotype. See [the validation design](independent_outcome_validation_cohort_feasibility.md). Until this is done, do not use the count for individual or clinical decisions.

## Running the example implementation

Keep `shared_score_22_research.py`, `score_22_research.py` and `unvalidated_22_shared_direction_score.csv` together and run:

```bash
python3 shared_score_22_research.py --dosages research_dosages.csv --output shared_direction_counts.csv
```

The dosage file is the same long-format `sample_id,variant,dosage_allele,other_allele,dosage` format described in the [README](score_22_research_readme.md); the output path must not exist. Do not put identifiable participant data into an unauthorized workspace. No participant-level scores were produced in this project.

A concrete three-allele example and the 22-marker reference allele inventory are in the [multiallelic representation audit](shared22_multiallelic_representation_audit.md). The software accepts preprocessed allele strings and cannot verify whether a source or sample is genuinely biallelic.
