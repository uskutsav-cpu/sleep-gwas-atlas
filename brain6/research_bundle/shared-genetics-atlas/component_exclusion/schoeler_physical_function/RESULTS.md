# Physical-function composite: sleep and frailty genetic correlations

The anonymously available [GCST90565837](SOURCE.md) baseline composite measures lung function, grip strength and height; higher values indicate better physical function. It contains no explicit sleep or exhaustion item, but it is not a frailty instrument with items removed. Using unconstrained LDSC, European LD scores, and the existing common-support inputs, we obtained the following exploratory estimates:

| Comparison with physical function | Genetic correlation (SE) | Logged p |
|---|---:|---:|
| Clinical insomnia (MVP) | −0.0900 (0.0274) | 0.0010 |
| Frailty index (FI) | −0.1909 (0.0232) | 1.6631 × 10⁻¹⁶ |
| Fried frailty score (FFS) | −0.4123 (0.0189) | 2.7126 × 10⁻¹⁰⁵ |
| Hospital frailty risk score (HFRS) | −0.0979 (0.0189) | 2.1272 × 10⁻⁷ |

Each passes the post-hoc four-comparison threshold p < 0.0125; this does not account for the broader project-wide search. These describe genome-wide shared association patterns, not causal effects or individual risk predictions. A negative estimate is consistent with alleles associated with more insomnia or frailty tending to align with those associated with lower physical-function scores. In particular, these results cannot isolate the influence of removing exhaustion: the phenotype, components and cohorts differ, and no matched score-deletion GWAS was available. A [post-hoc conservative contrast](conservative_contrasts/RESULTS.md) uses the sum of reported SEs as a covariance-independent approximation; no paired block-jackknife comparison was saved, and phenotype and cohort differences prevent attributing differences to exhaustion removal.

## Integrity and diagnostics

The complete 713,768,673-byte source (9,804,452 rows) matched the publisher MD5 `c10d4dca6efaa4dba2e39ae29add61ae` after four connection retries. GRCh37 coordinate plus exact non-palindromic allele pair matched 991,705 source variants to the previously prepared 991,707-marker set; all four LDSC analyses retained 991,705 allele-valid markers. No eligible duplicate IDs or ambiguous reference keys were reported. Matched per-row sample size was 405,979 throughout. Source matched rows were streamed, not retained as a separate file. See [receipt](stream-20261029T161051-7498.json), [protocol](PROTOCOL.md), [stream](stream.py), [run](run.sh), and [complete LDSC log](ldsc-20261029T161051-7498.log). An independent [prepared-input audit](common_inputs_receipt.json) confirms all four input files contain the same 991,707 unique SNPs with compatible non-palindromic allele pairs. The streamed source marker list was not saved, so the exact LDSC output marker sets were not independently exported; equal final counts alone cannot establish their identity.

Physical-function observed-scale SNP heritability was 0.3233 (SE 0.0171), with intercept 1.1924 (0.0523). The MVP and HFRS intercepts were 1.1358 (0.0109) and 1.1895 (0.0149); the physical-function cross-trait intercepts were −0.0017 (0.0111) with MVP, −0.1056 (0.0111) with FI, −0.2927 (0.0128) with FFS, and 0.0053 (0.0105) with HFRS. These diagnostics and possible UK Biobank overlap limit interpretation; an intercept does not by itself identify the source of bias. Inputs for clinical insomnia include approximate sample-size and p-derived Z assumptions documented in the main report. The earlier published insomnia MR analysis for the same composite and its overlap limitations are described in the [source review](SOURCE.md); MR and genetic correlation answer different questions.

## Relationship to the ML prototype

The new physical-function comparison is a post-hoc sensitivity outcome, not independent training or validation data for the existing predictor. Adding it as a feature or evaluation target after inspecting these results would change the model and risks leakage through UK Biobank overlap and dependent derived correlations. A future extension needs prespecified feature construction, source and cohort holdouts, paired uncertainty, and a genuinely independent validation cohort; do not train on AlphaGenome outputs under the accepted output terms.
