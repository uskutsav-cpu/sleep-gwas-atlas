# Scientific protocol and inference boundaries

## Study hierarchy

Keep the original 45-trait core, 396-test family, prior extension, Track B A/B/control and their original hashes immutable. Brain6 is a separately documented, **post-screen follow-up**. The six disorder identities are fixed by the user's request; measured sleep partners are selected from actual prior results.

A = snoring × parental lifespan, B = insomnia × ADHD, and CONTROL = insomnia × frailty are not renamed or overwritten. The default automation refuses new ADHD Pair-B PLACO computation. Import only verified original bytes, actual method, input fingerprint, run fingerprint, full row count and failure mapping. Method differences between original PLACO and PLACO+ must be disclosed before comparing or pooling outputs.

The other five disorder tracks are not assumed significant merely because they appeared in an old proposal. The measured 12×6 table includes every positive, negative, nonsignificant, invalid and sensitivity-only row. Preserve original family FDR. A valid `NO_ELIGIBLE_PAIR` decision is allowed instead of forcing a disease through a positive-discovery pipeline.

## Pair and hypothesis families

Record one primary and at most one secondary partner per disorder, with an explicit reviewer and reason. The selection lock records that the global results were already observed. It is not an original preregistration and it does not automatically confer selective-inference validity.

The default locus screen uses the conventional 5×10^-8 genome-wide reference threshold divided by the number of selected follow-up pairs. This is an additional between-pair Bonferroni factor, **not** a proof of exact error control for a selected family, and not a replacement for variant-level or local-region testing rules.

PLACO results expose within-pair BH and a conservative between-pair factor separately. Numerical failures occupy denominator slots for conservative accounting but are still output as `NA`, not assigned genuine null results. Whole-pair acceptance also depends on the explicit failure-rate gate. Extreme-Z exclusions are an explicit method-defined exclusion set, not silent missingness. Their trait-level signals must not vanish from the source audit.

Cell enrichment corrects across the supplied pair×cell tests. Choose the complete annotation/testing universe before interpretation; do not pass only cells that looked enriched. Conditional/local tests and MR directions require their own reviewed full-family corrections. The MR adapter does not silently manufacture a study-wide FDR from one pair's P values.

## Input biology and measurement

- Keep source-study identity, phenotype definition, clinical vs proxy status, ancestry, build, case/control numbers, variant-level N and access/provenance explicit.
- SE for an odds ratio must be on the log-odds scale before using log(OR). Never divide log(OR) by an OR-scale SE.
- Effective sample size, raw case-control sample size and total sample size are not interchangeable. Fine-mapping explicitly requires a reviewed locus-level N and allowed variation, rather than an unreported median replacement.
- Existing GRCh37/GRCh38 coordinates must match. No build is relabelled as another build. Actual liftOver would require source/target references, multimap/lost-variant accounting and allele revalidation; it is not automatically implemented here.
- Exclude unsupported ambiguous alleles and duplicates transparently. Both records of a duplicate are quarantined rather than selecting a favorable row. Dense locus work does not replace full variants with a HapMap3-only table.
- Genome-wide overlap/coverage and allele alignment must precede correlation/pleiotropy analysis. Matching only rsID is insufficient when position/build disagree.

## Native methods

PLACO+ parameter estimation uses genome-wide eligible Z scores and marginal P values. Do not estimate variance or trait correlation independently in small chunks, or restrict fitting to significant SNPs. Chunking divides expensive per-variant evaluation, not the statistical model. Source code version and SHA must be recorded.

Use signed genotype correlations, matching SNP order and counted alleles for SuSiE/coloc. r² is not signed r. Both trait effects are aligned to the same LD allele. Reference missingness is quantified; no lead SNP is quietly lost. Intact LD blocks are not arbitrarily thinned to fit memory. Memory/LD coverage failures remain explicit.

SuSiE must converge under the reviewed settings. Keep every credible set, its members, posterior inclusion probabilities, purity and diagnostics. Large GWAS–reference mismatch is a failed QC state, not a reason to flip alleles until H4 increases. Out-of-sample LD does not justify estimating residual variance as if it were in-sample LD.

Coloc runs at the signal level and retains H0, H1, H2, H3 and H4, plus prior sensitivity. A maximum H4 across many signals/priors is only a descriptive prioritization quantity, not independent repeated confirmation. Lack of a credible set in one trait is insufficient evidence, not proof of distinct causal variants. High statistical colocalization support is not proof of a biological mechanism.

LAVA's low local-univariate power gate is distinct from an execution error. Existing failed legacy LAVA output remains diagnostic-only. New local analyses use the prespecified region list, ancestry/reference checks and explicit sample-overlap handling. No automatic switch to a different method to rescue significance.

Genomic SEM requires genuine genetic and sampling covariance matrices. Off-diagonal sampling covariance cannot be invented from a pairwise rg table. Model fit is a review stage: the wrapper does not automatically declare a latent factor real or causal.

MR instruments are selected using the exposure GWAS and LD reference, independently of the outcome P values and PLACO results. Run both directions as separately declared analyses. Review strength, sample overlap, harmonization, heterogeneity, outliers, phenotype definition and estimator assumptions. An insignificant MR-Egger intercept is not proof that pleiotropy is absent; agreement of estimators is not causal proof.

## Molecular/cell evidence

Only sourced coding, eQTL/sQTL colocalization, chromatin links, expression and proximity evidence enter the normalized evidence tables. Proximity or expression alone does not identify a causal gene. The default H4=0.8 labels are descriptive report categories, not a universal evidence guarantee.

Competitive cell enrichment uses independently defined loci and matched control strata, not thousands of correlated SNPs treated as independent observations. Public expression/accessible-chromatin scores, covariates and controls must be established outside the enrichment function. This is not a replacement for a calibrated S-LDSC/MAGMA analysis or donor-level single-cell replication.

Use a broad, common cell-type comparison set where the resource supports it, then additional disease-relevant regional resources. Do not prefill neuronal/microglial/dopaminergic "positive" cells for particular disorders. Shared loci or genes across disorders need actual evidence and a consistent coordinate/identifier definition.

## Terminal states

- `PASS`: this method's declared computation/QC criteria passed; no automatic causal claim.
- `NO_SIGNAL`: a valid tested analysis did not yield an eligible signal at its threshold.
- `INSUFFICIENT_EVIDENCE`: underpowered/absent eligible signals or a scientific-review boundary; not a negative result.
- `FAILED_QC_NOT_CONSUMED`: output preserved, but disallowed as accepted downstream evidence.
- `FAILED_OR_BLOCKED`: data, software, numerical or resource obstacle, reported with its reason.
- `SYNTHETIC`: always separated from empirical data and forbidden in an unmarked production run.

A pipeline may complete an honest accounting of no signals or unresolved loci. Completing the project does not require producing six positive mechanisms.
