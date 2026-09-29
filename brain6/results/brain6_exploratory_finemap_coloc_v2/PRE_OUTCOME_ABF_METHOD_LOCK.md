# Brain6 exploratory single-signal ABF method lock (v2)

**Frozen before inspection of any candidate association posterior.** This is
separate from the v1 SuSiE-RSS attempt. The v1 100k-UKB BCOR LD matrices failed
their fixed numerical gate in all 32 scheduled inputs; the raw matrices had
out-of-range correlations up to 1.01956 and negative eigenvalues, and none
passed the locked 0.005-shrinkage Cholesky check. No SuSiE posterior or
SuSiE-coloc result was generated. v2 does not relabel that failed attempt.

## Eligible analyses and data

The 25 protected pair-specific PLACO candidate IDs and GRCh37 bounds remain
unchanged. For each candidate and each of its two traits, use the preserved
normalized dense SQLite GWAS source and source-card allele/effect convention.
The source files for insomnia, ADHD, MDD, bipolar disorder, schizophrenia,
and Parkinson disease contain binary log-odds beta and standard error; the
Dashti long-sleep released BOLT-LMM effect is not admitted to this binary ABF
branch because its model/effect-scale translation remains unresolved.

Only exact rsID and position agreement with the pinned European LAVA UKB v1.1
GRCh37 reference are accepted. Alleles are oriented to reference A1 through
exact, swap, complement, or complement-swap matching. A/T and C/G alleles,
multiallelic or invalid alleles, missing or nonfinite beta/SE, and SE<=0 are
excluded. Reference MAF must be >=0.01 and NOBS >=90,000. At least 500
oriented SNPs, at least 50% of the eligible reference SNP universe, and at
least 80% of valid source GWAS rows must remain. Both traits must satisfy
these gates, with at least 500 shared SNPs, for coloc. Input N is audited but
is **not** used by ABF because the source supplies beta and SE for each SNP;
v1's within-locus effective-N homogeneity filter is specific to the failed
SuSiE-RSS single-N approximation and is not applied here. Density failures
remain explicit `NOT_RUN_INPUT_QC` rows. The LD reference supplies only a
fixed marker universe and allele standard, not a covariance matrix for ABF.

## Single-signal fine-mapping

Run `coloc::finemap.abf` version 5.2.3 on each admitted trait-locus with
`type="cc"`, source log-odds beta, `varbeta=SE^2`, and per-SNP prior p1=1e-4.
The package's default binary effect prior SD=0.2 is retained. Report the
null-model posterior and each variant's unconditional `SNP.PP`; conditional
variant PIP is `SNP.PP / (1 - PP.null)`. Report a 95% credible set only when
the locus association posterior `1-PP.null` is >=0.80; otherwise status is
`NO_ASSOCIATION_SUPPORTED_CREDIBLE_SET`. If emitted, the credible set is
the smallest descending-conditional-PIP prefix totaling >=0.95. Report its
size and included variants. These are **single-causal-variant-model ABF**
posteriors, not SuSiE/multiple-signal GWAS fine-mapping. They are exploratory
for all unpromoted PLACO candidates and never constitute LAVA confirmation.

## Trait–trait colocalization

Run `coloc::coloc.abf` on the same shared oriented SNP intersection for each
pair passing input QC. Use `type="cc"`, beta, varbeta, and priors
p1=p2=1e-4, p12=1e-5. Report PP.H0, PP.H1, PP.H2, PP.H3, and PP.H4. Also
run fixed p12 sensitivity at 1e-6 and 5e-5; never select a favorable prior
after seeing results. Descriptive `ABF_MODEL_SUPPORT` requires PP.H4>=0.80,
PP.H4/(PP.H3+PP.H4)>=0.80, and both single-trait association posteriors
>=0.80 on the shared SNP universe. Even with that label, unknown multiple
signals and sample overlap preclude a claim of confirmed shared causality.
There is no familywise significance or independent-replication claim.

## Scientific boundaries

ABF assumes at most one causal variant for each trait in a tested interval.
No genotype-valid dense LD is available to prove that assumption or perform
reliable multi-signal conditioning here. High H4 can be misleading with
multiple signals. A null or low posterior is model-conditional, not evidence
that the PLACO candidate is false. No result is allowed to alter the frozen
LAVA family, promote a candidate, or become `CONFIRMATORY_FINE_MAPPING` or
`CONFIRMATORY_TRAIT_COLOC`. All 50 candidate-trait and all 25 pair-candidate
units receive result or explicit NOT_RUN rows. Numerical results remain
`EXPLORATORY_SINGLE_SIGNAL_ABF`.

Method sources: [original ABF coloc paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383),
[coloc fine-map documentation](https://chr1swallace.github.io/coloc/reference/finemap.abf.html),
[multiple-signal limitation and SuSiE extension](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009440).
