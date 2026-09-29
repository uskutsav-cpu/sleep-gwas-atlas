# Brain6 exploratory GWAS fine-mapping and trait–trait coloc: pre-outcome lock

This branch was fixed after inspecting only source metadata, variant counts,
alleles, coordinates, and sample-size distributions. No local association
statistics, fine-mapping posterior, or coloc posterior had been inspected. The
protected 25 pair-specific PLACO candidates and their original GRCh37 bounds
are taken verbatim from `candidate_evidence_25.tsv`; the 20 geographic groups
do not replace those pair-specific analysis units. Canonical LAVA status and
candidate promotion remain unchanged.

## Inputs and admission

* GWAS: the seven checksum-recorded dense normalized SQLite snapshots in
  `brain6-work/preparation-v1/normalize_<trait>/variants.sqlite`. Source
  allele and effect conventions are those in `config/gwas_schemas.tsv`.
* LD: LAVA UK Biobank v1.1 European GRCh37 reference (100,000 samples),
  accessed with pinned LAVA 0.1.5 `load.reference` and `read.ld`. This reference
  was selected for ancestry, build, dense marker coverage, and a pre-existing
  project lock, before any fine-mapping outcomes.
* Runtime: `susieR` 0.14.2 and `coloc` 5.2.3 in the isolated SSD R library,
  R 4.6.1 for posterior calculations; LAVA 0.1.5 under its existing R 4.3.3
  environment for LD extraction. The exact package/library manifest is saved
  with outputs.
* SNP keys must have identical rsID and position across GWAS and LD, valid
  A/C/G/T alleles, and unambiguous exact, swap, strand-complement, or
  complement-swap orientation. Effects and LD are oriented to LD A1. A/T and
  C/G SNPs are excluded even when frequency is available. Duplicate keys,
  contradictory coordinates, and nonfinite beta/SE/N are errors.
* LD reference SNPs require MAF >=0.01 and NOBS >=90,000. A trait-locus must
  retain at least 500 oriented LD SNPs, at least 50% of eligible reference
  SNPs in its interval, and at least 80% of valid GWAS rows. The 20% tail of
  normalized effective N around the within-locus median is excluded; at
  least 80% of otherwise eligible SNPs must survive this N-consistency gate.
  These gates are source and input QC, not association significance gates.
* Long sleep is excluded from formal GWAS PIP and trait–trait coloc here:
  its frozen binary-tail BOLT-LMM source lacks verified variant-level analyzed
  N and an exact mapping from the released model to a case-control RSS
  likelihood. Its valid dense rows are still audited. This exclusion cannot
  be reversed after seeing posteriors within v1.
* MDD and Parkinson trait-loci fail the prospective density requirement in
  the input audit (MDD reference coverage 14–20%; Parkinson 6–15%). Their
  missing posteriors are explicit input-QC failures, not null findings.

## Primary exploratory fine-mapping

For each admitted trait-locus, run `susieR::susie_rss` on source-derived signed
Z = beta/SE and allele-oriented LD. The effective N input is the median of
the existing trait-specific normalized N field after the N-consistency gate;
it is a binary-GWAS approximation, **not** claimed to be a per-SNP literal
participant count. Use L=10, 95% credible-set coverage, `min_abs_corr=0.5`,
`max_iter=1000`, `tol=1e-3`, and fixed residual variance. A constant 0.005
diagonal shrinkage is applied to reference LD before fitting, and the raw and
regularized LD matrices must be symmetric, finite, unit-diagonal, and positive
definite after shrinkage at Cholesky precision. The unregularized matrix is
checked for symmetry, finite entries, and unit diagonal; a singular raw
matrix is recorded rather than silently altered. No tuning by candidate
result is permitted.
Record convergence, every credible set, each variant PIP, credible-set purity,
and the fraction of the locus reference universe retained. A converged fit
without a credible set is reported as `NO_CREDIBLE_SET`, never as a fine-mapped
causal variant. All rows remain `EXPLORATORY_FINE_MAPPING` until the existing
frozen promotion workflow formally confirms a locus.

## Trait–trait coloc

Only the six insomnia–ADHD candidate loci have both traits meeting the dense
input gate at lock time. Both source files report log-odds effect and SE for
their named effect allele. Refit each trait with the same SuSiE parameters
on the identical SNP intersection, then run `coloc::coloc.susie` on those
paired fits, with priors p1=p2=1e-4,
p12=1e-5. Report **every** signal-pair comparison and H0–H4 posterior vector;
do not select a favorable pair as the region's sole result. Define exploratory
regional support only when both traits have converged 95% credible sets of
purity >=0.5 and at least one signal pair has PP.H4 >=0.80 and
PP.H4/(PP.H3+PP.H4) >=0.80. This is descriptive, not a corrected frequentist
discovery claim. Also calculate `coloc.abf` under its explicit one-causal-
variant-per-trait assumption as a **sensitivity analysis** over the same
variants and default binary effect SD=0.2, with priors p1=p2=1e-4 and p12=1e-5.
Report p12 sensitivity values 1e-6 and 5e-5 without switching the primary
prior. A single-signal ABF posterior alone never upgrades an unsupported
multi-signal coloc result. Failed or absent SuSiE credible sets yield explicit
`NOT_ESTIMABLE` primary coloc status.

## Scientific and provenance limits

The source GWAS may have cohort overlap and the LD is a reference rather than
in-sample matrix. Binary trait RSS uses an effective-N approximation. These
analyses cannot make the failed frozen LAVA endpoint pass, cannot promote a
PLACO candidate, and are not independent replication. Credible sets are
conditional on the admitted marker set, not on all variants in the region.
The scripts record versions, source hashes, exact SNP sets, QC failures,
parameters, and SSD output paths. Any changed threshold, prior, source, LD,
or method requires a new branch version before inspecting new results.

Method sources: [SuSiE-RSS documentation](https://stephenslab.github.io/susieR/reference/susie_rss.html),
[coloc-SuSiE methods](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009440),
[coloc dataset and prior documentation](https://chr1swallace.github.io/coloc/reference/check_dataset.html).
