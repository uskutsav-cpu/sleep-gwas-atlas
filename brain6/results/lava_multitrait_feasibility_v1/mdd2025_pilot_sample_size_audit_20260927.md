# MDD2025 LAVA trait-only pilot: sample-size semantics audit

Status: **SCIENTIFIC_INTERPRETATION_HOLD**. This is an audit of the isolated
111-locus diagnostic pilot, not a change to its frozen configuration or to any
canonical LAVA analysis. The pilot's numerical outcome remains unverified while
`/Volumes/Extreme SSD` is unmounted.

The frozen pilot materializer writes source `NEFF` into the LAVA sumstats `N`
column, while the input-info file supplies 412,305 cases and 1,588,397 controls.
The locked Howard 2019 MDD source card declares `n_semantics: total`; its
canonical comparison therefore uses a different sample-size definition.
The source describes `NEFF` as an effective sample size and also supplies
per-variant `NCAS` and `NCON`. Those fields have different meanings and must not
be interchanged without a method justification. The frozen pilot config SHA-256
is `d57fb8ecc0a9bc1198a5c41df5dbd74d7d750b32ed7c27501bb151a730b2cbbd`;
the decision-rule SHA-256 is
`85e8587e24ce148d7e8f919824f55e2b36b4aa4777389658fbe49c59d39b8089`.
The [PGC MDD2025 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11829167/)
describes effective sample size in terms of an equivalently powered balanced
case/control study. That makes it especially important to justify combining
source `NEFF` with the observed, unbalanced case fraction in LAVA.

There is also an independent family-level stop: canonical v3 has 3,720
`NOT_RUN` cells, of which 589 are MDD. Even if a replacement MDD GWAS made
all 2,495 MDD loci testable, the other traits would still contribute 3,131
`NOT_RUN` cells, exceeding the frozen 873-cell ceiling. This follows directly
from `lower_bounds.json`. A full MDD-only screen cannot reopen the seven-trait
family gate and should not be launched solely because this pilot looks better.

LAVA 0.1.5 [documents](https://github.com/josefin-werme/LAVA/blob/main/README.md)
the sumstats `N` as the number of samples and states that the input-info `N`
is not used for analysis. Its
[`process.locus` implementation](https://github.com/josefin-werme/LAVA/blob/main/R/input_processing.R)
passes the sumstats `N`, LD frequency, and input-info case proportion to
`process.binary`, then uses the mean locus `N` in sampling variance and local
heritability calculations. The
[`process.binary` implementation](https://github.com/josefin-werme/LAVA/blob/main/R/binary_processing.R)
sets the reconstructed number of cases to `N * case.prop` and genotype counts
to functions of `N` and reference allele frequency. Thus, the pilot is modeled
as a pseudocohort with effective `N` but the overall case fraction. That is not
the literal per-variant case/control sample represented by `NCAS`/`NCON`.

This does **not** establish that substituting `NCAS + NCON` is correct for a
heterogeneous meta-analysis. Source `NEFF` may incorporate cohort weighting or
other design effects not recoverable by summing the source counts. Nor does
the existing source card establish that `NEFF` plus the aggregate case fraction
is a valid LAVA 0.1.5 representation. The observed local-h² testability could
depend on this choice. The pilot therefore cannot, on its own, authorize a full
2,495-locus MDD2025 sensitivity run until this interpretation is resolved.

When the SSD returns, first audit the existing four worker shards and all
source/reference receipts without restarting work. Keep the pilot results and
decision rule intact. Record numerical outcomes as **diagnostic only**. Before
any new LAVA run, establish a separate, prospectively versioned sample-size
protocol from the source's methods and LAVA's binary model. It must define the
per-variant `N`, the corresponding case proportion, and why they jointly
represent this meta-analysis; compare the choice with both `NEFF` and
`NCAS + NCON` using the same predetermined loci. Do not choose the model by
which one produces more testable loci. Canonical v3, v2, and roundoff artifacts
remain frozen.
