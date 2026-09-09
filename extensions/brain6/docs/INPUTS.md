# Input contracts and native settings

## Canonical GWAS

Tab-separated, optionally gzip:

`SNP CHR BP A1 A2 BETA SE P N EAF INFO`

Positions are 1-based. A1 is the effect allele. Chromosomes are autosomes 1–22. BETA/SE must share the same effect scale. Frequencies refer to A1. Missing frequency/INFO is permitted only under explicit source policy; missingness must not be replaced with an invented frequency.

The source-card `column_map` maps each canonical field to its actual source column. `path`/`sha256` bind exact local bytes. `effect_scale` is `beta` or `odds_ratio`; the latter also needs `se_scale=log_odds`. `n_semantics` is explicit. The source-card code supports declared N from a mapped column, a constant, or case/control counts; inspect the template and `brain6/gwas.py` before registering a source that needs a custom adapter.

Real source cards require `access_permitted: true`, a real citation and source URI. This is permission/provenance information, not an application for controlled access. Do not populate those fields for data you are not authorized to use.

## Result matrix and decisions

Default matrix columns: `sleep_trait disease_trait rg se p fdr analysis_status`.

The current repository may use other names or labels. Copy `configs/brain6.json` and map the **observed** fields. The default expects exactly 396 rows and exactly 72 unique rows for the six chosen disorders across the twelve sleep phenotypes. It does not infer a valid QC status from a missing column.

Decision table: `disease_trait sleep_trait role reason`.

Roles: `PRIMARY`, `SECONDARY`, `NO_ELIGIBLE_PAIR`. The freeze command verifies eligibility against the original full-family FDR and demands a primary or explicit none entry for every disorder. `NO_ELIGIBLE_PAIR` is disallowed when eligible rows exist; a decision to exclude an eligible disorder for a different reason needs a separately documented protocol revision.

## Native job example

For SuSiE, make a copy of `configs/susie.template.json` and fill actual locus/LD paths, trait number, phenotype scale, sample size and rationale, appropriate case fraction or sourced quantitative SD, and reviewed numerical policies.

An example bindings file:

```json
{
  "inputs": {
    "pair_lock": "/absolute/path/brain6-pairs.lock.json",
    "source_card": "/absolute/path/trait.source.json"
  },
  "outputs": {
    "fit": "fit.rds",
    "pip": "pip.tsv",
    "credible_sets": "credible_sets.tsv",
    "diagnostics": "diagnostics.json",
    "semantic_status": "status.json"
  }
}
```

`make-job` automatically binds standard settings-file inputs and method wrappers. Additional custom source files should go into `inputs`. Use complete absolute paths in settings; jobs execute in isolated temporary directories, not your shell's working directory.

For R settings, add `expected_packages` after validating the production environment:

```json
{"expected_packages": {"coloc": "YOUR_TESTED_VERSION", "susieR": "YOUR_TESTED_VERSION"}}
```

The strings above are deliberately not real versions. Obtain actual versions using `Rscript scripts/native_environment.R`, run the native smoke/validation suite on that host, and archive a genuine `renv.lock` or equivalent environment snapshot. Exact version checks are available in `common.R` but this delivery does not fabricate a native lock for software it could not install.

### LDSC

Use the actual Python interpreter for the pinned original LDSC environment, not automatically the extension's Python. Register `ldsc_script` or `munge_script`, summary statistics, and the correct LD-score/weights prefix. Genomic SEM uses directory forms in its R wrapper. Binary liability-scale h² needs sourced population and sample prevalence. Observed-scale rg is not magically improved by inventing a prevalence.

### PLINK and signed LD

Native adapters target PLINK **1.9 syntax**, including `.clumped` output and `--r square`. Do not substitute a `plink2` binary and assume all formats/flags match.

A signed LD manifest requires:

```json
{
  "kind": "signed_r",
  "counted_allele": "BIM_A1",
  "ancestry": "EUR",
  "genome_build": "GRCh37",
  "synthetic": false,
  "matrix_path": "/absolute/path/signed.ld",
  "matrix_sha256": "VERIFIED_SHA256",
  "bim_path": "/absolute/path/locus.bim",
  "bim_sha256": "VERIFIED_SHA256"
}
```

The BIM defines exact matrix order and counted alleles. Blocks use `LOC CHR START STOP`; they are unique and nonoverlapping with 1-based inclusive bounds. The runner preserves the declared block and reports all missing/unalignable variants.

### PLACO import and collation

Legacy import JSON must include `reviewed`, `pair_id=insomnia__adhd`, terminal `status`, actual `method`, source `path`/`sha256`, `input_fingerprint`, `run_fingerprint`, `expected_rows`, and an explicit `column_map` into:

`SNP CHR BP Z1 Z2 P_PLACO status`

Statuses for collation are `TESTED`, `NUMERICAL_FAILURE`, `EXCLUDED_EXTREME_Z`. Do not map a legacy failed value to TESTED or drop failures. Compare the source's actual extreme-variant policy to the new one before pooling. Imported evidence remains legacy, not independent replication.

The completed chunk manifest uses `path sha256`. `collate-placo` requires the exact precomputed total row count and selected-family size. The output table is `results.tsv.gz` plus `statistics.sqlite` and `status.json`. A scientifically rejected collation can be audited but must not pass into `deep-followup`.

### LAVA

Bind `input_info`, `sample_overlap_file` (or a reviewed no-overlap justification), PLINK reference prefix and `loci_file`. All filenames in the input-info table are checksum-bound. Regions must match the declared genome build. The supplied adapter deliberately handles one locus at a time; this trades repeated startup/reference loading for bounded execution and checkpoints. Benchmark batching on the actual host before increasing parallelism.

### Genomic SEM

`mode=covariance`: parallel vectors of paths, trait names, sample prevalence and population prevalence, genuine reference/weights directories and jackknife block count. Continuous-trait prevalence entries are JSON `null`, to be represented as R missing values. Review binary sample-size/prevalence conventions carefully.

`mode=model`: real covariance RDS, a reviewed model string and model-review flag. The fit is saved and marked review-required; it is not auto-promoted based only on exit status.

### MR

Supply canonical exposure and outcome tables for independently LD-clumped, exposure-selected instruments. Set both `instruments_ld_clumped` and `instrument_selection_exposure_only` only after actually performing these steps. Declare overlap review. Neither Boolean is a substitute for an instrument-selection provenance ledger. Run the reverse direction separately and apply family-level testing across the declared directions.

### MATLAB pleioFDR

The launcher needs an official source checkout, a reviewed official text configuration, two preprocessed MAT files, and a matching reference MAT file. It copies **code only** into isolated staging and forces output/cache paths inside staging. The input `alignment_manifest` requires reviewed ancestry/build, exact reference row-order SHA, reference variant count and file records for all three MAT files. Both trait records must explicitly match the reference row-order SHA.

This manifest documents prior alignment; it does not create it. Preparing those native MAT arrays from raw GWAS and proving their reference-index alignment remains a separate data-preparation/validation task. Declare exact expected output filenames from the chosen upstream version. The launcher leaves scientific acceptance review-required even if MATLAB succeeds.

## Functional evidence

Locus table: `pair_id locus_id status trait_h4 pip_max`.

Molecular table: `pair_id locus_id gene_id evidence_type molecular_h4 tissue cell_type source`.

Allowed evidence types: `coding`, `eqtl_coloc`, `sqtl_coloc`, `chromatin_link`, `expression`, `nearest_gene`. H4/PIP values must be finite probabilities or explicit `NA`. The molecular table is not a list of guessed nearest genes. Each source must be traceable to a reviewed public/authorized resource and a version.

Competitive annotation table: `pair_id cell_type locus_id score selected stratum`.

Each pair×cell×locus is unique. `selected` is 0/1. `score` is a predeclared annotation score. Every stratum contains selected and matched control loci. Build strata and locus units before evaluating enrichment; not afterward to make a cell significant. The test is a one-sided matched competitive permutation and corrects across all supplied pair×cell tests.
