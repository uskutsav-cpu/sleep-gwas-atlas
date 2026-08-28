# sleep-gwas-atlas

An auditable pipeline for a locked 45-trait Sleep/Circadian Genetic Atlas.

- **Scope:** 12 sleep/circadian traits and 33 non-sleep traits.
- **Phase 0:** source verification, ancestry/build checks, harmonization, QC ledgers, and HapMap3 munging.
- **Phase 1:** LDSC SNP heritability, predefined QC, 396 sleep-by-disease genetic correlations, FDR, and a provenance-stamped heatmap.
- **Later phases:** MiXeR, LAVA, pleiotropic-locus discovery, Genomic SEM, fine-mapping, colocalization, molecular/cell-type/pathway integration, and cautious causal inference.

## Scientific status

The analysis scope is now frozen in `config/analysis_panel.tsv` at exactly 45 traits. `config/analysis_panel.lock.json` locks the ordered trait identities and domain counts, while every generated readiness or Phase 1 inclusion table records the full manifest SHA-256. A trait therefore cannot be silently swapped while preserving a 45-row count.

The source registry now provides evidence-backed public-source and verified
schema records for **45/45** selected traits. Six explicit substitutions resolve
the prior access or ancestry blockers without changing the locked trait IDs:
FinnGen R9 supplies Finnish/European MS, asthma, T2D, wide-IHD-as-CAD-proxy,
and melanoma endpoints, while Burren 2024 `GCST90435144` supplies the NFE
telomere-length stratum. Their phenotype and power differences are carried in
the manifest rather than hidden. Source verification does not by itself claim
that local materialization, harmonization, LDSC, or h² QC passed. On a fresh
clone with ignored raw data absent, the expected readiness summary is:

```text
source_verified: 45 / 45
harmonization_ready: 0 / 45
ldsc_ready: 0 / 45
liability_h2_ready: 0 / 45
phase1_pass: 0 / 45
```

The historic 86-row `config/traits.tsv` and `config/panel_45_selection.tsv` remain as provenance for the candidate-selection process. They are not production inputs. No real LDSC results are versioned on this branch; synthetic smoke-test artifacts are not scientific results.

The former access blockers and the evidence-backed public substitutions are
documented in `docs/remaining_source_access_blockers.md` and
`docs/public_eur_substitution_plan.md`. The pipeline did not submit a request,
sign a data-use agreement, or affirm a download checkbox on the user's behalf.

All 45 registered sources have now been exercised locally through
harmonization, HapMap3 munging, and h² estimation. All 12 sleep traits and 31
non-sleep traits pass the predefined LDSC h² gate; T2D fails the intercept gate
and melanoma fails the h²-Z gate. The locked 396-pair sleep×non-sleep family is
complete as 372 primary Phase 1 pairs plus 24 explicitly labelled
`QC_FAILED_SENSITIVITY` pairs. There are 153 primary associations at FDR<0.05
when correction is applied over all 396 tests (155 under the primary-only
372-test correction). Sensitivity rows are retained for completeness but are
excluded from primary inference.

The complete multivariable LDSC covariance gate is also finished locally: all
45 traits have S/Rg/I matrices, all 1,035 lower-triangle estimates, and a
1,035-by-1,035 sampling-covariance matrix. The raw estimated S matrix is not
positive semidefinite and the required 1,082-block jackknife is numerically
ill-conditioned, so those warnings are carried into model selection. HDL also
fails the intercept threshold in GenomicSEM even though the standalone Python
two-step h2 estimator narrowly passed it; the original Phase-1 result remains
unchanged and HDL is excluded from confirmatory SEM. See
`docs/full_covariance_results.md` for the audited diagnostics.

An odd-chromosome discovery/even-chromosome confirmation analysis was then run
for the 42 traits passing the multivariable LDSC input gate. None of ten
discovery-derived candidate models passed held-out fit and residual
admissibility together, so no factor model was promoted and factor GWAS/Q_SNP
remain intentionally unrun. See `docs/genomic_sem_validation_results.md`.

The production LAVA workflow is now pinned and contract-complete, including
all 112,275 local-univariate tests, overlap correction, all 396 eligible pair
definitions, per-locus checkpoints, and a strict result validator. The real
analysis has not started because the official recommended UK Biobank LD v1.1
reference requires 15 GiB uncompressed while the current volume has only about
4–6 GiB free. The older 1,000 Genomes reference is not substituted because LAVA
0.1.5 warns of local-h2 bias and type-I error inflation. See
`docs/lava_workflow.md`.

The insomnia liability conversion uses a rounded 30% frequent-complaint
prevalence from the primary study. The longevity liability conversion uses the
study's phenotype-defined 10% survival-tail prevalence, while AF uses the ESC
adult-population 3% approximation and ADHD preserves the original Demontis
GWAS 5% convention. Schizophrenia uses the primary study's 1% lifetime-risk
convention, and its release-specific half-effective-N field is explicitly
doubled for LDSC. Stroke uses a rounded 3% adult ever-stroke prevalence
supported by CDC surveillance. Crohn disease and ulcerative colitis use rounded
0.3% and 0.5% high-prevalence European approximations, respectively.
Rheumatoid arthritis preserves the published 1%
liability-scale convention used by Ha et al. 2021 for European and East Asian
RA genetic analyses. These ignored local artifacts are not counted
as present in the fresh-clone summary, are not
the canonical full-panel result tables, and do not satisfy any finish-line
acceptance gate.

## Locked panel contract

All production entry points read `config/analysis_panel.tsv`. The manifest records, per trait:

- phenotype definition;
- selected source GWAS and dataset version;
- registered source ID and raw filename;
- phenotype type and sample counts;
- ancestry and genome build;
- PMID/DOI when verified;
- population prevalence and citation when required; and
- source-verification declaration.

The separate `config/gwas_schemas.tsv` records the literal variant, allele,
effect, SE, P, frequency, INFO, and sample-size fields for each selected
source-and-trait pair. This distinction is required when one archive contains
multiple phenotypes with different effect columns.
A source can be identity/checksum verified while its schema remains pending;
that state cannot advance to `HARMONIZATION_READY`. A verified schema that
lacks usable rsIDs is recorded as `SCHEMA_VERIFIED_REQUIRES_VARIANT_MAPPING`:
the columns are understood, but harmonization requires the audited strategy in
`config/variant_mapping_plans.tsv` and the provenance-checked map reproduced by
`scripts/00_setup.sh`. The current mapper is intentionally restricted to the
1,184,461-variant intersection of the pinned EUR LDSC coordinates and HapMap3
allele list; it is not a genome-wide resolver.

Verified GRCh38 sources likewise require an exact row in
`config/liftover_plans.tsv`. Setup pins UCSC's hg38-to-hg19 chain by byte count,
MD5, and SHA-256. The harmonizer rejects unmapped or ambiguous points,
reverse-complements alleles on reverse-strand mappings, and records the chain
hash and losses in the QC ledger; build conversion is never inferred from an
rsID alone.

Unknown facts are written as `UNRESOLVED`, never guessed. Validate the contract with:

```bash
python3 scripts/00_validate_panel.py
```

Changing a trait identity or order requires an explicit, reviewed update to the lock file. Metadata can be curated without changing scope, but downstream artifacts retain the manifest hash used to produce them.

## Readiness model

`scripts/10_phase0_audit.py` derives one evidence-backed ledger instead of overloading a `TODO/CURATED` flag:

```text
SOURCE_VERIFIED
  -> HARMONIZATION_READY
  -> LDSC_READY
  -> LIABILITY_H2_READY
  -> PHASE1_PASS
```

The stages mean:

- `SOURCE_VERIFIED`: the selected public source, URLs, checksum, file mapping, and registry declaration agree.
- `HARMONIZATION_READY`: source verification plus resolved phenotype/publication, EUR ancestry, validated GRCh37/hg19 input or a checksum-valid registered liftover, required sample metadata, and materialized raw input.
- `LDSC_READY`: harmonization QC ledger and HapMap3 munged summary statistics exist.
- `LIABILITY_H2_READY`: LDSC-ready, with a cited matching population prevalence for binary phenotypes. For continuous traits this is the final-h² readiness gate.
- `PHASE1_PASS`: final-h² readiness plus the predefined LDSC h² QC verdict `PASS`.

Generate the ledger:

```bash
python3 scripts/10_phase0_audit.py \
  --out results/tables/trait_readiness.tsv \
  --strict
```

`--strict` fails when a manifest row claims `SOURCE_VERIFIED` without matching registry evidence. Later-stage blockers remain explicit rows rather than being silently dropped.

## Pinned environments

Exact top-level versions are recorded in `environment/tool_versions.tsv`.

The workflow/CI runtime uses Python 3.11.11, Snakemake 8.30.0, NumPy 1.26.4, pandas 2.2.3, and Matplotlib 3.9.4. PLINK 1.90b7.7 and R 4.3.3 are pinned for stages that introduce them. Create that environment with:

```bash
conda env create -f environment/workflow.yml
conda activate sleep-gwas-atlas-workflow
```

LDSC runs separately under Python 3.9.23 at the pinned CBIIT commit `6c673952cee74bd5c57aef1555a03b1c015399a0`. The setup downloads several GB of reference data, so inspect the source plan before running it:

```bash
bash scripts/00_setup.sh
export PYTHON_BIN=python3
export LDSC_PYTHON=.ldsc-env/bin/python
export LDSC_DIR=ldsc
```

Complete multivariable LDSC and Genomic SEM use a separate R 4.3.3
environment. `environment/genomicsem.yml` pins the conda-packaged R
dependencies; the setup script additionally pins `simsalapar` 1.0-13 and the
official GenomicSEM commit recorded in `environment/tool_versions.tsv`:

```bash
bash scripts/25_setup_genomicsem.sh
```

LAVA shares this R runtime. Its separate setup pins the exact LAVA commit and
all additional CRAN dependencies by source checksum; it extracts the official
2,495-locus definition but never downloads the large LD reference implicitly:

```bash
bash scripts/30_setup_lava.sh
bash scripts/32_download_lava_reference.sh  # storage preflight only
```

## Snakemake workflow

The default workflow validates the locked panel and produces a source/readiness ledger without downloading raw GWAS data:

```bash
snakemake --cores 1
```

Real-data work is opt-in. Add only already reviewed traits to `phase0_traits` in `config/workflow.yaml`, materialize their registered raw files, and invoke:

```bash
snakemake --cores 1 phase0
snakemake --cores 1 h2
snakemake --cores 1 phase1_rg
snakemake --cores 1 full_covariance
snakemake --cores 1 lava_inputs
```

`full_covariance` always reads all 45 ordered manifest traits, not the optional
`phase0_traits` subset. It runs GenomicSEM's multivariable LDSC with the pinned
EUR LD scores, applies documented liability conversions to binary traits, and
produces the 45×45 genetic covariance/correlation/intercept matrices plus the
1,035×1,035 sampling-covariance matrix required for overlap-aware Genomic SEM.
For 45 traits GenomicSEM requires 1,082 jackknife blocks; this is a long-running
production target and the resulting high-block-count diagnostic must remain an
explicit robustness warning.

The `lava` target is intentionally separate from downloading reference data.
It requires the complete official UK Biobank LD v1.1 chromosome files, resumes
through per-locus checkpoints, and publishes results only after the complete
locked test family can be collated and validated:

```bash
snakemake --cores 1 lava
```

After the run, validate the complete export and derive the downstream SEM input
ledger:

```bash
python3 scripts/26_validate_covariance.py
python3 scripts/27_genomicsem_trait_qc.py
```

The existing shell entry points remain available and are called by Snakemake:

```bash
bash scripts/02_munge.sh sleepdur breast_cancer
bash scripts/03_h2_qc.sh sleepdur breast_cancer
bash scripts/04_rg.sh --h2 results/tables/h2_summary.tsv
```

Binary traits require a cited, phenotype-matched population prevalence for final liability-scale h². `--observed-scale` remains an interim QC option, not a substitute for the final table.

## MiXeR terminology

`scripts/05_collate.py` reports `ldsc_n_eff_h2` and `ldsc_n_eff_h2_gt_12000` only as the CDG3 LDSC screening arithmetic. It never labels that calculation `mixer_pass`. Actual MiXeR eligibility must come from a successful univariate MiXeR model and its diagnostics before bivariate MiXeR is attempted.

## Reproducible smoke test and CI

The smoke test uses deliberately fake data and exercises:

- the exact 45-row panel contract;
- tabular and whitespace-delimited raw inputs;
- BETA and OR effect formats;
- harmonization/QC ledgers;
- all 45 synthetic h² logs;
- readiness-backed Phase 1 selection;
- the full synthetic 12 × 33 = 396 correlation family;
- report, metadata, and heatmap generation; and
- prominent synthetic-output watermarks.

```bash
PYTHON_BIN=python3 bash scripts/run_smoke_test.sh
```

GitHub Actions runs panel/source validation, Python and shell syntax checks, a Snakemake dry run, and the full synthetic smoke test on every push and pull request.

## Canonical Phase 0/1 outputs

```text
config/analysis_panel.tsv                    locked 45-trait analysis manifest
config/analysis_panel.lock.json              ordered trait/domain identity lock
config/public_gwas_sources.tsv               source URLs, versions, access, hashes
config/variant_mapping_plans.tsv             audited mapping strategy for eight sources
config/liftover_plans.tsv                    pinned hg38-to-hg19 plans for two sources
ref/hm3_grch37_variant_map.tsv.gz             reproducible ignored GRCh37 HapMap3 identity map
ref/hg38ToHg19.over.chain.gz                  checksum-pinned ignored UCSC chain
results/tables/analysis_panel_provenance.tsv manifest and trait-set fingerprints
results/tables/trait_readiness.tsv            per-trait readiness stages/blockers
data/harmonized/<trait>.qc.txt                filter-by-filter QC ledger
data/munged/<trait>.sumstats.gz               canonical HapMap3 LDSC input
results/tables/h2_summary.tsv                 h², SE, Z, intercept, ratio, QC
results/tables/phase1_inclusion.tsv           locked-panel Phase 1 eligibility
results/tables/rg_matrix.tsv                  rg, SE, P, and Benjamini-Hochberg FDR
results/figures/fig2_rg_heatmap.png           provenance-stamped Phase 1 heatmap
```

Generated raw data, intermediate files, logs, tables, figures, caches, and smoke-test artifacts are ignored by Git. A release must freeze selected outputs separately with checksums, configs, provenance, tool versions, and the producing commit hash.

The complete scientific finish line is tracked separately from infrastructure
progress in `docs/atlas_v1_finish_line.md`. Audit it at any time with:

```bash
python3 scripts/99_atlas_acceptance.py --report-only
```

Without `--report-only`, the command fails until every atlas-v1.0 scientific and
release gate is supported by real, non-synthetic artifacts.
