# sleep-gwas-atlas

An auditable pipeline for a locked 45-trait Sleep/Circadian Genetic Atlas.

- **Scope:** 12 sleep/circadian traits and 33 non-sleep traits.
- **Phase 0:** source verification, ancestry/build checks, harmonization, QC ledgers, and HapMap3 munging.
- **Phase 1:** LDSC SNP heritability, predefined QC, 396 sleep-by-disease genetic correlations, FDR, and a provenance-stamped heatmap.
- **Later phases:** MiXeR, LAVA, pleiotropic-locus discovery, Genomic SEM, fine-mapping, colocalization, molecular/cell-type/pathway integration, and cautious causal inference.

## Scientific status

The analysis scope is now frozen in `config/analysis_panel.tsv` at exactly 45 traits. `config/analysis_panel.lock.json` locks the ordered trait identities and domain counts, while every generated readiness or Phase 1 inclusion table records the full manifest SHA-256. A trait therefore cannot be silently swapped while preserving a 45-row count.

The source registry currently provides evidence-backed public-source records for **42/45** selected traits. All 45 selections now have registry and schema rows, but the final three deliberately remain pending: the selected IMSGC multiple-sclerosis discovery meta-analysis has no located full-statistics release, the exact DIAGRAM type-2-diabetes route requires affirmative terms acceptance, and the selected dbGaP melanoma archive requires Authorized Access. This is deliberately broader than the old binary `CURATED` count: source verification does not claim that ancestry/build/schema checks, local materialization, harmonization, LDSC, prevalence evidence, or h² QC have passed. On a fresh clone with ignored raw data absent, the expected readiness summary is:

```text
source_verified: 42 / 45
harmonization_ready: 0 / 45
ldsc_ready: 0 / 45
liability_h2_ready: 0 / 45
phase1_pass: 0 / 45
```

The historic 86-row `config/traits.tsv` and `config/panel_45_selection.tsv` remain as provenance for the candidate-selection process. They are not production inputs. No real LDSC results are versioned on this branch; synthetic smoke-test artifacts are not scientific results.

The exact evidence and user/institution actions needed for the last three source
records are documented in `docs/remaining_source_access_blockers.md`. The
pipeline never submits a request, signs a data-use agreement, or affirms a
download checkbox on the user's behalf.

Thirteen registered sources have also been exercised locally end to end. The
ignored Jansen 2019 insomnia, Campos 2020 snoring, Yengo 2018 BMI, Deelen 2019
longevity, Timmers 2019 parental-lifespan, Evangelou 2018 SBP, Nielsen 2018
atrial-fibrillation, Demontis 2023 ADHD, Trubetskoy 2022 schizophrenia, Mishra
2022 stroke, de Lange 2017 Crohn-disease and ulcerative-colitis strata, and
Ishigaki 2022 rheumatoid-arthritis inputs produced 6,077,635/10,862,567,
7,168,629/11,010,158, 1,973,592/2,336,269, 1,175,095/8,856,352,
6,663,125/9,085,648, 5,964,514/7,088,067, 10,246,131/12,149,979,
5,692,669/6,774,224, 6,341,702/7,659,767, 1,176,288/7,511,476,
1,144,234/9,570,787, 1,144,294/9,588,016, and 9,659,407/13,297,690 retained
harmonized variants, respectively. All thirteen were HapMap3-munged and passed
the predefined LDSC h² gate in an explicitly named checkpoint; the gated Phase
1 checkpoint contains the 2 locally ready sleep traits × 11 locally ready
non-sleep traits = 22 pairs.
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
