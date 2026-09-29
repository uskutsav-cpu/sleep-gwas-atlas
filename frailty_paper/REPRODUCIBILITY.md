# Frailty package reproducibility

Run commands from the repository root. The available package-level commands are:

```sh
make -C frailty_paper help
make -C frailty_paper audit-frozen
make -C frailty_paper manifest
make -C frailty_paper verify-manifest
make -C frailty_paper audit-manifest-metadata
make -C frailty_paper audit-hlma-portal-crosswalk
make -C frailty_paper validate-plan
make -C frailty_paper plot-frozen-fi
# Explicit preview only when the runtime differs from the pinned environment:
make -C frailty_paper plot-frozen-fi-preview
make -C frailty_paper checklists
make -C frailty_paper import-review-exports
make -C frailty_paper validate-review
make -C frailty_paper audit-bibliography
make -C frailty_paper audit-manuscript-claims
make -C frailty_paper test
FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper preflight
make -C frailty_paper smoke
make -C frailty_paper physical-component-qc
make -C frailty_paper fi-correction-sensitivity
make -C frailty_paper latent-correction-sensitivity
make -C frailty_paper supplementary-tables
```

`audit-frozen` verifies provenance hashes for the locked sleep-atlas outputs,
reproduces the read-only 12-row FI result extract idempotently, and audits the
separate discovery-extension outputs. It does not rerun either GWAS analysis.
The extension audit may report missing core-checkpoint artifacts; that is a
documented property of this checkout, not a passing full reproduction.

`manifest` rebuilds `manifests/all_acquired_resources.tsv` from already
verified local source files and retained PubMed XML. It performs no network
requests and does not download any data.

`audit-hlma-portal-crosswalk` consumes the saved public HLMA listing API
snapshot and rebuilds the 23-object link metadata and 46-row snATAC-fragment
archive-path crosswalks. It does not request any archive bytes. OMIX IDs and
portal URLs are link evidence only; donor-label conflicts and release-metadata
discrepancies remain explicitly flagged, and per-file CNP mapping, object
identity, genome build and reuse terms are not inferred.

`verify-manifest` checks every manifest path, byte count, and SHA-256 without
changing the inventory or source files. It rejects duplicate resource IDs and
paths that escape the repository root. When data are stored on a separately
mounted volume through repo symlinks, pass that exact approved root with
`FRAILTY_EXTERNAL_STORAGE_ROOT=/mounted/path make verify-manifest`; all other
external targets remain rejected.

`audit-manuscript-claims` checks selected headline counts and effect estimates in
the provisional manuscript against their source tables, verifies the exact text
reported, and writes a source-hashed JSON record under `analysis/`. It does not
validate scientific eligibility, screening decisions, causal interpretation, or
complete manuscript accuracy.

`preflight` combines resource-manifest verification, required-metadata completeness
checking, frozen analysis-plan validation, cohort-overlap validation, review-decision
validation, bibliography and selected quantitative-claim integrity checking, other
registered source audits, and the package test suite. `audit-manifest-metadata` separately checks all 22 required
columns, nonblank metadata, unique resource IDs, positive file sizes, and SHA-256
syntax. Explicit unknown values are allowed; empty values are not. Set `FRAILTY_EXTERNAL_STORAGE_ROOT` when registered files
resolve through repository symlinks to a mounted data volume. It does not
download data, change screening decisions, run GWAS analyses, or regenerate
results; tests may create and remove their own temporary fixtures.

## External derived-output workspace

When the internal volume has less than the locked 20 GiB workspace minimum,
route the large harmonized and munged files to a mounted volume with at least
20 GiB free. The setup helper refuses to move non-empty output directories;
it preserves the empty tracked placeholders beside the repository and links
the output paths to the external workspace.

```sh
export FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'
export FRAILTY_ANALYSIS_WORKSPACE="$FRAILTY_EXTERNAL_STORAGE_ROOT/analysis-workspace"
make -C frailty_paper storage-configure \
  FRAILTY_EXTERNAL_STORAGE_ROOT="$FRAILTY_EXTERNAL_STORAGE_ROOT" \
  FRAILTY_ANALYSIS_WORKSPACE="$FRAILTY_ANALYSIS_WORKSPACE"
make -C frailty_paper storage-check FRAILTY_ANALYSIS_WORKSPACE="$FRAILTY_ANALYSIS_WORKSPACE"
mkdir -p "$FRAILTY_ANALYSIS_WORKSPACE/tmp"
export TMPDIR="$FRAILTY_ANALYSIS_WORKSPACE/tmp"
export REF_DIR="$FRAILTY_ANALYSIS_WORKSPACE/reference"
export CONDA_PKGS_DIRS="$FRAILTY_ANALYSIS_WORKSPACE/conda-pkgs"
export LDSC_ENV_DIR='.ldsc-env' LDSC_DIR='ldsc'
bash scripts/00_setup.sh
export LDSC_PYTHON="$LDSC_ENV_DIR/bin/python"
```

The storage check must pass immediately before any analysis stage that writes
large outputs. Keep `TMPDIR` and `REF_DIR` on the same mounted volume for
temporary and LD-reference data; package downloads are cached there through
`CONDA_PKGS_DIRS`. The small pinned LDSC source checkout and Python environment
remain in the repository's ignored local directories. Outputs outside
`data/harmonized` and `data/munged` still need their own capacity check. The
helper verifies existing raw-data links and preserves the local empty output
directories as `.local-preserved` before linking.

For a separate `frailty_v1` primary-FI run, keep generated data and logs in its
versioned workspace and keep the small sealed result table in the repository:

```sh
export FRAILTY_V1_WORKSPACE="$FRAILTY_ANALYSIS_WORKSPACE/frailty_v1"
export HARMONIZED_DIR="$FRAILTY_V1_WORKSPACE/harmonized"
export MUNGED_DIR="$FRAILTY_V1_WORKSPACE/munged"
export LDSC_LOGDIR="$FRAILTY_V1_WORKSPACE/logs"
export H2_OUT="$PWD/frailty_paper/results/frailty_v1/h2_summary.tsv"
export TMPDIR="$FRAILTY_ANALYSIS_WORKSPACE/tmp"
export LDSC_PYTHON="$PWD/.ldsc-env/bin/python"
export LDSC_DIR="$PWD/ldsc"
export REF_DIR="$FRAILTY_ANALYSIS_WORKSPACE/reference"
bash scripts/02_munge.sh frailty
bash scripts/03_h2_qc.sh frailty
```

Each stage checks the actual harmonized, munged and temporary destinations for
the 20 GiB minimum before writing. Keep each `frailty_v1` run in its own
versioned directories; do not point it at the legacy frozen-atlas paths.

`frailty_paper/.venv/bin/python scripts/22_download_latent_frailty_catalog.py`
resumes the eight accession-indexed latent-frailty Catalog downloads. It keeps
unresolved phenotypes labeled by accession, retains the Catalog YAML and
`md5sum.txt` sidecars, and promotes each file only after the official MD5
matches. Use `--accession GCST...` to limit a run or `--jobs 1` for serial
transfer.

`scan-gwas-sources` streams each acquired GWAS summary-statistics gzip to EOF,
records its header and row count, and counts missing/impossible basic numeric
fields. It writes `manifests/gwas_source_scan.tsv`, never filters or rewrites
source data, and does not replace harmonization, duplicate-variant, allele,
build or phenotype-specific QC.

`validate-plan` parses the frozen YAML and checks the plan checksum, inherited
input checksums, primary FI identity, 12-trait panel shape, and all-396
multiplicity contract. `test` runs the package unit/integration checks,
including tamper-detection tests for the lock validator. `smoke` combines
plan validation with the read-only frozen-result audits. These commands use
the package virtual environment (`.venv`); install its documented local
dependencies before running them in a fresh checkout.

`plot-frozen-fi` creates a descriptive forest plot from the existing 12-row
frozen sleep–FI extract and checks the Python/Matplotlib versions against
`environment/tool_versions.tsv`. It does not rerun LDSC or independently
reprocess the missing harmonized/munged inputs. If the versions do not match,
the command stops. `plot-frozen-fi-preview` explicitly permits an unpinned
preview and records the actual and required versions and output checksums in
`analysis/preliminary_frozen_fi_global_rg_forest.provenance.json`. The current
preview is not publication-ready; rerender it under the pinned environment
before any publication use.

`physical-component-qc` rebuilds the five-row structural QC summary from the
acquisition registry, full-stream/basic-value scan, and duplicate-key audit.
It verifies size, hash, row-count, required-field and uniqueness agreement,
then marks each component's analysis eligibility `BLOCKED` while build,
ancestry/QC, effect model, and exact generating-study provenance remain
unresolved. It performs no filtering or harmonization.

`fi-correction-sensitivity` verifies the frozen FI extract against its audit
manifest and reports Bonferroni-adjusted P values using the same locked 396-pair
family as a secondary comparison with the inherited BH q-values. It does not
recompute LDSC, change the primary correction, or independently reprocess the
missing original harmonized/munged files.

`latent-correction-sensitivity` verifies the fixed 12 × 7 latent-factor RG
family and reports Bonferroni-adjusted P values within its locked denominator
of 84, alongside the existing secondary-family BH q-values. It records both
correction outcomes and preserves unknown exact participant overlap. BH remains
the prespecified primary correction for this sensitivity family; the target
does not rerun LDSC or change any thresholds.

`supplementary-tables` reruns the latent correction summary and rebuilds the
currently available source-linked supplementary tables, status index and hash
provenance. Existing non-identical outputs are protected; review the source
changes before explicitly choosing the builder's `--replace-outputs` option.

`checklists` deterministically builds full item-level tables for PRISMA 2020
(including subitems), PRISMA 2020 Abstracts, PRISMA-S, and STREGA. Each row
includes a status, planned manuscript section, evidence/gap, next action, and
manuscript page/line field. Page/line values remain `NOT_DRAFTED` until a
manuscript exists.

After authorized Embase, Scopus, Web of Science, or PsycINFO RIS/CSV exports
are placed in `review/manual_exports/` with the database named in each
filename, `import-review-exports` records each original file's SHA-256 and
normalizes bibliographic fields without modifying the source exports. The
review-record builder then merges normalized records with PubMed and
deduplicates on shared PMID, DOI, or normalized-title keys, retaining every
source-record identifier in the audit and screening queue. If reviewer data
already exist, the builder refuses to overwrite them.

`validate-review` is a read-only check of dual screening: decision vocabulary,
distinct reviewers, completion of both votes, disagreement adjudication,
title/abstract include routing, full-text locations, and reasoned exclusions.
It reports the current counts without changing screening decisions.

Frailty harmonization, QC, heritability, replication, LAVA, pleiotropy,
fine-mapping, colocalization, functional follow-up, final figures, and
manuscript targets are intentionally not exposed as runnable commands yet.
Their required frailty inputs, source permissions, overlap and build decisions,
and upstream result gates are incomplete; current status and blockers are
recorded in `STATUS.md`, `BLOCKERS.md`, and `config/analysis_plan_decisions.tsv`.
Large acquisition remains paused until at least 20 GiB are available as
required by the project plan. `scripts/09_require_storage.sh` enforces that
gate before either core or frailty-specific GWAS downloads and at the top-level
collection entry point. A blocked storage check exits before download commands
run. Existing core and extension analyses remain immutable.

The commands above depend on the repository's Python environment. Setup is
documented in `scripts/00_setup.sh`; that script installs packages and may use
the network, so these audit commands do not invoke it automatically.

The workflow environment is pinned to Python 3.11.11 in
`environment/tool_versions.tsv`. `make -C frailty_paper setup-env`
(`FRAILTY_PYTHON=...` and optionally `FRAILTY_VENV_DIR=...`) checks that version
before creating or updating the selected virtual environment, then installs the pinned
`requirements-pipeline.txt`, `requirements-workflow.txt`, and
`requirements-acquisition.txt` direct dependencies.
Set `FRAILTY_PYTHON` to the executable for Python 3.11.11 when it is not the
default `python3`. Set `FRAILTY_VENV_DIR` to an empty target when building a
separate clean environment. The setup script refuses to layer a different Python
version over an existing target. Check the selected target's `bin/python --version` before treating a local
test run as pinned. The data-validation tests for figure inputs do not render
figures. The current Figure 4 can be refreshed from a lock-matched full-family
audit with `scripts/52_refresh_figure4_current_audit.py`; it uses the Python
standard library plus `rsvg-convert` for SVG, PDF and PNG outputs. Figure 5
and the other plotting scripts continue to use Matplotlib from the pinned
requirements.

## Exact snapshot for the tested Python workflow environment

The full Python package suite was also run in a scratch Conda environment on macOS Apple silicon. A platform-specific recreation snapshot is recorded in `environment/frailty_workflow_py311_osx_arm64_conda_explicit.lock` (88 Conda package builds with SHA-256 values) and `environment/frailty_workflow_py311_osx_arm64_pypi.lock` (45 PyPI package versions). `analysis/pinned_workflow_environment_validation_2026-09-23.json` records both lock-file hashes, the hashes of code/config/requirement inputs, and the current 70/70 test log hash. S16 checks the snapshot hashes before building.

To recreate that Python workflow environment on macOS Apple silicon from the repository root:

```sh
conda create --yes --prefix work/conda-envs/frailty-py311-replay --file environment/frailty_workflow_py311_osx_arm64_conda_explicit.lock
work/conda-envs/frailty-py311-replay/bin/python -m pip install -r environment/frailty_workflow_py311_osx_arm64_pypi.lock
cd frailty_paper
../work/conda-envs/frailty-py311-replay/bin/python -m unittest discover -s tests -v
```

This snapshot is scoped to `osx-arm64`; it is not a portable lock for other operating systems. Conda package artifacts are SHA-256 pinned. PyPI package versions are exact-pinned but their downloaded artifact hashes are not recorded. The snapshot excludes R and other separately installed system binaries and does not rebuild the canonical `frailty_paper/.venv`; the standard setup script remains the direct-requirements bootstrap described above.

A clean Python 3.11.11 venv was additionally created with:

```sh
make -C frailty_paper setup-env \
  FRAILTY_PYTHON=../work/conda-envs/frailty-py311/bin/python \
  FRAILTY_VENV_DIR=work/conda-envs/frailty-paper-py311-clean-2026-09-23
```
Its `pyvenv.cfg` confirms `include-system-site-packages = false`. The original setup snapshot passed `pip check`, all 70 package tests, and the complete `preflight` when selected with `PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python`; current 94-test suite and pip-check results are recorded separately below.
The resolved 69-package version freeze, test/setup/preflight logs and hashes are
recorded in `environment/frailty_workflow_py311_clean_venv_freeze_2026-09-23.txt`
and `frailty_paper/analysis/clean_venv_setup_validation_2026-09-23.json`.
PyPI artifact hashes remain unavailable in the freeze; the active canonical
`frailty_paper/.venv` remains Python 3.13.13 and was not replaced. The isolated venv remains available at `work/conda-envs/frailty-paper-py311-clean-2026-09-23`. A 94/94 test suite and `pip check` passed on 2026-09-24 as a historical snapshot; the latest integrated Python 3.11.11 preflight passed 97/97 (see the captured log below). The original setup validation remains a historical 70-test snapshot. Earlier evidence is in `analysis/clean_venv_py311_unittest_2026-09-24_94.log`, `analysis/clean_venv_py311_pip_check_2026-09-24.log`, and `analysis/clean_venv_setup_validation_2026-09-23.json`.


Phase 16 sensitivity status and the combined Supplementary Table 15 are rebuilt with:

```sh
make -C frailty_paper phase16-sensitivity-summary
make -C frailty_paper supplementary-tables
```

The sensitivity matrix is a source-linked status synthesis across seven current headline conclusions and 17 planned domains. It reports completed, partial, unresolved, blocked and not-justified states; it does not create estimates or treat unavailable tests as null findings.

The complete 2026-09-23 preflight logs are historical 91/91-test snapshots in the isolated Conda Python 3.11.11 environment (`work/conda-envs/frailty-py311`) and project Python 3.13.13 environment. The 350-row, 94/94 integrated pinned preflight is a historical snapshot at `analysis/preflight_reverification_2026-09-24_350_94.log`. The 351-row integrated preflight including the source-abstract and 28-claim manuscript audits passed 97/97 tests under isolated Python 3.11.11; the 06:21–06:26 capture is retained at `analysis/preflight_post_abstract_audit_py311_2026-09-24.log`. A fresh preflight from source HEAD `6b906ba` at 07:01–07:03 UTC also passed all 351 resource checks and 97/97 tests; its log is `analysis/preflight_current_head_py311_2026-09-24.log`. Both runs used a mixed worktree, with screening decisions still at zero and participant intersections unknown. The standalone 94/94 suite logs are historical snapshots: `analysis/clean_venv_py311_unittest_2026-09-24_94.log` and `analysis/project_python313_unittest_2026-09-24_94.log`. Pinned log: `frailty_paper/analysis/preflight_phase1_manifest_py311_2026-09-23.log` (SHA-256 `af89de6109b1cc699778511a75b17b39f8dd0531745f3edb4ff3143e8c23cb50`); project log: `frailty_paper/analysis/preflight_phase1_manifest_2026-09-23.log` (SHA-256 `5a681ff2ac8d56e85eb6c52314421fb7aa6cad53ca2457b553888ed360926ab3`). The 3.11 environment passes `pip check`; the canonical `frailty_paper/.venv` remains Python 3.13.13. Earlier 85-, 81-, 77-, 74- and 70-test logs are historical snapshots. `make -C frailty_paper audit-bibliography` regenerates a hash-linked JSON check for Pandoc citation keys, duplicate BibTeX entries, missing cited entries and unused references.


The latest complete pinned preflight is `frailty_paper/analysis/preflight_current_head_py311_2026-09-24.log` (97/97 tests; 351/351 resources; 935 missing-abstract source checks; 28 manuscript claims passed; zero review decisions; exact participant intersections unknown). It ran from source HEAD `6b906ba` in a mixed worktree and is engineering/integrity evidence, not a scientific replay. After the S16 builder and regression expectation update, the standalone pinned Python 3.11.11 suite also passed 97/97; log: `frailty_paper/analysis/current_head_py311_unittest_2026-09-24.log`. The earlier full preflight is `frailty_paper/analysis/preflight_phase1_manifest_py311_2026-09-23.log` (SHA-256 `af89de6109b1cc699778511a75b17b39f8dd0531745f3edb4ff3143e8c23cb50`); the project Python 3.13.13 preflight is `frailty_paper/analysis/preflight_phase1_manifest_2026-09-23.log` (SHA-256 `5a681ff2ac8d56e85eb6c52314421fb7aa6cad53ca2457b553888ed360926ab3`). Both verified 347/347 external-root manifest files, the 22-column metadata audit, locked plan, overlap ledger, review queue, bibliography, Q_SNP and reporting locators, and all 91 tests. Reproduce the pinned run with `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' PYTHON='../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python' make -C frailty_paper preflight`. This is engineering/integrity evidence; it does not resolve the documented scientific blockers or constitute a full scientific replay.

The preflight also runs `audit-latent-qsnp`, which reproduces the source-method check in `frailty_paper/analysis/latent_factor_qsnp_pruning_audit.tsv` and its provenance JSON. It uses the Foote et al. Bonferroni Q_SNP threshold (5e-8/7) and ±1 Mb flanks. The seven acquired factor files contain no variants below that threshold, so this exact exclusion rule removes no rows; no input or existing estimate is rewritten. `make -C frailty_paper audit-latent-qsnp` runs the audit independently.


The regenerated S7 table carries a pair-level indicator-overlap annotation: insomnia × general factor is flagged as a direct construct indicator/part-whole comparison; sleepiness rows note the related TIR indicator with unresolved loading. The 84-row table is source-linked in `paper/supplementary_tables/supplementary_tables.provenance.json`.


`make -C frailty_paper audit-reporting-locators` checks every checklist row for a current line reference or explicit absent-section status, verifies referenced files and line bounds, and writes `frailty_paper/analysis/reporting_checklist_locator_audit_2026-09-23.json`. The current audit passes for 92 rows and 93 line references. This mechanical check does not assess whether each cited passage semantically satisfies the reporting item. Checklist regeneration preserves the existing `manuscript_page_line` value by `(framework, item_id)` so routine rebuilding does not erase reviewed locators.


The Phase 1 metadata audit at `analysis/phase1_manifest_metadata_audit_2026-09-23.md` is a historical 347-row snapshot; the current 350-row metadata and full file-integrity checks pass in `analysis/acquisition_manifest_reverification_2026-09-24_350.log`. The latest complete 2026-09-23 preflight logs both record 91/91 tests and pass the integrated checks; the current 2026-09-24 pinned preflight also passes with the 350-row manifest and 94 tests (`analysis/preflight_reverification_2026-09-24_350_94.log`). Historical standalone 94-test suite records are in `analysis/project_python313_unittest_2026-09-24_94.log` and `analysis/clean_venv_py311_unittest_2026-09-24_94.log`. The latest pinned integrated preflight is the separately captured 97-test run noted above. The older 85-test preflight predates the metadata audit.

## Current test entry point — 2026-09-26

The Brain6 power-audit regression tests in `tests/test_audit_brain6_sleep_power_pilot.py` use standard-library `unittest` so package discovery does not require pytest. `make -C frailty_paper test` passes all 157 discovered tests in the package venv; captured output is `frailty_paper/analysis/make_test_unittest_2026-09-26_0656.log`. The original pytest-based run is retained as historical validation evidence only.
