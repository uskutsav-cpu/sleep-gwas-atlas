# Brain6 v0.2: six-disorder deep-analysis extension

An additive, executable extension for `uskutsav-cpu/sleep-gwas-atlas`.

**Scope:** ADHD, major depressive disorder, schizophrenia, bipolar disorder, Alzheimer's disease and Parkinson's disease. Select sleep partners from the **actual completed atlas**, preserve the original screening FDR, and reuse compatible existing evidence without modifying Track B A/B/CONTROL.

**Delivery status:** Python components and synthetic software integration have been executed and tested. Native R, PLINK, LDSC, MATLAB and Snakemake production analyses have **not** been executed in the delivery environment. This package is not a finished empirical study, a verified live-repository merge, or a publication-ready scientific release. See `docs/IMPLEMENTATION_STATUS.md` and the delivery validation directory.

## New in v0.2

This release builds on the recovered v0.1 package, rather than replacing the live atlas. **206 Python tests pass and four real-native smoke tests skip because their runtimes are missing in this delivery environment.** These are extension tests, not a combined test of your inaccessible current repository.

New: exposure-only MR candidate selection and verified clump/outcome alignment; native MR family accounting; cohort-aware sourced replication; full-cis QTL indexing; resumable molecular QTL → SuSiE/coloc orchestration; PIP-weighted BED cell scores; and a complete stage-coverage audit. Null six-disorder panels are supported without inventing partners. Native source changes invalidate cache fingerprints; receipt inventories and strict JSON validation are hardened. The outer installer supports only checksum-matched v0.1 upgrades; edited local files remain protected.

Read **[the v0.2 runbook](docs/V02_RUNBOOK.md)** for the new commands and exact evidence boundaries, and **[upgrade notes](docs/UPGRADE_NOTES.md)** before reusing existing artifacts. Native code is supplied, but no R/PLINK/LDSC/LAVA/PLACO/SuSiE/coloc/GenomicSEM/MATLAB empirical computation was executed here.

## Start here

This directory is self-contained. It does not replace the repository's root README, scripts, manifests, analysis panel, Snakefile, or scientific outputs.

```bash
cd extensions/brain6
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
python -m brain6 --help
```

The above installs Python dependencies, not R, MATLAB, PLINK or the original genetics programs. The Python tested-version snapshot is in the delivery validation folder. The package's dependency ranges are compatibility declarations, **not** a fully reproduced native environment lock.

### Run the synthetic demonstration

```bash
python -m brain6 demo --out work/SYNTHETIC-demo
python -m brain6 heatmap --ranking work/SYNTHETIC-demo/SYNTHETIC_ranking.tsv \
  --out work/SYNTHETIC-demo/SYNTHETIC_heatmap.png
```

The demonstration normalizes two small artificial GWAS-like tables, quarantines duplicates, harmonizes both allele orientations, joins and indexes them, makes deterministic PLACO input chunks, and extracts a 72-row brain subset from an artificial 396-row table. It **does not** run or imitate PLACO, LDSC, SuSiE, coloc, LAVA or an empirical MR study. Diagnostic IVW arithmetic is explicitly labelled synthetic.

Use a new output directory for each demonstration. Existing artifacts are immutable, not silently overwritten.

## Apply to the existing project

Use the delivery's `apply_brain6.py` from the existing repository root. Without flags it only checks the repository identity, payload checksums and path conflicts. `--apply` adds the extension on a dedicated feature branch. `--commit` and `--push` are separate explicit actions.

The installer never uses `git add .`, force-push, reset, stash, or automatic replacement. Existing staged changes prevent automatic committing. Unstaged and untracked work outside the added paths is preserved. Commits use **your existing Git identity**.

Do not copy an entire old project archive over your active checkout. These are additive source files, not a replacement for current data or code.

## Real-data entry point: inspect, rank, then select

From this directory:

```bash
mkdir -p work
python -m brain6 audit --repo ../.. --out work/local-checkout-audit.json
cp configs/brain6.json work/brain6.reviewed.json
```

Inspect the real matrix's header and QC labels. Set `matrix_columns`, `matrix_qc_column` and `primary_qc_values` in your copied configuration to match it. The code does not guess that an unknown status means primary QC passed.

```bash
python -m brain6 rank \
  --matrix ../../results/tables/rg_matrix.tsv \
  --config work/brain6.reviewed.json \
  --out work/brain6-measured-ranking.tsv
cp configs/pair_decisions.template.tsv work/pair-decisions.reviewed.tsv
```

Choose one primary, optionally one secondary, measured and QC-eligible sleep partner per disorder. Keep an explicit `NO_ELIGIBLE_PAIR` entry where no pair passes the selection gate; do not force an Alzheimer's or Parkinson's positive result. Then:

```bash
python -m brain6 freeze \
  --matrix ../../results/tables/rg_matrix.tsv \
  --config work/brain6.reviewed.json \
  --decisions work/pair-decisions.reviewed.tsv \
  --reviewer 'Your actual reviewer name' \
  --out work/brain6-pairs.lock.json
```

**Important correction to earlier planning:** chronotype–schizophrenia, chronotype–bipolar, insomnia–Alzheimer's and duration–Parkinson's were examples in a proposal, not verified highest-ranked result rows. Their partners remain `UNRESOLVED` in this delivery. Insomnia–ADHD is an existing legacy pair, not new independent evidence.

The freeze records that selection occurs **after** the global screen. A new family threshold of `5e-8 / number_of_selected_pairs` does not solve all selective-inference questions. It does not replace the original 396-test FDR or retrospectively alter the A/B threshold.

## Data preparation: reuse the same engine across all selected pairs

Create reviewed source cards by copying `configs/source.template.json`; absolute paths are recommended. Record the original file SHA-256, source/build/ancestry, effect scale and sample-size meaning. Put trait IDs and source-card paths in a JSON registry:

```json
{
  "insomnia": "/absolute/path/insomnia.source.json",
  "adhd": "/absolute/path/adhd.source.json",
  "mdd": "/absolute/path/mdd.source.json"
}
```

Include **every selected sleep and disorder trait**, not only the three illustrative entries above.

```bash
python -m brain6 assemble-preparation \
  --pair-lock work/brain6-pairs.lock.json \
  --sources work/sources.reviewed.json \
  --output-root work/preparation-v1 \
  --runtime-out work/preparation.runtime.json --reviewed
python -m brain6 pipeline --runtime work/preparation.runtime.json
```

Source normalization is shared across pairs. Disk-backed joins, 20,000-row batches and on-disk sorting avoid a single giant pandas merge. Free-space checks remain mandatory; lower RAM use does not remove the need for adequate disk.

The generated DAG deliberately skips new PLACO chunk creation for `insomnia__adhd`. Import its verified canonical evidence using `import-legacy`, after checking actual method and source fingerprints. Do not describe a PLACO result as PLACO+ merely to match a newer branch.

## Native statistics

A reviewed **settings JSON** describes a method. A **job JSON** binds every file and declares outputs. Jobs run synchronously in isolated transaction directories and stop on missing engines, invalid data, output schema failures or excessive numerical failures.

```bash
python -m brain6 make-job --method susie --name example_locus_trait1 \
  --settings work/susie.reviewed.json --bindings work/susie.bindings.json \
  --out-root work/job-specs --reviewed
python -m brain6 run-job \
  --job work/job-specs/example_locus_trait1.job.json --root work/native-v1
```

For the command above, `bindings.json` contains `inputs` (additional file paths) and `outputs` (relative output filenames including `semantic_status: status.json`). Standard method input paths, PLINK BED/BIM/FAM members and LD-score chromosome files are automatically discovered and hashed.

See `docs/INPUTS.md` for an exact example and `configs/` for each adapter's settings. Templates contain deliberate unresolved values. Passing `--reviewed` is an assertion by the operator, not an automated substitute for scientific review.

### Genome-wide PLACO+ with chunk resumption

After preparation, use the *actual* paired artifact and chunks directories, the new pair lock, and a reviewed local copy of the official PLACO R source with its independently recorded SHA-256:

```bash
python -m brain6 plan-placo \
  --pair /absolute/path/pair_insomnia__mdd \
  --chunks /absolute/path/chunks_insomnia__mdd \
  --pair-lock work/brain6-pairs.lock.json --pair-id insomnia__mdd \
  --source /absolute/path/PLACO_v0.2.0.R --source-sha256 YOUR_VERIFIED_SHA256 \
  --out work/placo-mdd-v1 --reviewed
python -m brain6 execute-placo --plan work/placo-mdd-v1/plan.json
```

`insomnia__mdd` here is command syntax, **not an automatic selection**. The pair must actually occur in the reviewed lock. The planner estimates parameters using genome-wide inputs once, calls the official per-variant method in bounded chunks, and collates using the complete eligible denominator. Numerical failures keep `NA` P/Q values. BH is not recalculated independently for each chunk.

### Clumping through fine-mapping and colocalization

Copy `configs/followup.template.json` and bind the accepted PLACO collation artifact, pair index, official PLINK reference, reviewed LD blocks, trait-specific sample sizes, R-package versions and numerical policies.

```bash
python -m brain6 deep-followup --settings work/followup.reviewed.json
```

This sequential runner performs native PLINK clumping, groups leads into intact blocks, generates signed LD, aligns **both** traits to the counted alleles, runs trait-specific SuSiE and multi-signal coloc, and retains every failed/underpowered locus in a final evidence ledger. All per-signal coloc outputs remain available. A maximum-H4 summary is a descriptive convenience, not a multiplicity-corrected shared-mechanism claim.

### Other adapters and evidence integration

- `munge`, `h2`, `rg`: original LDSC CLI, no forced cross-trait intercept, physical/QC checks.
- `lava`: native per-locus LAVA with sample-overlap inputs; low local h² and execution errors remain different states. No reuse of rejected old LAVA results.
- `genomicsem`: genuine multivariate LDSC S/V/I and reviewed `usermodel`; never a pairwise-rg matrix masquerading as sampling covariance.
- `mr`: native TwoSampleMR IVW, weighted median and MR-Egger, heterogeneity and leave-one-out. Instruments must be selected on the **exposure**, not PLACO or outcome significance.
- `pleiofdr`: isolated MATLAB launcher around official code. Already aligned, checksum-bound reference/trait MAT exports and a MATLAB license are required. Its output remains review-required, not automatically accepted.
- `integrate`, `cell-enrichment`, `cross-disorder`: combine **normalized, sourced evidence tables**; matched independent-locus competitive enrichment; non-duplicated locus counts. These do not create raw single-cell sequencing data or download every QTL resource automatically.
- `forest`, `heatmap`, `manhattan`, `qq`, `report`: figures and receipts. Manhattan display binning never changes the inference dataset.

## Workflow engines

The Python DAG runner works without Snakemake. An optional separate Snakefile exposes the same tasks:

```bash
python -m pip install -e '.[workflow]'
snakemake --snakefile workflow/Snakefile --cores 1 \
  --resources production_slots=1 \
  --config runtime=/absolute/path/preparation.runtime.json
```

Snakemake was not installed or executed in this delivery environment. The provided GitHub workflow adds Python testing and a native R syntax-check job when run by your repository; no successful GitHub Actions run is being claimed here.

## Why a run can stop

`BLOCKED_BY_SOFTWARE`, missing data, mismatched builds, ambiguous alleles, insufficient LD coverage, no eligible primary pair, insufficient power, failed numerical methods and excessive failure rates are meaningful outcomes. A stopped or negative result is not fixed by lowering thresholds after inspecting it. Keep the original artifact and make any revised analysis a separately reviewed run.

**The final scientific deliverable is a complete evidence ledger and defensible conclusions—not six guaranteed mechanisms.**
