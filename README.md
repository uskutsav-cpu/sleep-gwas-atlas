# sleep-gwas-atlas

An auditable Phase 0/1 pipeline for the Sleep/Circadian Genetic Atlas:

- **Phase 0:** source-GWAS curation, hg19/EUR validation, harmonization, and HapMap3 munging;
- **Phase 1:** LDSC SNP heritability, power/confounding QC, sleep-by-disease genetic correlations, FDR correction, and a provenance-stamped heatmap.

The registry currently contains **86 candidate traits** (20 sleep/circadian and 66 disease/trait). It is intentionally not described as 90: no extra traits have been invented simply to meet an earlier documentation claim.

## Current scientific status

The code path has a synthetic end-to-end smoke test. **No real GWAS summary statistics, EUR LD reference panel, or real LDSC results are stored in this repository.** Every registry row is currently `TODO`, so the empirical Phase 0/1 analysis is not yet complete.

That distinction is deliberate. The pipeline will not:

- treat an unsourced population prevalence as a liability-scale assumption;
- silently lift hg38 sources to hg19;
- convert a curation status into an h2 QC result; or
- turn synthetic logs into a presentable figure or report.

## Set up the real analysis environment

This downloads the maintained Python 3 LDSC implementation, Python dependencies, and EUR reference data. It requires network access, several GB of storage, and should be started only after reviewing the source/data-access plan.

```bash
bash scripts/00_setup.sh
export PYTHON_BIN=.venv/bin/python
export LDSC_PYTHON=.venv/bin/python
export LDSC_DIR=ldsc
```

The setup follows the Python 3 `ldsc39` branch from CBIIT; the legacy upstream LDSC repository itself now points users to that maintained implementation.

## Phase 0: curate before downloading

`config/traits.tsv` is the pipeline registry. A trait should be changed to `CURATED` only after its source file and metadata are verified. For each binary trait, record a citation for `pop_prev` in a `pop_prev_citation` column before requesting liability-scale h2. Keep the source's exact ancestry subset and genome build in the registry as well.

```bash
$PYTHON_BIN scripts/10_phase0_audit.py \
  --out results/tables/phase0_curation_audit.tsv

# After a curated hg19/EUR raw file is placed in data/raw/:
bash scripts/02_munge.sh insomnia mdd
```

The harmonizer writes one `data/harmonized/<trait>.qc.txt` ledger per trait. It requires rsIDs, autosomal hg19 coordinates, effect alleles, effect size, SE, P, and an auditable sample-size rule; it removes strand-ambiguous SNPs, MHC, low INFO/MAF variants when present, duplicates, and low effective-N variants.

## Phase 1: h2 gate, rg, and Figure 2

```bash
# Final h2 for binary traits requires a cited pop_prev_citation.
bash scripts/03_h2_qc.sh insomnia mdd

# An observed-scale binary h2 is allowed only for interim QC / rg preparation.
bash scripts/03_h2_qc.sh --observed-scale insomnia mdd

# Uses only traits that are both CURATED and h2 PASS; writes the inclusion table.
bash scripts/04_rg.sh --h2 results/tables/h2_summary.tsv
```

`05_collate.py` keeps low h2 Z and high LDSC intercept as separate QC failures. The gate is h2 Z >= 4 and intercept <= 1.20. `04_rg.sh` creates `results/tables/phase1_inclusion.tsv` rather than asking anyone to overwrite curation state with `PASS`.

## Reproducible smoke test

After installing the Python requirements (or with any environment that already has NumPy, Pandas, and Matplotlib):

```bash
PYTHON_BIN=.venv/bin/python bash scripts/run_smoke_test.sh
```

The test accepts tabular and whitespace-delimited raw inputs, BETA and OR effect formats, runs h2 and rg parsing across the complete registry, exports the proposal metadata table, generates the report and heatmap, and watermarks all fake output under `results/_smoketest/`.

## Outputs

```text
config/traits.tsv                         registry and curation state
scripts/10_phase0_audit.py                source-metadata readiness audit
data/harmonized/<trait>.qc.txt            filter-by-filter QC ledger
data/munged/<trait>.sumstats.gz           HapMap3 LDSC input (not committed)
results/tables/h2_summary.tsv             h2, Z, intercept, QC reason
results/tables/phase1_inclusion.tsv       CURATED + h2-PASS selection record
results/tables/rg_matrix.tsv              rg, SE, P, Benjamini-Hochberg FDR
results/figures/fig2_rg_heatmap.png       provenance-stamped Phase 1 heatmap
```

Generated raw data, intermediate files, logs, tables, figures, and smoke-test artifacts are ignored by Git so an orphan result cannot accidentally become versioned evidence.
