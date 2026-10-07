# LD and native fine-mapping reproduction

Run from the repository root on the new research branch. The historical
source/result directories are read only. Current runtime: conda Python3.13
with NumPy2.4.3/SciPy1.18.1; R4.6.1 with isolated exact CRAN package sources.
Inspect `R_installed_packages.tsv` and `R_runtime_sessionInfo.txt` for the
full installed lock. No package was installed in the user's global R library.

## Small self-contained LD rebuild: no source download

The three versioned `chr*_genotype_factors.npz` files contain exact ALT
0/1/2 dosages, sample/population identities and variant order. Full matrices
are ignored local intermediates, about210MiB combined. Each full float64
matrix requires at most182MiB and the sequential rebuild stays within the
initial8GiB RAM budget.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/rebuild_signed_LD_from_factors.py
```

This independently rebuilds all three full matrices and requires exact
SHA256 agreement with the original empirical reference outputs. Results are
`factor_rebuild_validation.json` and the ignored local matrices under
`work/ld_genotypes_research_v1/`. NumPy/BLAS version can affect final-bit
hashes; a different build must also compare scientific matrix values and
record a new environment, never claim exact hash agreement silently.

## Repeat original genotype-source acquisition if needed

Inspect `LD_PRE_OUTCOME_PROTOCOL.md`, `LD_SOURCE_MANIFEST.json`, source reuse
policies, memory/disk availability and exact byte-range plan first. Current
original bytes/header segments total26,443,591 plus754,117metadata bytes.
The acquisition script verifies exact206ranges and refuses full-file fallback.
The original retained v5b TBIs were identical to the new official TBIs.

```sh
/opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/acquire_regional_ld.py metadata
/opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/acquire_regional_ld.py download
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/reconstruct_signed_ld.py
/opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/additional_prior_lead_comparison.py
```

The additional comparison requires its recorded primary Ensembl mapping
JSONs, acquired under the separate published-lead protocol; those remote
mapping receipts include exact source URLs and SHA256. No expanded genotype
acquisition is required. Initial reconstruction code content is retained as
`reconstruct_signed_ld_initial_executed.py` and matches its original execution
hash. The current reproducer redirects large matrices to ignored local work;
all scientific calculations are unchanged.

## Native scientific runtime

The exact source-package acquisition manifests bind30package names/versions,
CRANMD5, sourceSHA256 and byte counts. Source tarballs total34,293,244bytes,
and the final isolated library is64MiB (150MiB library cap). The initial
incorrect archive URLs and binary4.6metadata404 are retained as failed
access attempts; they were not substituted for runtime validation.

```sh
/opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/acquire_R_dependencies.py
Rscript brain6/translational_psychiatry_research_v1/research_v1/ld/install_exact_R_sources.R
```

Do not query a live index and silently adopt new versions. The source/package
manifest is the scientific runtime lock; a changed source/runtime requires a
new receipt and version before a new experiment. The source-index preparation
helper exists for a fresh explicit runtime experiment, not implicit upgrades.

## Native chr5 experiment: source inputs must remain local

Use the parent discovery acquisition/extraction script and original archive
manifests to restore exact insomnia/ADHD sources first. Their hashes and
regional hashes are in `discovery/acquisition_receipts.json`. ADHD's original
README prohibits public redistribution of its result file: raw regional files,
source-derived signedZ/beta/SE diagnostics and nativeRDS objects remain
ignored local inputs. Do not put them in a public branch or deliverableZIP.
Publication-style aggregate/posterior outputs do not substitute for external
investigator permission requirements or participant-level access.

```sh
OPENBLAS_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/prepare_finemap_inputs.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript brain6/translational_psychiatry_research_v1/research_v1/ld/native_RSS_consistency.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript brain6/translational_psychiatry_research_v1/research_v1/ld/native_finemap_coloc.R
Rscript brain6/translational_psychiatry_research_v1/research_v1/ld/export_native_BFs.R
/opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/ld/validate_native_coloc.py
```

The first stage fails checksum/source gates rather than accepting incomplete
slices. The RSS stage uses the exact pinned API adapter and requires finite
complete diagnostics. The fit stage checks all frozen source/consistency
gates and executes all8fixedfits/12coloc/3ABFsensitivities. The final Python
calculation independently recomputes H0–H4 from native Bayes factors; it is
numerical validation of the posterior equation, not independent disease
replication or validation of the RSS model assumptions.

Native raw recovery paths are `work/ld_genotypes_research_v1/C_finemap_inputs/`:
`<trait>_<model>_native_fit.rds`, `<model>_p12_<prior>_coloc_susie.rds`,
`same_universe_p12_<prior>_coloc_abf.rds`, compact`<trait>_RSS_consistency.rds`
and local signedLD/beta/SE/Z tables. Their hashes are in
`native_analysis_provenance.json`. CompactRSSRDS retain every kriging table
entry and critical scalar; matrix factors regenerate the eigenanalysis.

## Regression and empirical validation remain separate

```sh
/opt/miniconda3/bin/python -m unittest discover -s brain6/translational_psychiatry_research_v1/research_v1/ld -p 'test_ld_research.py' -v
Rscript brain6/translational_psychiatry_research_v1/research_v1/ld/test_native_R_adapter.R
```

Ten source-free regressions and three native synthetic API checks passed.
Actual genotype reconstruction, native GWAS fits, source identity and
independent posterior recomputation are separate numerical/scientific
records; successful tests do not establish a new psychiatric mechanism.
