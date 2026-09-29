# Exploratory branch runtime

Observed on 2026-09-28: Python 3.13.5; `pysam==0.23.3`; `pyliftover==0.4.1`; R 4.6.1; `coloc==5.2.3`. The temporary Python packages were installed under `/private/tmp/brain6-pysam` and the R package under `/private/tmp/brain6-r-lib`; these temporary folders are not part of the source receipt. The scripts require equivalent pinned packages when rerun.

Run order from repository root, with `/Volumes/Extreme SSD/brain6-work/brain6-exploratory-functional-v1/` writable and official EBI regional access available:

```sh
PYTHONPATH=/private/tmp/brain6-pysam python3 brain6/scripts/brain6_exploratory_functional_v1_extract_qtl.py
PYTHONPATH=/private/tmp/brain6-pysam python3 brain6/scripts/brain6_exploratory_functional_v1_extract_gwas.py
python3 brain6/scripts/brain6_exploratory_functional_v1_prepare_coloc.py
BRAIN6_R_LIB=/private/tmp/brain6-r-lib Rscript brain6/scripts/brain6_exploratory_functional_v1_run_coloc.R
python3 brain6/scripts/brain6_exploratory_functional_v1_build_status.py
```

All scripts refuse to overwrite completed source slices or result tables. The first two scripts preserve per-slice hash receipts and validate existing slices on resume. To reproduce in a new environment, point the script constants to an isolated new versioned output directory, while retaining the same frozen source files, chain and protocol.
