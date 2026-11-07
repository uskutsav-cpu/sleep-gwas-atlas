# Reproducibility of key external inputs

The [curated external inventory](external_inputs.tsv) records **102 files** used by major analyses in this atlas. At this check, all 102 were readable on the mounted volumes and matched the SHA-256 recorded in earlier receipts. The inventory covers:

- 58 processed GWAS, LD-score and LDSC code files for the two main genome-wide rechecks;
- 32 CDC NHANES XPT files used by the provisional phenotype pilot;
- raw archived Brain6 insomnia and ADHD, Atkins FI and Ye FFS copies;
- the archived ADHD-derived munged file, four three-frailty common-support inputs, two brain-cell QTL matching tables, and the prior UKB HFRS candidate extraction.

This is a **curated list, not all dependencies of all analyses**. The [machine-readable status](external_check.json) and [audit code](extended_audit.py) record exact absolute paths, roles, expected and observed hashes, receipt links and caveats. Hashes authenticate equality to the earlier locally recorded copies, not publisher origin, correct phenotype coding, cohort independence, statistical calibration or clinical validity. In particular, the archived ADHD-derived file's name does not authenticate its unlocated raw input or ingest process. The current FI copy differs from the unavailable legacy frozen-atlas FI hash. The two brain QTL tables are processed subsets, not independent source authentication or a study-matched LD panel.

The files reside on the external SSD and mounted disk images; they are **not included in the atlas bundle**. To repeat this check after mounting the same volumes, from the workspace root run:

```sh
python3 outputs/shared-genetics-atlas/data_manifest/extended_audit.py
```

The command reads and hashes the inputs and rewrites only `external_inputs.tsv` and `external_check.json` in the atlas. It reports missing, unreadable and mismatched files but does not reacquire them. A clean result does not ensure every R/Python package, system library, reference index, remote source, model output or internet endpoint needed by each analysis is present. Follow each linked analysis protocol; several original run scripts use fixed output prefixes and can overwrite prior outputs if run unchanged. Select a new output path for reruns.

Important evidence gaps remain: no full per-SNP sleep/fatigue-only-excluded FI GWAS was identified in the checked public sources; the published broader mental-wellbeing-excluded FI files were not found in the checked release; the ADHD-derived raw source and ingest receipt remain unverified; and suitable study-matched LD is not established for the requested full local and multi-signal analyses. These are not missing-file assertions about every possible source. The [selected atlas manifest](../manifest.tsv) and [bundle integrity manifest](../BUNDLE_INTEGRITY.md) cover different sets of files and should be retained together with this inventory.
