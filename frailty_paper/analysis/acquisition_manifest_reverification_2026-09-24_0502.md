# Current acquisition-manifest re-verification

**Run time:** 2026-09-24 05:02 UTC

**Branch / starting HEAD:** `frailty-paper-v1` / `77dc9c47468e554af7a588815ec1aabafe3335da`

**Current manifest:** `frailty_paper/manifests/all_acquired_resources.tsv`

**Manifest SHA-256:** `3724ee5e878a5d6b7cb46f27c5c8ebb3f74229898ee659de23f62f0f3c2f8fba`

## Commands and results

```text
make -C frailty_paper PYTHON=.venv/bin/python audit-manifest-metadata
MANIFEST_METADATA_OK rows=351 required_columns=22

make -C frailty_paper PYTHON=.venv/bin/python verify-manifest \
  FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'
RESOURCE_MANIFEST_OK rows=351 files_verified=351
```

At the run, the internal volume had 7.8 GiB free and the external analysis volume had 1.7 TiB free. The external data volume is available; new large derived outputs should still use it because internal free space remains below the repository's 20 GiB gate.

This verifies the current worktree manifest's required metadata columns and each registered resource's recorded path, byte size, and SHA-256 against the configured external root. It does not validate source eligibility, phenotype equivalence, exact cohort overlap, or statistical suitability. The manifest itself is modified in the mixed worktree and was not committed by this verification record.
