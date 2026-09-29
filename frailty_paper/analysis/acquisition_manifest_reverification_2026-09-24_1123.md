# Acquisition manifest re-verification — 2026-09-24 11:23 UTC

This is a fresh read-only integrity check against the mounted acquisition volume. It was run on branch `frailty-paper-v1` at source HEAD `9fd48e9bfb4ffde75eaa72e969532197d32e7e51`, with the repository worktree still mixed by unrelated Brain6 and root workflow/environment edits.

## Command

```sh
FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' \
  make -C frailty_paper \
  PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python \
  verify-manifest audit-manifest-metadata
```

## Result

```text
RESOURCE_MANIFEST_OK rows=351 files_verified=351
MANIFEST_METADATA_OK rows=351 required_columns=22
```

The verifier checked each manifest path, byte count and SHA-256 against the configured external storage root. The metadata auditor checked all required manifest columns and row constraints. The exact manifest and tool hashes at execution were:

| Input | SHA-256 |
|---|---|
| `manifests/all_acquired_resources.tsv` | `cb7e26c1192ae73c52acf45ce69f79a8a8aff928ab365fa919ecc8cb71894685` |
| `scripts/16_verify_acquired_resource_manifest.py` | `4be937ef1aca90c877d0aa93eb76649d5692b9f8545a0cdf58dc5875a84d51ed` |
| `scripts/45_audit_acquisition_manifest_metadata.py` | `4c03f2c50098bd3e7507806800aa15c82a242ceaeb94387c0416ce9f0adeedb1` |

This verifies registered file identity, integrity and required metadata at this snapshot. It does not establish that any GWAS passes phenotype/source eligibility, nor that downstream scientific analyses are valid.
