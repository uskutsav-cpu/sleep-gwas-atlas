# Acquisition manifest re-verification — 2026-09-24 13:58 UTC

Fresh read-only verification against the mounted external acquisition volume. Branch `frailty-paper-v1`; source HEAD `0b034f65ce339bfabbe0ff5e34ee30f62bc8032c`. The mixed worktree contains unrelated Brain6/root changes, left outside this check.

## Command

```sh
FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' \
  make -C frailty_paper \
  PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python \
  verify-manifest audit-manifest-metadata
```

## Result

The verifier passed all **351/351** registered resources by path, byte count and SHA-256. The metadata audit passed all **351 rows** and **22 required columns**. Captured command output: `acquisition_manifest_reverification_2026-09-24_1358.log` (SHA-256 `dc1b2932200affe6eb9660924e7afcd759576cb7c62454e8c4f2985030d9589e`).

| Input | SHA-256 |
|---|---|
| `manifests/all_acquired_resources.tsv` | `cb7e26c1192ae73c52acf45ce69f79a8a8aff928ab365fa919ecc8cb71894685` |
| `scripts/16_verify_acquired_resource_manifest.py` | `4be937ef1aca90c877d0aa93eb76649d5692b9f8545a0cdf58dc5875a84d51ed` |
| `scripts/45_audit_acquisition_manifest_metadata.py` | `4c03f2c50098bd3e7507806800aa15c82a242ceaeb94387c0416ce9f0adeedb1` |

This establishes current registered-file identity and metadata integrity only. It does not establish phenotype eligibility, valid downstream analysis, or completion of manual review steps.
