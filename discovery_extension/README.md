# Novelty-Enriched Phenome Discovery Extension

This namespace is scientifically and operationally separate from the locked
`atlas-v1.0` core. The core contains 45 traits and the prespecified 396
sleep-by-non-sleep tests. Nothing under this directory changes the core panel,
its GWAS identities, its QC thresholds, or its multiple-testing family.

The extension is exploratory and will be described as:

> The prespecified 45-trait core atlas was complemented by a separately
> defined novelty-enriched phenome-wide discovery extension.

## Immutable core boundary

`core_checkpoint.json` records the producing commit plus byte-level SHA-256
hashes for the locked panel, source/schema manifests, Phase-1 QC tables, the
complete 396-pair result family, and its figure. Run the checkpoint verifier
before and after every extension stage:

```bash
python3 discovery_extension/scripts/00_verify_core_checkpoint.py
```

The verifier fails closed if any pinned core artifact changes, if the current
commit no longer descends from the checkpoint commit, or if the core result is
not the exact 12 x 33 Cartesian product.

## Namespace contract

Extension-only configuration, manifests, logs, results, figures, provenance,
and acceptance gates live below this directory. Synthetic fixtures, if ever
introduced, must live under `synthetic/` and must never be copied into real
results.

