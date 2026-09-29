# Clean committed-snapshot preflight — 2026-09-24

## Scope

Validated the committed tree at `e2c8c7cbdc75b7345dca4335420351897e2d5b95` (`Record PubMed reporting audit commit hash`) in an isolated directory at `/private/tmp/frailty-clean-e2c8c7c-20260924`. The snapshot was assembled from Git archives of the committed `frailty_paper/`, `frailty_analysis/`, `environment/`, `config/`, `results/`, `scripts/` and `data/` paths. It did not use uncommitted repository files. The isolated Python 3.11.11 environment was `/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03/work/conda-envs/frailty-paper-py311-clean-2026-09-23`.

Runtime inputs unavailable in the archive were mounted read-only from `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1`: raw, harmonized, munged and frailty-paper data. The 311 PubMed XML files were hard-linked into the temporary snapshot only after each file's size and SHA-256 matched the committed `review/pubmed_source_files.tsv` manifest. The blank review outputs were rebuilt there; the resulting output-manifest SHA-256 matched the committed source manifest (`7301215e3d9fa40cfd479b0afbc4cc73b3cd3e8c1c4c227bda22e6055480379c`) and its content matched the live checkout.

## Results

- Acquired-resource integrity: 351/351 files verified; metadata audit: 351 rows and 22 required columns.
- Frozen analysis-plan validator: pass (12 sleep traits; 396-pair multiplicity family).
- Cohort-overlap validator: pass (34 rows); exact participant intersections remain unknown.
- Review structure: pass (56,092 title/abstract records; zero decisions; no full-text records).
- Bibliography, missing-abstract audit, manuscript quantitative-claim audit (28/28), and reporting locators (92 rows, 93 line references, zero errors): pass.
- Latent-factor Q_SNP audit: 32,035,589 source rows, zero significant Q_SNP variants, zero excluded rows; audit passed for seven factors.
- Unit suite: 112/112 passed under Python 3.11.11.

Full command output: [`preflight_clean_committed_snapshot_py311_2026-09-24.log`](preflight_clean_committed_snapshot_py311_2026-09-24.log), SHA-256 `0b83ac7c929f541b21836b92dc9a421969e72301188413acfdd6a27a85aa7c64`.

## Limits

This is a clean committed-code snapshot with explicitly mounted external inputs, not an end-to-end scientific replay. It does not resolve unscreened records, licensed database searches, unavailable alternate frailty GWAS, exact cohort overlap, replication, or the failed/incomplete LAVA family gate. Passing software and provenance checks do not establish scientific readiness or change those statuses.
