# Frozen discovery-extension audit

Read-only, checksum-based review. No extension analysis was rerun.

## Core checkpoint

Checkpoint commit: `22df92c8d6894f93c50561ae8a86acf2f48f50cb`; ancestor of current `HEAD`: **True**.
Pinned artifacts verified under their exact or documented append-only rules: **8**.
Pinned artifacts missing: **10**; mismatches: **0**.
The core checkpoint is **not fully re-verifiable from this checkout** because required files are absent. This qualifies the older extension status report, which recorded the checkpoint as passing at generation time.

Missing pinned core artifacts:
- `results/tables/analysis_panel_provenance.tsv`
- `results/tables/source_readiness.tsv`
- `results/tables/phase0_readiness.tsv`
- `results/tables/trait_readiness.tsv`
- `results/tables/h2_summary.tsv`
- `results/tables/phase1_inclusion.tsv`
- `results/tables/rg_primary_phase1.tsv`
- `results/tables/rg_qc_failed_sensitivity.tsv`
- `results/tables/phase1_summary_report.md`
- `results/figures/fig2_rg_heatmap.png`

Pinned core mismatches:
- None

## Frozen extension outputs

Reported outputs matching their provenance hashes: **9/9**.
Result row counts: 100 h2 traits; 1200 rg pairs; 217 replication candidates; 23 top replicated discoveries.
Frozen report records 603 extension-FDR-significant pairs among 1200 planned tests.
Extension readiness: `GLOBAL_DISCOVERY_AND_REPLICATION_COMPLETE_DOWNSTREAM_BLOCKED`.

Checksums:

- `counts` — `discovery_extension/results/final_extension_counts.tsv`: **PASS**
- `replication_figure_caption` — `discovery_extension/figures/extension_replication_summary.caption.txt`: **PASS**
- `replication_figure_pdf` — `discovery_extension/figures/extension_replication_summary.pdf`: **PASS**
- `replication_figure_png` — `discovery_extension/figures/extension_replication_summary.png`: **PASS**
- `report` — `discovery_extension/final_report.md`: **PASS**
- `top_novel_discoveries` — `discovery_extension/results/top_novel_discoveries.tsv`: **PASS**
- `100-trait rerun h2` — `discovery_extension/results/ldsc/extension_trait_readiness.tsv`: **PASS**
- `1200-pair rg family` — `discovery_extension/results/ldsc/extension_rg_matrix.tsv`: **PASS**
- `217-pair replication family` — `discovery_extension/results/replication/replication_results.tsv`: **PASS**

## Reuse boundary

The existing extension reports qualified global discovery and independent replication. Local correlation, pleiotropy, fine-mapping, colocalization, mechanism and causal interpretation remain blocked upstream; the extension report encodes these as `NA_BLOCKED_UPSTREAM`, not as zero findings. This exploratory 100-trait family remains separate from the 45-trait core atlas and from the current frailty analysis.
