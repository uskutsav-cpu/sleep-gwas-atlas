# Figure 2 provenance refresh and visual QA

**Review time:** 2026-09-24 05:11 UTC

**Branch / starting HEAD:** `frailty-paper-v1` / `145dd14500365accab3dd86de94c34e74a4be040`

## Provenance refresh

The prior Figure 2 provenance hash for `frailty_paper/analysis/replication_resource_audit.tsv` did not match the current committed source table. The plot was regenerated from the current source inputs with the pinned Python 3.11.11 environment:

```text
work/conda-envs/frailty-py311/bin/python frailty_paper/scripts/38_plot_global_rg_matrix.py \
  --repo /Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03 \
  --outdir frailty_paper/analysis --replace-outputs
```

The script reports 96 estimate pairs, with 9/12 primary FI pairs passing the inherited all-396 BH family, all 9/12 also passing the same-family Bonferroni sensitivity, and 48/84 latent-factor pairs passing the separate fixed 84-pair BH family. The latent-factor family remains sensitivity-only. No estimates, correction families, or interpretation rules were changed.

After regeneration, all 11 provenance checks pass: six input hashes, three rendered-output hashes, the source-data TSV hash, and the plotting-script hash. The source-data file contains 96 unique sleep-by-frailty estimates and four summary rows. Current input state records no eligible independent pairwise sleep–frailty replication.

## PDF and visual inspection

`figure2_global_rg_matrix.pdf` is a one-page vector PDF (1425 × 967.5 pt); its 4000-pixel PNG and SVG are retained alongside it. The PDF was rendered to a 2850 × 1935 PNG at 144 dpi and visually reviewed. The title, axis and trait labels, cell values, within-family significance markers, two distinct correction-family notes, replication status, cohort-overlap caveat, and non-causal interpretation are legible and unclipped. No layout defect was observed.

This is visual and provenance QA for a provisional global-rg figure, not a final-paper or independent-replication claim. The figure correctly labels FI values as frozen-atlas reuse and latent-factor results as sensitivity-only; local sharing, component decomposition, and mechanistic evidence are not implied.
