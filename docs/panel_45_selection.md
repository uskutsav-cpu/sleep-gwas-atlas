# Final 45-trait panel: selection contract

`config/panel_45_selection.tsv` records how the final Phase 0/1 traits were
selected. The authoritative production contract is now
`config/analysis_panel.tsv`, validated against
`config/analysis_panel.lock.json`. It has exactly 45 distinct trait IDs:

- 12 sleep/circadian traits, 6 neuro/psychiatric traits, 6 immune traits,
  9 metabolic/cardiovascular traits, 6 aging traits, and 6 cancer traits.
- 24 sleep/circadian, aging, and cancer rows originated in the Codex curation stream.
- 21 neuro/psychiatric, immune, metabolic, and cardiovascular rows originated
  in the parallel curation stream.

Selection means a trait is in scope; it does **not** mean it is source-verified,
harmonization-ready, LDSC-ready, or eligible for Phase 1. The 86-row historic
registry is retained only as candidate-selection provenance. Production scripts
must not read it.

The panel deliberately prefers common, interpretable disorders and traits
(for example Alzheimer’s disease, major depression, asthma, type 2 diabetes,
coronary artery disease, stroke, breast cancer, and sleep duration) over
duplicated phenotype definitions or poorly powered placeholders. Every source
still has to meet the existing EUR, GRCh37/hg19, public-access, header, and
metadata gates before the derived readiness stages can advance.
