# Frozen validation protocols

Each new experiment has a separate protocol and run identity. Historical thresholds and outcomes are unchanged.

- LD: `research_v1/ld/LD_PRE_OUTCOME_PROTOCOL.md` and acquisition manifests lock three GRCh37 regions, 503 European samples, biallelic nonpalindromic variants, MAF ≥0.01, source/pedigree checks and numerical gates. No clipping or PSD projection was permitted.
- Native chr5: `FINEMAP_PRE_FIT_PROTOCOL.md` locks the full chr5:103447968–104447968 window, source/N/frequency/coverage gates, strong association in both traits, RSS consistency and fixed model/prior sensitivities. Eight SuSiE fits, 12 coloc-SuSiE comparisons and three same-universe ABF checks were executed. Posteriors are exploratory; the independent four-hypothesis family remains separate.
- Clinical fixed variants: `sources/FIXED_VARIANT_FINGEN_PROTOCOL.md` and its manifest were frozen before three public API outcomes. Four pair slots retain alpha 0.05/4. Conjunction P is the maximum of the two source P values. Missing fields remain NOT_ESTIMATED; no favorable SNP search or dense-region inference was admitted.
- Calibration: `statistics/calibration_protocol.json` and `calibration_freeze.json` lock 24 scenarios ×10,000 truth-known Gaussian-summary replicates with exact AR(1) LD. These are method diagnostics, not empirical GWAS validation.
- Global contrasts: `global_contrasts/protocol.json` locks all 72 retrospective comparisons, conservative sampling-covariance bounds, Bonferroni 72 and marginal h² QC.
- Molecular: `molecular/sensitivity_protocol.json` locks four retained contexts ×five shared-variant priors ×three QTL sdY scales. A separate reviewer protocol locks the 60-condition nonpalindromic stress analysis.

Schema/version/serialization failures and amendments are retained. Corrections changed field access, API handling or serialization, not hypotheses, priors, variant rules, gates or denominators. The original frozen files remain available.
