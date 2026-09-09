# Downstream analysis contract

The downstream policy freezes every remaining atlas-v1.0 layer before its
results are available. It does not convert missing results into evidence.
Fine-mapping begins only at primary cross-method shared loci, uses signed
European GRCh37 LD, permits multiple signals with SuSiE-RSS, and requires
convergence plus LD/sumstat diagnostics. Colocalization uses signal-specific
`coloc.susie`, a locked prior-sensitivity grid, and distinguishes shared from
distinct signal models.

Molecular, regulatory, cell-type, and pathway evidence must retain exact
dataset releases, local snapshots, checksums, biological context, and explicit
negative or access-blocked searches. No landing page, nearest-gene assignment,
TWAS association, or regulatory overlap is treated as causal evidence.

Causal inference is last, bidirectional, and requires strong instruments,
heterogeneity/pleiotropy/directionality checks, and compatibility with
colocalization before any causal wording. The robustness matrix covers every
predeclared family, including alternate LD, boundaries, priors, models,
cohorts, GWAS, multiplicity, SEM leave-one-out, and local-correlation choices.

The two upstream-complete canonical tables can already be generated without
implying downstream completion:

```bash
python3 scripts/51_build_atlas_core.py
python3 scripts/51_build_atlas_core.py --validate-only
```

This publishes only `traits.tsv`, `trait_pairs.tsv`, and their checksum-bound
provenance. The integrated-atlas acceptance gate remains blocked until every
other layer has exact terminal coverage and every canonical evidence row, if
any, is validated. A genuinely empty evidence table is permitted only when its
complete locked task family records a real null or justified non-applicability.
