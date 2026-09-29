# Supplementary material — provisional

This supplement accompanies a provisional analysis report. It contains only
tables supported by acquired, source-linked data available at this audit. It is
not a complete or frozen supplement: licensed literature exports and screening
are missing, additional frailty definitions have not passed their source/QC
gates, and no frailty-specific independent replication or downstream locus
analysis is available.

## Available tables

- **Supplementary Table 1 — GWAS metadata.** Thirty-two GWAS source records,
  with accessions, sample descriptions, ancestry/build, cohort and overlap
  notes, effect encoding, file size, checksum and access/license notes. This
  table includes candidate resources that are acquired but not analysis-
  eligible; source presence must not be read as phenotype validation.
- **Supplementary Table 4 — SNP heritability.** Twenty current estimates for
  the 12 sleep/circadian traits, primary FI and seven latent frailty factors.
  This table does not include eligible physical-component, Fried or HFRS h²
  estimates.
- **Supplementary Table 5 — available global genetic correlations.** One
  combined, family-labeled table of 168 estimates: 12 frozen-atlas sleep × FI
  pairs, 84 secondary sleep × latent-factor pairs and 72 read-only frozen-atlas
  sleep × aging-context pairs. The original q-values and denominators are
  retained. The inherited FI/aging family has denominator 396; the latent
  sensitivity family has denominator 84. The rows are not pooled or compared
  across families.
- **Supplementary Table 7 — latent-frailty results.** The 84 secondary LDSC
  estimates with family q-values, overlap status and claim limits. These are
  sensitivity analyses, not independent replication. Pair-level columns flag
  the insomnia/general-factor part-whole dependence (insomnia is an input
  indicator) and the related but nonidentical tiredness/lethargy indicator for
  daytime sleepiness. These source-measure notes do not resolve exact cohort
  intersections.
- **Supplementary Table 8 — replication-resource audit.** Candidate-source,
  phenotype-match, access and pairwise eligibility assessments. No audited
  candidate currently establishes an independent sleep–frailty rg pair.
- **Supplementary Table 15 — multiple-testing sensitivity (partial).** The
  table combines 12 frozen sleep × FI results corrected with same-family
  all-396 Bonferroni and 84 latent-factor sensitivity results corrected with
  same-family 84-pair Bonferroni, retaining their primary BH q-values. All nine
  FI pairs passing inherited all-396 BH also pass all-396 Bonferroni. In the
  latent-factor family, 48/84 pass BH and 37/84 pass Bonferroni. The latent
  family remains sensitivity-only with exact participant overlap unknown.
  A companion 119-row conclusion-by-sensitivity matrix reports each current
  headline conclusion against the planned Phase 16 checks and distinguishes
  completed sensitivities from partial convergence, unresolved overlap,
  blocked sources, and analyses not justified by upstream gates. It is a
  status synthesis, not a set of new estimates; unavailable checks are not
  null results. Alternate frailty definitions, component estimates,
  cohort/ancestry sensitivities and locus-level sensitivity tests remain
  unavailable, unresolved or not justified; this is not the full sensitivity
  suite. See
  [`table_s15_sensitivity_conclusion_matrix.tsv`](supplementary_tables/table_s15_sensitivity_conclusion_matrix.tsv).
- **Supplementary Table 16 — software, resources and versions (partial).**
  Pinned tool policy is paired with observed software versions from the FI h²
  reproduction, LDSC logs, figure provenance and renderer records. Pinned but
  unused/not-evidenced tools and plotting-version mismatches are labeled.
  Gated analyses and final commands are not covered yet.

Machine-readable tables are in `supplementary_tables/`. Their row counts,
status and file names are in `supplementary_tables/supplementary_table_status.tsv`;
the accompanying provenance JSON hashes every input and output. Working table
legends are in `table_legends.md`.

## Tables not available or not justified

Supplementary Tables 2–3 cannot be populated: the 56,117-record PubMed queue
has zero screening decisions, no studies have been included or appraised, and
licensed Embase, Scopus, Web of Science and optional PsycINFO exports are
absent. Supplementary Table 6 is blocked because the five physical-component
candidate sources lack adequate phenotype/build/effect/model provenance. Fried
and HFRS estimates are unavailable. Supplementary Tables 9–14 (LAVA,
PLACO/shared loci, fine-mapping, colocalization, molecular-QTL and cell-type
evidence) are not justified without eligible frailty-specific independent
replication and locus evidence. Supplementary Table 16 remains partial because
software and run-version details are distributed across manifests, provenance
files and logs; the current Table 16 is a partial version inventory and must be
reconciled with final commands before the package freeze.

These unavailable tables are documented in the status index, not represented
by empty estimate rows or inferred null findings. The supplement must be
revised after the relevant source, review and analysis gates are completed; no
current output is a final result freeze.
