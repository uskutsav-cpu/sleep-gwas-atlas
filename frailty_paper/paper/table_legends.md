# Supplementary-table legends — working package

The available tables below are source-linked drafts, not a complete frozen
supplement. `supplementary_tables/supplementary_table_status.tsv` tracks all 16
planned tables and identifies which outputs are unavailable or not justified.

## Supplementary Table 1. GWAS metadata

One row per each of the 32 GWAS summary-statistics resources in the verified
acquisition manifest: 12 sleep/circadian traits, the primary FI, six aging or
adjacent-context traits, five physical-component candidates and eight Catalog
accessions (seven latent frailty factors plus pneumonia). Fields retain source
accession, publication, sample description, ancestry, build, effect encoding,
cohorts, overlap notes, file size, SHA-256 and license/access notes. A registered
source is not thereby considered scientifically eligible. See
`supplementary_tables/table_s1_gwas_metadata.tsv`.

## Supplementary Table 2. Systematic-review studies — not available

No studies have been included because the 56,117-record PubMed queue is
unscreened, licensed database exports are absent, and full-text screening has
not begun. No included-study table is generated. Update only after the required
database exports have been imported, screened and adjudicated.

## Supplementary Table 3. Risk-of-bias results — not available

No included studies have been appraised. The JBI tool registry and blank
item-level templates are scaffolds, not risk-of-bias results. Do not infer low
risk from missing appraisals.

## Supplementary Table 4. SNP heritability

Twenty current h² rows cover the 12 sleep/circadian traits, primary FI and seven
latent frailty factors. This is a partial table: it does not contain eligible
physical-component, Fried or HFRS h² analyses. Retain source-specific scale,
diagnostics and gate status. See
`supplementary_tables/table_s4_snp_heritability.tsv`.

## Supplementary Table 5. Available global genetic correlations

The 168 currently available rows comprise 12 frozen-atlas sleep × FI pairs, 84
secondary sleep × latent-factor pairs and 72 read-only frozen-atlas sleep ×
aging-context pairs. The table preserves each correction family and its
denominator. FI q-values inherit the all-396 BH family; latent q-values use the
fixed 84-pair secondary BH family; aging-context q-values inherit the all-396
family. Latent results are sensitivity-only, and aging-context rows are not
replication. No physical-component, Fried or HFRS results are included. See
`supplementary_tables/table_s5_global_rg.tsv`.

## Supplementary Table 6. Physical-component results — blocked

Five candidate component GWAS files are structurally readable but lack adequate
deposit-level phenotype/build/effect/model provenance. No sleep-pair estimates
are generated or represented as zero. Reconsider only after source eligibility
and all frozen QC gates pass.

## Supplementary Table 7. Latent-frailty-factor results

The 84 sleep × seven latent-factor LDSC estimates, including SE, Z, P, q, family
denominator, cohort-overlap status and claim limits. These are secondary
sensitivity results and not independent replication. See
`supplementary_tables/table_s7_latent_factor_rg.tsv`.

## Supplementary Table 8. Replication-resource audit

Source-level candidate audit and pairwise eligibility assessment. Candidate
frailty-side datasets do not establish an eligible independent sleep–frailty
rg pair; unknown participant overlap is not treated as zero. The HFRS source
row distinguishes its published, hit-conditioned FinnGen-to-UKB lead lookup
from project-level sleep–frailty replication and records two ST1 workbook-to-article
count differences. See
`supplementary_tables/table_s8_replication_audit.tsv`.

## Supplementary Tables 9–14. Downstream genetic and molecular evidence — not justified

No frailty-specific LAVA family, shared-locus set, fine-mapping, trait–trait
colocalization, locus-specific molecular-QTL or cell-type analysis is available.
The corresponding tables are not generated because upstream independent
replication/locus evidence and eligible inputs are absent. Existing atlas/Track
B results are not substituted for frailty-specific results.

## Supplementary Table 15. Multiple-testing sensitivity (partial)

Reports 12 primary sleep × FI comparisons with the inherited all-396 BH values
and Bonferroni adjustment over the same 396 tests, plus 84 secondary sleep ×
latent-factor comparisons with their prespecified 84-pair BH q-values and a
same-family 84-pair Bonferroni sensitivity. Nine of 12 FI pairs retain
significance under both corrections; 48 of 84 latent-factor pairs pass BH and
37 pass Bonferroni. BH remains the primary correction for both families. The
latent estimates remain sensitivity-only because exact participant overlap is
unknown; this table does not establish independent replication. It covers
multiple-testing correction only, not the full planned suite of alternate
phenotype, cohort-overlap, ancestry, LD, or locus-level sensitivities. See
`supplementary_tables/table_s15_correction_sensitivity.tsv` and the
119-row conclusion-by-sensitivity status matrix in
`supplementary_tables/table_s15_sensitivity_conclusion_matrix.tsv`. The
matrix records blocked and untested domains explicitly; it does not treat
them as negative results or create new estimates.

## Supplementary Table 16. Software, resources and versions — partial

The current table combines pinned tool-version policy with observed versions
from the FI h² reproduction, LDSC logs, figure provenance and renderer
manifests. It distinguishes pinned-but-not-evidenced tools from tools recorded
in current runs and explicitly labels plotting-version mismatches. The table is
not a complete software audit for analyses that remain gated and must be
reconciled with final commands before freeze. See
`supplementary_tables/table_s16_software_resources_versions.tsv`.
