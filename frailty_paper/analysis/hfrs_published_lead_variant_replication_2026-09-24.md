# HFRS published lead-variant replication cross-reference

Generated: 2026-09-26T10:08:07+00:00

## Source and method

This descriptive audit reads Supplementary Tables 1 and 2 from the publisher workbook cited by Mak et al., Nature Aging 5:1589–1600 (2025), doi:10.1038/s43587-025-00925-y. The article says these tables contain FinnGen variants crossing P<5×10⁻⁸ and reports independent leads at r²<0.01; it reports checking those variants in UK Biobank. UKB result availability here requires non-missing beta, SE, and P cells. Nominal significance is P<0.05; genome-wide significance is P<5×10⁻⁸, matching the article's reported cutoffs. The source SHA-256 is `123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625`. Per-variant rows are in `hfrs_published_lead_variant_replication_2026-09-24.tsv`; machine-readable provenance is in `hfrs_published_lead_variant_replication_2026-09-24.json`.

## Results

| Table | FinnGen hits | UKB available hits | All-hit P<0.05 / P<5×10⁻⁸ | Independent leads | UKB available leads | Lead P<0.05 / P<5×10⁻⁸ |
|---|---:|---:|---:|---:|---:|---:|
| ST1, HFRS | 1588 | 1261 (article 1262) | 688 / 73 (article 688 / 73) | 53 | 36 | 13 / 2 (article 14 / 2) |
| ST2, HFRS without dementia | 492 | 435 (article 435) | 118 / 21 (article 118 / 21) | 42 | 26 | 10 / 1 (article 10 / 1) |

The public ST1 workbook has 1,261/1,588 hits with UKB beta, SE, and P present, compared with 1,262 reported in the article. Its all-hit nominal and genome-wide p-value counts match (688 and 73). Among lead rows, the workbook has 13/36 with P<0.05, compared with 14/36 in the article, and both report 2 with P<5×10⁻⁸. ST2 all-hit and lead availability and p-value counts match the article. These one-row differences are preserved as public-workbook-to-article discrepancies. No p-value was reconstructed or recategorized to force agreement. Article-reported counts are at lines 108–110 of the published Results.

## Interpretation boundary

This reproduces a subset of the Mak et al. source study's own FinnGen-to-UKB lookup. It does not establish independent replication of any sleep–frailty relationship in this project: UKB is shared by several project-side sources, and exact pairwise participant intersections remain unknown. The workbook is hit-conditioned and is not a substitute for full HFRS summary statistics; do not use it for LDSC, LAVA, or genome-wide sleep–HFRS tests. HFRS plan decision D10 remains open pending exact endpoint and authorized source verification. No effect-direction comparison is made because this audit does not independently verify the effect-allele convention for both beta columns.
