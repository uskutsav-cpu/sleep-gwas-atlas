# Physical-function lookup at prior chromosome 5 sleep–ADHD variants

A partial indexed query of the [harmonized physical-function GWAS](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/harmonised/) found one exact GRCh38 coordinate, rsID and allele-pair match at each previously prioritized variant. Higher baseline scores represent better physical function, and each row reports N = 405,979.

| Variant (GRCh38) | Effect / other allele | Beta (SE) | Source p |
|---|---|---:|---:|
| rs2431108, chr5:104612267 | C / T | −0.004367 (0.001120) | 9.61×10⁻⁵ |
| rs77960, chr5:104628884 | A / G | −0.004396 (0.001121) | 8.86×10⁻⁵ |

Both rows are nominally associated, but neither reaches the conventional genome-wide threshold 5×10⁻⁸. These variants were selected after prior sleep–ADHD analyses and may be correlated, so a two-row numerical correction would not cover the project's selection or wider testing. Their signs concern the reported effect alleles and physical-function score only. A marginal association does not establish that the same variant causes both physical-function and sleep/ADHD signals, nor does it validate AlphaGenome predictions or identify a causal gene. The physical-function measure is not a Fried score without exhaustion, and overlap with other UK Biobank sources is possible.

The [protocol](PROTOCOL.md), [query](query.py), and [receipt](receipt.json) preserve the selected rows, reference alleles and query. This read only portions of the harmonized GRCh38 file; its whole-file checksum was not verified. The separate original GRCh37 file was fully verified for the [genome-wide comparison](../RESULTS.md). A regional colocalization analysis would require complete allele-harmonized local data, suitable priors, signal/QC checks and careful overlap assessment; these two p-values are not colocalization.

An [original-release cross-check](original_check/RESULTS.md) confirmed that both selected rows use opposite effect alleles between the two releases, with consistent beta sign flips and complementary frequencies.
