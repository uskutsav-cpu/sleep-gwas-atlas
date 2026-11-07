# Physical-function lookup at two prior Brain6 variants

A partial indexed query of the public [harmonized GCST90565837 release](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/harmonised/) found one exact GRCh38 coordinate, rsID and allele-pair match for each preselected chromosome 3 variant. The trait is the objective baseline physical-function composite, where higher values mean better function; N = 405,979 was recorded in both rows.

| Prior variant (GRCh38) | Effect / other allele | Beta (SE) | Source p |
|---|---|---:|---:|
| rs62264764, chr3:117920728 | A / G | 0.0016133 (0.00147953) | 0.27553 |
| rs62264767, chr3:117923158 | C / A | 0.00161314 (0.00147961) | 0.27560 |

Neither row reaches even nominal p < 0.05 in this selected lookup. This is not evidence of no effect, and does not establish shared causal variants, sleep–frailty replication, or gene mediation. The two variants were selected previously and may be correlated; a much wider project search has been performed. This is a single-variant check, not a regional colocalization or a test of all variants in the locus. The UK Biobank-based physical-function phenotype is not a matched exhaustion-deleted frailty score and may overlap other cohorts.

The [protocol](PROTOCOL.md), [query](query.py), and [receipt](receipt.json) retain the full rows, index hash, release URL and software version. The indexed query read only selected portions of the harmonized GRCh38 release, so its publisher's whole-file MD5 was **not** verified here. The earlier whole-file checksum relates to the separate original GRCh37 release. Harmonization fields are preserved without assigning undocumented semantics. Reproduction requires an HTSlib-enabled `pysam`, network access and the saved index; there was no use of AlphaGenome output as training data.
