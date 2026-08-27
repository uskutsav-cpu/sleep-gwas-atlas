# Public immune GWAS acquisition

This log records the locally exercised immune release selected for
`atlas-v1.0`. Raw and derived data remain ignored and are not redistributed.
The ten-trait, 16-pair checkpoint is explicitly partial and does not
replace the locked 45-trait, 396-pair analysis.

| Trait | Selected release | Verified file evidence | Local checkpoint |
|---|---|---|---|
| Crohn disease | [de Lange 2017 IIBDGC](https://doi.org/10.1038/ng.3760), `GCST004132` European Crohn stratum | The direct 312,502,048-byte GRCh37 gzip has SHA-256 `ddc4258b51f6da59749cbb798d4e387b5d1fab1746f3101d31c7bfa72301c872` and 9,570,787 literal rows for 12,194 cases and 28,072 controls. The source has coordinate-encoded `MarkerName`, alleles, beta, SE, and P, but no rsID, INFO, or frequency. | The pinned coordinate-plus-allele map retained 1,144,234 unambiguous HapMap3 rows; all remained after munging. With `K=0.003`, liability h² is 0.2508 (SE 0.0266), with Z 9.43 and intercept 1.1144. |

## Crohn disease identity and liability checkpoint

The verified archive is hard-linked to `data/raw/crohn.txt.gz`, avoiding a
second local copy. Because the source has no rsID, the production path assigns
identities only when the GRCh37 coordinate and unordered allele pair agree
with the checksum-pinned EUR HapMap3 map. Of 9,570,787 source rows, 8,057,115
were outside that restricted map and 369,438 candidate rows had a conflicting
allele pair. The 1,144,234 retained effects had zero coordinate conflicts or
ambiguous multi-rsID mappings. The QC ledger explicitly records that source
INFO and frequency fields are absent rather than silently treating them as
available filters.

The source counts imply constant effective N 34,004.9. For liability scaling,
the manifest uses `K=0.003` as a rounded high-prevalence European approximation.
This is anchored to the population-based South Limburg estimate of 331 per
100,000 discussed by [Carbonnel and Boutron
2017](https://doi.org/10.1093/ecco-jcc/jjx072) and the underlying [South Limburg
study](https://doi.org/10.1093/ecco-jcc/jjx055); it is not represented as a
global prevalence.

The 1,144,234-SNP input passed the predefined h² gate. In the 16-pair
partial checkpoint, snoring×Crohn disease has r_g -0.0322 (SE 0.0295,
P 0.2751, partial-family FDR 0.3144) and does not pass. Insomnia×Crohn disease
has r_g -0.0708 (SE 0.0312, P 0.0231, FDR 0.03984) and passes this partial
family, but is not a substitute for the final 396-pair FDR. The Crohn, combined
IBD, and ulcerative-colitis strata from this study are not independent, so
future cross-trait interpretation must account for their sample overlap.
