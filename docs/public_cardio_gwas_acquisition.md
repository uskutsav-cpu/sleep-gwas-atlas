# Public cardiovascular GWAS acquisition

This log records the locally exercised cardiovascular releases selected for
`atlas-v1.0`. Raw and derived data remain ignored and are not redistributed.
The eleven-trait checkpoint is explicitly partial and does not replace the
locked 45-trait, 396-pair analysis.

| Trait | Selected release | Verified file evidence | Local checkpoint |
|---|---|---|---|
| Systolic blood pressure | [Evangelou 2018](https://doi.org/10.1038/s41588-018-0205-x), OpenGWAS `ieu-b-38` | The 216,210,855-byte GWAS-VCF archive has SHA-256 `3cdead05ac3daa88265c8fd8fa53f5f89e4f0e8945dfc6b2f8fc5a3bafaa260c`. Named VCF fields produced 7,088,067 rows from 7,088,083 records after excluding 16 missing rsIDs. | Harmonization retained 5,964,514 rows; 1,107,872 munged effects. Observed-scale h² 0.1325 (SE 0.0048), intercept 1.1535. |
| Atrial fibrillation | [Nielsen 2018](https://doi.org/10.1038/s41588-018-0133-9), `GCST006061` | The 190,312,717-byte archive has SHA-256 `3ef0f55a29ba0c065df14e55cfb4fb4351ac4da90e0da28e3438426e63e4879f`; its exact member produced 12,149,979 rows. | Harmonization retained 10,246,131 rows; 1,211,889 munged effects. With the cited adult `K=0.03`, liability h² 0.3160 (SE 0.0350), intercept 1.0123. |
| Any stroke | [Mishra 2022 GIGASTROKE](https://doi.org/10.1038/s41586-022-05165-3), `GCST90104539` European stratum | The direct 283,517,495-byte GRCh37 gzip has SHA-256 `15384ea2ff7e03faaf66a2a0c315fd7d2f71c0d39f76906278131d6db695c624` and 7,511,476 literal rows. Its schema has coordinates and alleles but no rsID or INFO. | The pinned coordinate-plus-allele map retained 1,176,288 unambiguous HapMap3 rows with zero identity conflicts; all remained after munging. With `K=0.03`, liability h² 0.0967 (SE 0.0080), intercept 1.0731. |

## Stroke liability and correlation checkpoint

The prevalence convention is deliberately explicit. CDC BRFSS surveillance
estimated age-standardized self-reported ever-stroke prevalence at 2.9% among
US adults in 2020–2022 ([Imoisili et al. 2024](https://doi.org/10.15585/mmwr.mm7320a1)).
The manifest rounds this to `K=0.03` as an adult-population approximation; it
does not present that value as age-specific European or global lifetime risk.

The source carries no rsID, so the production path assigns identities only
when GRCh37 coordinate and unordered allele pair agree with the
checksum-pinned EUR HapMap3 map. The resulting 1,176,288 effects passed the h²
gate (Z 12.09). In the 18-pair partial checkpoint, snoring×stroke has r_g
0.0598 (SE 0.0358, P 0.0953, partial-family FDR 0.1410) and does not pass.
Insomnia×stroke has r_g 0.1382 (SE 0.0414, P 0.0008, FDR 0.002400),
insomnia×SBP has r_g 0.0572 (SE 0.0208, P 0.0059, FDR 0.01328), and
insomnia×atrial-fibrillation has r_g 0.0646 (SE 0.0288, P 0.0249, FDR
0.04075); all three pass the partial-family threshold. None substitutes for
the final 396-pair FDR.
