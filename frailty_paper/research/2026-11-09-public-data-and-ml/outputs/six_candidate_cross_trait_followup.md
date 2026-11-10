# Cross-trait follow-up for six expanded-screen candidates

The six candidates from the exploratory MVP-insomnia/HFRS screen were checked in a second clinical insomnia GWAS, FinnGen R13, and aligned across the project's FI, FFS, UKB-insomnia and HFRS sources. **rs817509**, in addition to the previously highlighted rs12536335, has the same estimated sign and unadjusted p<0.05 in all six source results. This is a descriptive, post-selection pattern, not six independent replications or proof of a causal shared variant.

| Variant; positional label | FI p | FFS p | UKB insomnia p | MVP insomnia p | FinnGen insomnia p | FinnGen HFRS p |
|---|---:|---:|---:|---:|---:|---:|
| rs12536335; FOXP2 | 2.63 × 10⁻⁶ | 0.0017 | 7.65 × 10⁻⁷ | 4.57 × 10⁻⁸ | 0.00625 | 7.85 × 10⁻⁹ |
| rs817509; RP11-384F7.1 | 0.0223 | 0.0099 | 3.78 × 10⁻⁷ | 2.86 × 10⁻⁵ | 0.00272 | 1.21 × 10⁻⁷ |
| rs2155292; NCAM1 | 2.76 × 10⁻⁸ | 0.54 | 7.25 × 10⁻⁵ | 3.89 × 10⁻⁶ | 0.00116 | 1.10 × 10⁻⁹ |
| rs28528997; MAD1L1 | 8.04 × 10⁻⁵ | 2.3 × 10⁻¹² | 0.000292 | 3.24 × 10⁻⁷ | 0.484 | 2.16 × 10⁻⁶ |
| rs2815764; NEGR1 | 0.000136 | 1.2 × 10⁻⁵ | 0.00511 | 1.13 × 10⁻⁷ | 0.667 | 2.30 × 10⁻⁶ |
| rs1909650; PCDH8P1 | 6.32 × 10⁻⁷ | 0.097 | 9.02 × 10⁻¹⁰ | 1.58 × 10⁻⁵ | 0.102 | 0.000386 |

All six have matching point-estimate signs across the six results after alignment to the HFRS effect allele, but weak estimates can change sign and nominal p-values are not multiple-testing-adjusted. The FFS effect-allele convention is inferred, and five of these rsIDs have additional alleles in broader reference data. The FI and FFS GWAS use overlapping UK Biobank resources, and FinnGen R13 insomnia and the FinnGen HFRS GWAS may share participants. The results from the same ancestry or cohort should not be counted as independent evidence merely because they concern different outcomes.

For rs817509 G versus A, the Finnish insomnia log-odds beta was **+0.0249** (SE 0.00832; approximate 95% CI +0.0086 to +0.0412), while its HFRS beta was **+0.0829** (SE 0.0157) on the continuous hospital score scale. They concern different outcomes and units; their magnitudes cannot be compared or added. FinnGen's raw insomnia annotation names `AC068633.1`, while the earlier GRCh37 Ensembl nearest-body lookup names `RP11-384F7.1`; these are source-specific annotations, not proof of a functional target. No effort was made here to establish whether the identifiers refer to the same transcript or locus. The marker has a multiallelic reference rsID and requires variant-representation checks and regional fine-mapping before biological interpretation.

rs28528997 and rs2815764 have nominal evidence in UKB and MVP insomnia, but their FinnGen insomnia intervals include zero (p=0.484 and 0.667). This does not demonstrate no Finnish effect or a difference between populations: the estimates require a direct, appropriately specified heterogeneity comparison, and the source definitions and ascertainment differ. rs2155292 has a strong HFRS association but weak FFS evidence (p=0.54). The FI meta-analysis includes sleep/tiredness-related items; FFS contains exhaustion; HFRS can include broad codes for drowsiness or exhaustion. No comparison removes that overlap.

## Methods and traceability

The complete checksum-verified FinnGen R13 insomnia file was scanned: 21,327,010 rows, source SHA-256 `353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049`. All six candidates had exactly one row matching rsID, mapped GRCh38 position, and the specific allele pair. R13 ALT beta was flipped where necessary to match the HFRS effect allele. The lookup and approximate confidence intervals are in `six_joint_hits_finngen_insomnia_followup.csv` (SHA-256 `939cbe1afb3db91357767e175b03fb898237d195a6681b7b56e6fa51f146bb2e`). Source p-values were retained; no inferential test or score was refit. The final aligned six-trait matrix is `six_candidates_six_trait_direction_matrix.csv` (SHA-256 `6098541a0393c14aeea769b495533f48d1d47fa10b2e798a49d3e18eb4f24fa9`). Source alleles and the original expansion criteria are in `song124_mvp_hfrs_followup.md`.

The earlier candidate selection, comparison across many markers, overlapping cohorts, distinct phenotypes, unverified Song build/allele semantics and limited source QC prevent confirmatory claims. Proximity or consistent signs do not establish a common causal variant or gene; independent replication, LD-aware fine-mapping and colocalization are needed. No individual-level outcomes or validated 15–30-marker risk score were obtained in this follow-up.
