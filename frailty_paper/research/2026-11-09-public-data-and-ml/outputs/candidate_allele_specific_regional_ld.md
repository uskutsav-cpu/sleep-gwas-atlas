# Allele-specific regional LD among selected candidates

I used the same 503 European 1000 Genomes phase 3 reference individuals as the chromosome 5 analysis, reading exact positions from the indexed public GRCh37 VCF. Each pair below had the candidate's listed two-allele contrast in that reference VCF; all 1,006 haplotypes had phased calls. Correlation r is between the listed allele indicators, not between an individual's disease risks.

| First variant and counted allele | Second variant and counted allele | Reference r | Reference r² |
|---|---|---:|---:|
| rs113851554 T | rs78206187 G | +0.270 | 0.073 |
| rs35881094 G | rs17550821 A | −0.366 | 0.134 |
| rs11152363 A | rs1942262 A | +0.383 | 0.146 |
| rs11152363 A | FFS lead rs660010 G | +0.090 | 0.0081 |
| rs1942262 A | FFS lead rs660010 G | +0.050 | 0.0025 |

The TCF4-region pairs with the published lead rs660010 previously returned empty Ensembl LD API results. This **direct reference VCF calculation** gives low pairwise r² for the specific contrasts in this panel, rather than treating the empty results as zero. It does not establish conditional independence in the GWAS, prove the candidate and lead signals are distinct, or rule out other linked variants. The lead G allele matches Ye et al.'s clumping table's A1, whose effect convention remains inferred.

The within-candidate estimates differ somewhat from the earlier Ensembl rsID-level values (for example, rs113851554–rs78206187: 0.073 versus 0.095). The exact reason is not established; the services may use different reference/sample or calculation conventions. Keep both source-specific values rather than treating them as interchangeable. The other two pairs show modest LD, so the initial 15-marker set should not be assumed independent. Signed correlations are relative to the counted alleles in the table; a negative r does not imply that either allele is protective.

Reference VCF matching does not verify the original study coding or resolve all multiallelic rsID issues. The candidates are post-selected and the results are not a shared-signal test, fine-mapping or a validated risk score. The per-pair summaries, source coordinates, allele metadata and failures are in `candidate_regional_allele_specific_ld.json`; code is `work/reference_candidate_ld.py`. No individual genotypes were saved. See `chr5_allele_specific_ld_review.md` for the separate chromosome 5 result.
