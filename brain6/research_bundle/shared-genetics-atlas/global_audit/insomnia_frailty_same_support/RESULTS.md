# Registry-insomnia correlations with FI, FFS and HFRS on shared SNP support

A new LDSC execution compared FinnGen R13 registry insomnia with three frailty definitions on the same 1,134,953 input markers. The FI and FFS estimates reproduce the saved native same-support results to the logged precision. The HFRS association is a same-source sensitivity; HFRS and insomnia both use FinnGen, so they may share participants.

| Frailty outcome | r_g (SE) | p | Bonferroni p across these 3 tests |
|---|---:|---:|---:|
| FI | 0.4263 (0.0389) | 5.64×10⁻²⁸ | 1.69×10⁻²⁷ |
| FFS | 0.2445 (0.0371) | 4.47×10⁻¹¹ | 1.34×10⁻¹⁰ |
| HFRS | 0.6171 (0.0326) | 7.80×10⁻⁸⁰ | 2.34×10⁻⁷⁹ |

All four inputs have identical ordered SNP and allele keys, their three frailty hashes match the archived support record, and each regression retains 1,134,953 valid-allele SNPs. The S-file hash is recorded, but the support record did not independently specify it. The correction shown is for this bounded, post-hoc three-pair family, not the full 16- or 396-test project families.

For comparisons, 200 delete-block r_g estimates were reconstructed from covariance and both heritability deletion vectors. The shared insomnia deletion vectors were exactly equal across all three pairs, and the FI/FFS deletion estimates matched the earlier native outputs within 10⁻¹⁰. Paired jackknife differences and the three-contrast Bonferroni correction are:

| Contrast | Difference | Paired SE | Unadjusted 95% CI | Bonferroni p across 3 contrasts |
|---|---:|---:|---:|---:|
| HFRS − FI | 0.1908 | 0.0419 | 0.1086 to 0.2730 | 1.61×10⁻⁵ |
| HFRS − FFS | 0.3726 | 0.0424 | 0.2895 to 0.4557 | 4.62×10⁻¹⁸ |
| FI − FFS | 0.1818 | 0.0406 | 0.1023 to 0.2614 | 2.25×10⁻⁵ |

These are post-hoc statistical comparisons, with the HFRS point estimate rounded to four decimals. The confidence intervals are not multiplicity-adjusted. They are conditional on the shared support, LDSC model and input choices; they do not show that one frailty construct is superior, that effects are causal, or that differences persist in independent samples. In particular HFRS has an elevated univariate LDSC intercept **1.1877 (SE 0.0144)** and the R13–HFRS genetic-covariance intercept is **0.1790 (SE 0.0072)**. Their FinnGen overlap, general-European LD reference, HFRS ascertainment and use of total rather than verified per-variant N warrant caution. The intercept alone cannot identify the precise cause or shared-participant count. FI includes insomnia and tiredness; FFS includes exhaustion. This is not a sleep-item-exclusion analysis or a phenotypic correlation.

`correlations.tsv`, `paired_contrasts.tsv`, `receipt.json`, `paired_receipt.json`, the complete logs, delete values and scripts provide the calculations. The [single-pair protocol](PROTOCOL.md) and [paired protocol](PAIRED_PROTOCOL.md) record the analysis sequence and limitations. The earlier [UKB-only FI sensitivity](../../frailty_sensitivity/FRAILTY_COMPARISON.md) answers a different comparison.

## Input-lineage follow-up

The S input was subsequently checked row by row as an ordered exact subset of the previously pinned R13 insomnia file: 1,134,953 rows match SNP, alleles, N and Z, with both gzip files read to completion and the parent SHA-256 matching its frozen receipt. FI, FFS and HFRS are likewise ordered exact numeric/allele subsets of the archived three-frailty support, whose own and upstream hashes match its previous receipt. These [lineage checks](insomnia_lineage_receipt.json) and `frailty_lineage_receipt.json` strengthen the file-identity audit; they do not validate source phenotype definitions, underlying HFRS harmonization, per-variant sample sizes, sample independence or reference ancestry fit.

A subsequent [alignment receipt](hfrs_alignment_receipt.json) independently compares all 1,137,927 HFRS records in the parent three-frailty support with the pinned pre-alignment munged file. Every record uses the same ordered allele pair, signed Z and N as its upstream record; none required a swap or strand complement in this particular intermediate step. Because the final 1,134,953 HFRS input is an exact subset of that parent, this rules out a sign change introduced during these two subset/merge steps. It does **not** establish that the original raw HFRS effect allele or the earlier LDSC munging was correct, or that per-variant N, ancestry or sample overlap are resolved. `verify_hfrs_alignment.py` reproduces the check.

The subsequent [filtered-source munging audit](hfrs_munging_receipt.json) independently reconstructed all 1,134,953 selected HFRS signed Z scores from the filtered source p-values and beta directions. Every allele pair was directly aligned and the largest difference from the stored Z was 0.000500, consistent with three-decimal rounding. N values also matched the filtered file. This validates the filtered-to-final transformation, not the original source effect-allele semantics, participant overlap, actual per-variant sample sizes or ancestry fit.

## Original HFRS GWAS source check

The complete [GCST90707291 HFRS discovery source](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90707001-GCST90708000/GCST90707291/) was scanned and its MD5 matched the frozen catalog metadata. Of 1,161,680 filtered HFRS rows, 1,161,676 had a unique full-string rsID, chromosome and nonpalindromic-allele match with the original source and matching beta, SE, p and effect-allele frequency after alignment. The other four rsIDs appeared as individual tokens in comma-separated source identifiers; an [alias check](hfrs_alias_receipt.json) found a unique matching source row and matching values for each. None of those four entered the 1,134,953-marker comparison.

The original file labels `effect_allele` and `other_allele`; this confirms that the project's filtered beta and frequency agree with that **published source convention** after any allele reversal. Because the filtered table does not contain base-pair positions, this audit is not an independent GRCh38-to-GRCh37 liftover check. The supplied N=500,737 is the study-level discovery total and is not confirmed per-variant N. Exact sample overlap, cohort and model assumptions, and LD-reference fit remain open. Source bytes, counters and examples are recorded in [the full-scan receipt](hfrs_raw_alignment_receipt.json); [protocol](HFRS_RAW_PROTOCOL.md) and scripts preserve the checks.

A [coordinate audit and support sensitivity](coordinate_sensitivity/RESULTS.md) matched 1,134,416 of the 1,134,953 HFRS markers to original GWAS coordinates after liftover. Excluding the 537 unresolved markers from every trait changed each reported r_g by at most 0.0003; this does not resolve ancestry, overlap or source-model uncertainty.

A [covariance-intercept stress test](intercept_sensitivity/RESULTS.md) shows that fixed values around the fitted intercept change the conditional r_g estimate. These runs are not a confidence interval or an ancestry correction; the unconstrained result remains primary. The FinnGen workflow lists Finnish LD scores, but its specified cloud resource was inaccessible under the current environment’s access policy.
