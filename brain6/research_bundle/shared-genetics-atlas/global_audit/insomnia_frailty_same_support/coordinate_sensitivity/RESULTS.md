# HFRS coordinate audit and shared-support sensitivity

Using the project's GRCh37 EUR LD-score positions, the UCSC hg19-to-hg38 chain and the original HFRS GWAS's GRCh38 records, **1,134,416 of 1,134,953 selected markers** have a unique exact lifted chromosome, position, rsID and allele-pair match. For 420 markers the chain gave no unique mapping under the check; for 117 it returned an alternate-contig coordinate rather than the source's primary-chromosome position. These 537 are unresolved by this mapping procedure, not proven invalid. The original GWAS source passed its publisher MD5, and the [complete list](../hfrs_coordinate_unresolved.tsv), [audit receipt](../hfrs_coordinate_receipt.json), and [protocol](../HFRS_COORDINATE_PROTOCOL.md) retain status and examples. The check does not validate reference-sequence bases, LD-model fit or participant independence.

As a post-hoc sensitivity, the same 537 rsIDs were removed from **all four** trait inputs before rerunning LDSC. All three regressions then used 1,134,416 markers. The point estimates changed little at the log's four-decimal precision:

| Insomnia comparison | Original r_g (SE) | Coordinate-resolved r_g (SE) |
|---|---:|---:|
| FI | 0.4263 (0.0389) | 0.4261 (0.0389) |
| FFS | 0.2445 (0.0371) | 0.2444 (0.0369) |
| HFRS | 0.6171 (0.0326) | 0.6168 (0.0324) |

This checks stability to one explicit support restriction; it does not establish equivalence with a confidence interval, independently replicate the results or address the HFRS elevated intercept and shared FinnGen recruitment. The excluded markers need not be false or unimportant—one has |HFRS Z|=5.02. The paired jackknife was subsequently recomputed on the restricted support, as shown below. Both FI sleep/fatigue content and FFS exhaustion remain. The [protocol](PROTOCOL.md), `support_receipt.json`, `verification.json`, scripts, restricted-input hashes and full logs document the calculation.

## Paired contrasts after coordinate restriction

The 200 delete-block calculations passed the shared-insomnia-vector and per-pair SE checks. Restricted-support differences and paired SEs remain close to the earlier estimates:

| Contrast | Original difference (paired SE) | Restricted difference (paired SE) | Restricted Bonferroni p, 3 contrasts |
|---|---:|---:|---:|
| HFRS − FI | 0.1908 (0.0419) | 0.1907 (0.0421) | 1.76×10⁻⁵ |
| HFRS − FFS | 0.3726 (0.0424) | 0.3724 (0.0422) | 3.53×10⁻¹⁸ |
| FI − FFS | 0.1818 (0.0406) | 0.1817 (0.0405) | 2.22×10⁻⁵ |

The restricted point estimates are rounded in the LDSC logs, so differences and p-values are approximate. The correction covers this post-hoc three-contrast family; it is not a correction for all analyses in the atlas. These same-sample, same-model comparisons do not resolve HFRS's elevated intercept, population differences or shared ascertainment. Full results, unadjusted intervals and reproducibility hashes are in `paired_contrasts.tsv` and `paired_receipt.json`.
