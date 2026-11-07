# Physical-block sensitivity of conditional enrichment

All 1,138,262 common-support SNPs had a unique rsID and allele-pair match in the pinned GRCh37 HapMap3 map. For the post-hoc comparison of insomnia p≤0.01 before and after conditioning on frailty p≤0.001, delete-one-block jackknife estimates were:

| Conditioning trait | Increase in insomnia-p≤0.01 fraction | 5 Mb blocks: descriptive 1.96-SE interval | 10 Mb blocks: descriptive 1.96-SE interval |
|---|---:|---:|---:|
| FI | 3.80 percentage points | 1.68–5.92 points | 1.75–5.86 points |
| FFS | 3.30 percentage points | 1.35–5.25 points | 1.33–5.27 points |

The 5 Mb partition has 562 occupied physical windows; the 10 Mb partition has 290. The results are similar across these two block sizes, which is a useful sensitivity to clustering but not proof that all LD or sources of dependence have been captured. These are normal approximations for the descriptive difference, not calibrated tests of genetic sharing or false-discovery control. Physical blocks vary in size and marker density and may have long-range LD across boundaries; the thresholds were selected after reviewing data. Sample overlap, phenotype ascertainment, FI sleep items, rounded Z and other GWAS artifacts remain possible contributors. The fractions are not evidence for a specific variant or causal pathway, nor a published-style conjFDR analysis.

The reverse conditioning direction, standard errors, block counts and the largest single-block contribution are in `summary.tsv`; all underlying per-block numerator and denominator counts are in `block_counts.tsv`. A separate `verify.py` reproduced all eight differences and jackknife standard errors from those counts and recorded hashes in `verification.json`. The [protocol](PROTOCOL.md) specifies the join and calculation; [the original descriptive comparison](../RESULTS.md) retains the full thresholds and simplified Gaussian reference controls.
