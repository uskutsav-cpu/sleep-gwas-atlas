# Observed marginal Z-score tails in the PLACO+ input

The observed Z-score distributions in the parameter-estimation inputs have more large marginal values than a **single centered normal curve with the same fitted variance**. For example, among approximately 1.14 million retained rows, the FFS marginal has **783** values with |Z|>5; that single-normal reference gives an expected count of about **51**. The corresponding FI comparison is **306 vs 28**, and the insomnia marginal has approximately **112–118 vs 11**, depending on the paired subset.

| Pair and marginal | Retained rows | Observed |Z|>5 | Single-normal reference count |
|---|---:|---:|---:|
| Insomnia–FI: insomnia | 1,138,254 | 118 | 10.6 |
| Insomnia–FI: FI | 1,138,254 | 306 | 27.9 |
| Insomnia–FFS: insomnia | 1,138,229 | 112 | 10.5 |
| Insomnia–FFS: FFS | 1,138,229 | 783 | 50.9 |

**This is not a type-I-error test.** The observed inputs contain genuine single- and potentially two-trait associations, correlated SNPs, and a selection rule excluding doubly strong associations. The fitted variance comes from those same data. PLACO+ specifies a mixture of null components; it does not require the entire selected marginal to be one centered normal distribution. Consequently, excess counts do not identify false positives, invalidate PLACO+, estimate a false-discovery rate or validate any of the 82 diagnostic crossings. They do show why checking the true one-trait-effect distribution and study-specific sampling assumptions matters before promoting those crossings.

The [protocol](PROTOCOL.md), [all thresholds](marginal_tails.tsv), [analysis](analyze.py), [receipt](receipt.json), and [independent full rescan](verification.json) record counts for |Z|>3, 4, 5 and 6 and input hashes. The data are rounded Z-scores on restricted HapMap3 support; rare expected counts and LD prevent treating count ratios as independent-SNP significance tests. The separate [model-consistent null calculation](../model_null/RESULTS.md) and [fixed-mean stress test](../RESULTS.md) remain conditional on their distinct assumptions. No candidate or failed LAVA gate was reclassified.
