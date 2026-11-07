# Resampling the smaller LD reference for three selected pairs

For the three largest r² differences identified in the preceding post-selected comparison, resampling the 503 European-labeled 1000 Genomes participants produced the following within-panel percentile intervals:

| GRCh38 variant pair | 1000 Genomes r² | ADSP reported r² | 1000 Genomes bootstrap 2.5–97.5% |
|---|---:|---:|---:|
| chr3:117835650 T>G / 117858571 G>A | 0.706 | 0.434 | 0.525–0.858 |
| chr3:117835650 T>G / 117890455 T>G | 0.672 | 0.432 | 0.472–0.836 |
| chr3:117520255 T>C / 117698054 C>T | 0.417 | 0.188 | 0.101–0.700 |

For the first two pairs, the ADSP value falls below the small-reference bootstrap interval; the third falls inside it. This shows that simple individual resampling within the smaller panel does not encompass the first two observed differences under this particular calculation. It does **not** establish statistical significance or an ancestry effect: these were the three largest differences selected after inspection, reference participants may be related or structured, the larger panel also has sampling uncertainty, and the publisher's LD statistic and the dosage-correlation r² may not be identical estimands. Both panels differ from the disease and QTL samples. No full-matrix or causal conclusion follows.

We retrieved the same 503 European-labeled samples and GRCh37 variant records from the [1000 Genomes reference](https://www.internationalgenome.org/data) and verified that the reconstructed dosage r² matched the saved matrix. Each row summarizes 5,000 jointly resampled individuals with seed 20261028; no replicate was monomorphic. The fractions of resamples at or below the respective ADSP values are 0.0024, 0.0104 and 0.0932, but they are **descriptive resampling summaries, not p-values** for a selected between-panel comparison.

A later [within-European population check](POPULATION_RESULTS.md) found low allele counts and many monomorphic resamples in some small groups, underscoring how unstable subgroup r² can be. It does not establish the reason for the ADSP difference.

[Protocol](BOOTSTRAP_PROTOCOL.md), [code](bootstrap.py), [complete results](bootstrap.tsv) and [receipt](bootstrap_receipt.json) record the procedure. The original panel identifier file and marker-table hashes matched prior receipts; the remotely queried VCF and index were not whole-file hashed. No AlphaGenome output was used. This does not establish the cause of LD–Z residuals or repair failed LAVA QC.
