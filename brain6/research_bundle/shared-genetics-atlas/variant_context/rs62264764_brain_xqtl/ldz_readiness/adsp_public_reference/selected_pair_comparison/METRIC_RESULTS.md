# LD metric sensitivity within the smaller reference

For the three pairs with the largest previously observed between-reference differences, I recalculated LD in the **same 503 1000 Genomes participants** using both individual genotype-dosage correlations and phased haplotype correlations:

| GRCh38 variant pair | 1000 Genomes dosage r² | 1000 Genomes haplotype r² | ADSP reported r² |
|---|---:|---:|---:|
| chr3:117835650 T>G / 117858571 G>A | 0.706 | 0.678 | 0.434 |
| chr3:117835650 T>G / 117890455 T>G | 0.672 | 0.640 | 0.432 |
| chr3:117520255 T>C / 117698054 C>T | 0.417 | 0.357 | 0.188 |

The two LD calculations differ by 0.028–0.060 within this panel; for these three selected pairs, the ADSP values are lower than both. This does **not** identify the cause of the panel differences or establish an ancestry effect. The exact phase/genotype mode used for the public ADSP file was not independently confirmed; emeraLD supports both phased and unphased inputs ([method paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6298049/)). Sampling, relatedness, ascertainment, phase quality and population composition also differ. The pairs were selected after seeing large discrepancies, so they cannot be generalized to the region or interpreted as a formal statistical test.

All calls used in the haplotype calculation were phased diploid genotypes in the same sample order; dosage r² reproduced the original saved LD matrix for all three pairs. [Protocol](METRIC_PROTOCOL.md), [code](metric_sensitivity.py), [values](metric_sensitivity.tsv), and [receipt](metric_receipt.json) record the check. The remote VCF and index were not whole-file hashed. This does not construct a study-matched LD matrix, establish gene or trait causality, or repair the failed LAVA gates.
