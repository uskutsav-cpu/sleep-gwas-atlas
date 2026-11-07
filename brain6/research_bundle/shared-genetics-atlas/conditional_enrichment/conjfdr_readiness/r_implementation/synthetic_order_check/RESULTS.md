# Synthetic chromosome-order check for cfdr.pleio

The inspected R package (`cfdr.pleio` 0.0.0.9100, commit `76d5085e6d3f3ca9576d5d7564d2acf11bcfd021`) was installed with its missing dependencies on the SSD work image. A six-variant, two-chromosome synthetic test of the actual `init_data` method shows that the stored trait/variant rows can differ from the order of the chromosome-specific reference matrices:

| Order | Variant IDs |
|---|---|
| Trait and variant table returned by `init_data` | rs5, rs1, rs4, rs3, rs6, rs2 |
| Expected genomic order | rs2, rs6, rs4, rs3, rs1, rs5 |
| Chromosome reference matrices concatenated in chromosome order | rs4, rs6, rs2, rs5, rs1, rs3 |

In the inspected source, `frank` returns ranks, which are then used directly as row indices; ranks need not be the sorting permutation. Subsequently, pruning fills contiguous chromosome blocks using the order of the chromosome-specific matrices. A mismatch can therefore attach pruning decisions to the wrong trait rows. This synthetic example proves the mismatch is possible in this version; it does not quantify its impact on any real dataset. No real study data were run through this package.

The [reproducer](check.R) creates a new reference and output directory, initializes the unmodified installed package, and asserts that the trait rows match the saved variant table but differ from the chromosome matrix order. The [run log](run.log) records the result. The example is a software diagnostic only: it does not test full conjunctional FDR calibration, sample-overlap correction, the real reference or memory use. The issue must be corrected and tested with genomic-order, allele/variant alignment, pruning-index and full-pipeline checks before this implementation can be considered for a study run. The distinct, unimplemented sample-overlap correction remains unresolved. See the [readiness assessment](../RESULTS.md).

A [local correction and 22-chromosome synthetic pruning test](patch_check/RESULTS.md) address this specific ordering defect in a separate copy. Full-method validation and the independent overlap problem remain open.
