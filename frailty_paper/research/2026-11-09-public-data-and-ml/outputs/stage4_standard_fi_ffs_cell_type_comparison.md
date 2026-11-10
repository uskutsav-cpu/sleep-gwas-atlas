# Standard-model FI and FFS cell-type comparison

Using the same MAGMA v1.10 SNPwise-mean gene model, the same 1000 Genomes European reference and six FUMA expression matrices, we tested 102 cell-group expression properties for each trait. Each property model conditions on the matrix's average expression and MAGMA's gene-level covariates. BH adjustment was performed separately across all 102 groups for each trait; its calculation was independently checked with R `p.adjust`. The group selection and this analysis are exploratory. Full results, species, gene counts and regression estimates are in [the cell-type table](stage4_magma_cell_type_exploratory.csv).

| Human cortex expression group | Mak FI p | Mak FI BH q | Ye FFS p | Ye FFS BH q |
| --- | ---: | ---: | ---: | ---: |
| hybrid (source label) | 7.28 × 10⁻⁷ | 7.43 × 10⁻⁵ | 0.00172 | 0.0744 |
| neurons | 6.49 × 10⁻⁶ | 0.000331 | 0.00147 | 0.0744 |

The FI analysis had two groups below q = 0.05 and FFS had none. These are separate regressions: differing significance thresholds do **not** demonstrate that FI and FFS differ in their cell biology. The source label “hybrid” is preserved rather than interpreted as a precisely identified cell type. Some groups and matrices overlap biologically; these findings do not implicate independent cell types or a particular candidate gene. The human cortex expression measures are source-processed averages, not a donor-aware analysis of frail individuals. Mouse data are cross-species. There is no matched insomnia analysis in this standard-model comparison, and this result does not establish shared insomnia–frailty enrichment or mediation. FI and FFS differ in phenotype construction, samples and input characteristics; power differences are not controlled by this side-by-side display.

For context, a separate approximate Brown-1df FI analysis also identified the two cortex groups and mouse hypothalamus Glu3 at q < 0.05; its results are not the standard-model FI results presented above. See [FI Brown review](stage4_mak_fi_brown_magma_cell_type_review.md), [FFS review](stage4_ye_ffs_magma_cell_type_review.md), and [MAGMA preparation and quality checks](stage4_magma_preparation.md). The matched Brown analyses for additional insomnia and FI sources are still running.
