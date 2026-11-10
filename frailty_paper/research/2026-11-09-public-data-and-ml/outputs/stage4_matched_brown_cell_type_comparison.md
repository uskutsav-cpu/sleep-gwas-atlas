# Matched-model cell-type comparison: insomnia, FI and FFS

I used MAGMA v1.10 Brown 1-df gene results with the same 1000 Genomes European LD reference and six FUMA processed expression matrices for Jansen insomnia, Atkins FI, Mak FI and Ye FFS. Each trait was tested separately for positive association with 102 cell-group expression properties, conditioning on each matrix's average expression and MAGMA's internal gene covariates. These exploratory matrices span brain, immune, vascular and muscle systems, but do not represent every cell type or constitute individual-level cell experiments. Full results are in [the cell-property table](stage4_magma_cell_type_exploratory.csv).

| Study and phenotype | Human cortex neurons p | BH q across 102 planned groups | Groups with q < 0.05 |
| --- | ---: | ---: | ---: |
| Jansen insomnia | 0.00033556 | 0.0342 | 1 of 102 estimable |
| Atkins FI | 0.0047452 | 0.1328 | 0 of 101 estimable |
| Mak FI | 0.0000023837 | 0.0001216 | 3 of 102 estimable |
| Ye FFS | 0.00097979 | 0.0600 | 0 of 102 estimable |

The neurons group passed each separately adjusted threshold for Jansen insomnia and Mak FI. Mak FI's other two q < 0.05 groups were the human cortex source-labeled **hybrid** group and mouse hypothalamus **Glu3**. Ye FFS had no group below 0.05 under either its Brown or default gene model. Atkins FI had no estimable group below 0.05. These are separate tests of gene-property association, not a joint insomnia–frailty test; overlap does not identify shared variants, causal neurons or a variant→gene→cell mechanism. Different phenotype definitions, GWAS sample sizes and coverage, source cohorts, expression models and power affect the comparison. Some UK Biobank participants may overlap, and an FI measure includes a sleep-related item. Matrix and model choices were exploratory; mouse labels are cross-species and the hybrid label is the source's name, not a validated clinical cell identity.

One Atkins PBMC group, `CD4+_T_Helper2`, was **untestable**: MAGMA reported collinearity with the conditioned covariates. Its p-value and q-value are missing, not zero or a nonsignificant result. For this trait, correction retains all 102 planned hypotheses in the BH denominator, including the failed test, and reports q only for 101 estimable groups. Other four-source analyses have 102 estimable groups. Every trait/model BH calculation, including the missing-test case, was independently checked using R `p.adjust` with the planned family size.

The source gene analyses differ in sample-size information and GWAS input filters; see [MAGMA preparation](stage4_magma_preparation.md) and [the four-source gene comparison](stage4_candidate_brown_four_trait_gene_comparison.md). Jansen's saved analysis log records completion and 19,007 output genes even though its console progress capture ended early; MAGMA successfully read its intermediate gene file for all six cell-property analyses. These checks support computation integrity, not biological causality. The original default-model comparison is [here](stage4_standard_fi_ffs_cell_type_comparison.md), and FFS model sensitivity is [here](stage4_ye_ffs_magma_model_sensitivity.md).

Each q-value corrects its specified trait/model family, not every exploratory family, model choice and locus analysis in the project.
