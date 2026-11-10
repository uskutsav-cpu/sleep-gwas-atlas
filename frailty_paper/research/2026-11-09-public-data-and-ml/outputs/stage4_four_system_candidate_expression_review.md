# Candidate-gene expression across four systems

The [expression table](stage4_fuma_candidate_gene_expression.csv) maps stable candidate-gene IDs to seven publicly available processed [FUMA scRNA matrices](https://github.com/vufuma/FUMA_scRNA_data/tree/dd526163ea80af1a80a6cdc80db167144500694b), covering brain, immune cells, muscle and vasculature. These are descriptive source-provided cell-type expression values, **not enrichment or differential-expression tests**. They allow the evidence network to record which selected candidate genes were measured in which reference cell types.

| Processed matrix | Biological system | Source species | Cell-type columns | Candidate IDs present |
|---|---|---|---:|---:|
| GSE67835 human cortex without fetal samples | Brain | Human | 7 | 143 |
| GSE87544 hypothalamus | Brain | Mouse; mapped by source to human IDs | 46 | 130 |
| PBMC 10x 68k | Immune | Human | 11 | 164 |
| GSE98816 brain vasculature | Vascular | Mouse; mapped by source to human IDs | 15 | 136 |
| GSE168408 adult human prefrontal cortex level 2 | Brain | Human | 15 | 154 |
| Mouse Cell Atlas muscle | Muscle | Mouse; mapped by source to human IDs | 17 | 136 |
| Tabula Muris FACS limb muscle | Muscle | Mouse; mapped by source to human IDs | 6 | 130 |

The source matrices were selected from available references for coverage, not preregistered. The source provides processed expression values, and values across datasets and species are not necessarily directly comparable. Mouse–human ortholog mapping is inherited from the source; exact homolog mapping for each gene has not been independently reviewed. Cell-type labels are the source labels and should not be equated automatically with clinical populations. Dataset-specific donor numbers, processing, sparsity and batch effects require further checking. These data do not represent every hypothalamic, immune, muscle, or vascular subtype, and the mouse datasets are not human evidence.

The 182 stable IDs, including provisional IDs matched to proxy-predicted target labels, came from the current provisional gene table; candidates without stable IDs there are not included. [Coverage by dataset](stage4_fuma_candidate_expression_coverage.json) records which IDs are absent. The table includes 16,227 measured gene–cell-type entries, including zeros; an absent gene or a zero is not proof of biological absence. The [network](stage4_variant_gene_phenotype_network_scaffold.csv) includes the 10,445 entries with positive source values as gene-to-cell expression edges, identified by source dataset and species. Rare or low-level positive values are retained rather than treated as significant. The matrix [receipt](stage4_fuma_matrix_receipt.json) records the archive members and local hashes.

A path from a variant to a gene to a cell type to a phenotype joins different evidence types. It does **not** show that the variant acts through that gene in that cell, or that the cell mediates the phenotype. GWAS/QTL colocalization and properly controlled cell-type enrichment remain necessary; separate [MAGMA analyses](stage4_standard_fi_ffs_cell_type_comparison.md) are available. The separate [human muscle donor summaries](stage4_human_musc_candidate_expression_review.md) retain donor information that these aggregate source matrices do not supply.

The GSE168408 adult prefrontal-cortex matrix was added after the original six datasets. Its [secondary MAGMA analysis](stage4_adult_pfc_neuron_subtype_review.md) uses a separate 60-test correction family.
