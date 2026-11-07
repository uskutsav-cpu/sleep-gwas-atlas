# Chromosome 5 brain xQTL follow-up

The [NIAGADS brain xQTL atlas](https://xqtl.niagads.org/variant/rs4916723) returned **16 HMT-filtered molecular associations** for rs4916723: six bulk expression, four single-nucleus expression, two splicing and four methylation records. This marker belongs to the **ADHD-prioritized chr5:87.19–88.19 Mb region**, a different region from the stronger shared insomnia–ADHD colocalization at chr5:103.45–104.48 Mb. In the existing Brain6 single-signal analysis its conditional PIP is 0.884 for ADHD but 0.000243 for insomnia; the region's shared-signal PP4 is 0.138 at the stated default prior. Neither PIP establishes causality, and these data do not make this a shared insomnia locus.

Selected measured expression associations at rs4916723 (chr5:88558577 A>C, GRCh38; effects for ALT C):

| Molecular target and context | Study | Beta (SE) | Nominal p | Per-target FDR |
|---|---|---:|---:|---:|
| TMEM161B-AS1, oligodendrocyte progenitor cells | ROSMAP_CUIMC1_2_MIT | −0.1779 (0.02655) | 4.45×10⁻¹¹ | 1.41×10⁻⁹ |
| LINC02060, astrocytes | ROSMAP_CUIMC1_2_MIT | +0.1812 (0.02473) | 6.99×10⁻¹³ | 1.89×10⁻¹¹ |
| TMEM161B-AS1, bulk dorsolateral prefrontal cortex | ROSMAP | −0.1051 (0.01512) | 8.29×10⁻¹² | 2.54×10⁻¹⁰ |
| AC091826.2, bulk dorsolateral prefrontal cortex | ROSMAP | −0.2859 (0.03866) | 3.91×10⁻¹³ | 1.46×10⁻¹¹ |

The portal also reports a LINC02060 single-nucleus association in excitatory neurons, additional TMEM161B-AS1 bulk and splicing associations, and CpG methylation associations. The target of an mQTL is the indicated CpG, not expression of the gene used to annotate it. Different cell-type coefficients are not a formal interaction test; reported betas are not fold changes. Some datasets share ROSMAP donors and results should not be counted as independent replications. These associations arise in aging-brain data and are not evidence that the same effect occurs during ADHD development.

The portal's [methods](https://xqtl.niagads.org/about) define GRCh38 coordinates, ALT-allele effects, target-level BH FDR and cross-target hierarchical multiple testing (HMT). The download is HMT-filtered, so the listed FDRs do not correct for Brain6's post-hoc choice of region or all follow-up searches. These are single-variant molecular associations, not GWAS–QTL colocalizations; linked variants or multiple signals could explain them. They do not assign an ADHD effector gene, identify a mechanism, measure chromatin accessibility, validate an AlphaGenome prediction, or demonstrate insomnia–ADHD sharing. Earlier GTEx brain TMEM161B-AS1 QTL records and this follow-up supply supporting molecular context, not independent proof of causality. Gene names, proximity to MEF2C and cell annotations do not establish a pathway.

For the distinct strong shared-signal region, the portal's filtered query for [rs2431108](https://xqtl.niagads.org/variant/rs2431108) and [rs77960](https://xqtl.niagads.org/variant/rs77960) returned no records on 27 October 2026. That is **not** evidence of no brain association, no expression or no testing. The filtered search cannot replace a dense, tissue-matched QTL analysis or the existing [testis QTL follow-up](../chr5_dense_qtl/RESULTS.md).

The [source table](rs4916723.tsv), empty-response records for [rs2431108](rs2431108_no_records.html) and [rs77960](rs77960_no_records.html), [verification](verify.py), and [hash receipt](verification.json) preserve the lookup. Downloaded from the portal's per-rsID download endpoint; its [download page](https://xqtl.niagads.org/download) states that portal queries return HMT-significant associations. The checks verify row consistency, allele labels and beta/SE arithmetic, not statistical calibration, gene mediation or independent locus replication. A next step is harmonized dense cell-specific eQTL and GWAS summary data, with study-appropriate LD and multi-signal colocalization before perturbation experiments. Failed LAVA QC remains unresolved.

A [dense brain-cell expression–GWAS coloc pilot](dense_coloc/RESULTS.md) favors distinct dominant signals for the selected expression traits and ADHD, while insomnia results are sensitive to the shared-signal prior. Multiple signals, selected support and other model limitations remain.
