# Brain molecular QTL follow-up for rs62264764

Checked 27 October 2026 using the [ADSP FunGen xQTL atlas](https://xqtl.niagads.org/variant/rs62264764), release 1.0. This is an association lookup in aging brain tissue, not a test of sleep or ADHD causation. The [portal methods](https://xqtl.niagads.org/about) describe GRCh38 coordinates, ALT-allele effects, Benjamini–Hochberg FDR per molecular target, and hierarchical multiple testing (HMT) across targets. The variant search returns HMT-significant associations, not the complete set of tested results.

For rs62264764 (chr3:117920728 G>A, GRCh38), the portal returned nine records: seven methylation QTLs and two single-nucleus expression QTLs. Both expression records involve ENSG00000239268, labeled AC092691.1 by the atlas and LINC03051 in the current [Ensembl lookup](https://rest.ensembl.org/lookup/id/ENSG00000239268?content-type=application/json); the earlier GTEx label was RP11-384F7.2. This is a distinct transcript from LSAMP. Gene naming alone does not establish its function.

| Study and cell type | Target | ALT A effect (beta, SE) | Nominal p | Per-target FDR |
|---|---|---:|---:|---:|
| ROSMAP_CUIMC1_2_MIT, oligodendrocytes | ENSG00000239268 | +0.2883 (0.0489) | 6.01×10⁻⁹ | 3.59×10⁻⁷ |
| ROSMAP_CUIMC1_2_MIT, excitatory neurons | ENSG00000239268 | −0.1191 (0.0216) | 5.15×10⁻⁸ | 9.94×10⁻⁷ |

These are study-reported expression association coefficients; their units should not be interpreted as fold changes. Opposite signs do not establish a statistically significant difference between cell types without a formal interaction test. The FDR values and the portal's HMT filter do not correct for the post hoc selection of this locus in Brain6. They are not evidence that expression mediates either trait.

The strongest returned methylation association was with CpG cg16634544 at chr3:117801805–117801807 in ROSMAP dorsolateral prefrontal cortex (beta +0.02014, SE 0.00202, FDR 9.87×10⁻²⁰). The methylation target is a CpG; the portal's mapped gene annotation is not an expression measurement.

The [rs62264767 lookup](https://xqtl.niagads.org/variant/rs62264767) returned the same nine target/study records and identical reported association statistics, with its ALT C as effect allele. These are not nine independent confirmations. In a separate [1000 Genomes EUR reference check](../../replication/rs62264764_mvp/RESULTS.md), the two alleles had identical observed haplotypes (r²=1 in 503 participants); this does not establish the LD in the xQTL cohort. Portal fine-mapping returned no results for either variant, which does not imply that fine-mapping evidence is absent from all sources.

This finding adds cell-specific expression and methylation associations; it does **not** measure a chromatin-accessibility QTL or validate the AlphaGenome DNase prediction. It does not provide GWAS–QTL colocalization, a causal variant, an effector gene for either disease, or an independent replication of the sleep–ADHD association. Bulk GTEx frontal-brain results and the single-nucleus results are different assays and samples; the earlier GTEx LSAMP result remains unchanged. Follow-up requires harmonized dense brain-cell QTL and GWAS summary statistics, suitable LD and multi-signal colocalization, then experimental perturbation.

Reproducibility: the source downloads [rs62264764.tsv](rs62264764.tsv), [rs62264767.tsv](rs62264767.tsv), [Ensembl response](ensembl_gene.json), [verification script](verify.py), and [verification output and hashes](verification.json) are preserved here. The source portal association records are evidence of measured molecular associations, not validation of a causal model. See also the portal's [2026 study preprint](https://doi.org/10.64898/2026.05.21.26353713).

A [dense oligodendrocyte expression–GWAS pilot](dense_coloc/RESULTS.md) favors distinct dominant LINC03051 and insomnia/ADHD signals under a single-signal model. It cannot exclude secondary sharing or establish an effector gene.

A subsequent [excitatory-neuron regional pilot](excitatory_coloc/RESULTS.md) also favors distinct dominant LINC03051 and disease signals under the single-signal model. Both cell-type results remain compatible with untested secondary sharing.

A subsequent [LD–Z readiness check](ldz_readiness/RESULTS.md) shows close reference allele frequencies but notable residuals for both cell datasets and the previously analyzed GWAS. Study-matched LD and model validation remain necessary for multi-signal interpretation.

A [same-marker FinnGen insomnia sensitivity](excitatory_coloc/finngen_sensitivity/RESULTS.md) is dominated by the eQTL-only hypothesis under the tested single-signal model. Its prior-sensitive posterior does not resolve causal sharing or identify a gene.

A [study-produced LSAMP fine-mapping and brain-cell GWAS follow-up](lsamp_finemapping/RESULTS.md) identifies two molecular credible sets; the single-signal expression–GWAS comparisons favor distinct dominant signals. These results do not establish an effector gene or exclude secondary sharing.

A [developmental ATAC peak lookup](../rs62264764_developmental_caqtl/RESULTS.md) found a peak overlapping the AlphaGenome summary window but no candidate-variant rows in the selected caQTL files. This is measured chromatin context, not a measured effect of either candidate allele.

A [public forebrain-organoid MPRA supplement lookup](../organoid_mpra_2024/RESULTS.md) found no reporter-coordinate match under either GRCh37 or GRCh38 positional assumptions; it does not test the candidate allele effect or rule out activity.
