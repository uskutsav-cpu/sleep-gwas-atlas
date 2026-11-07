# Public FOXP2 promoter experiment: context for Brain6

**Finding:** A [2026 peer-reviewed FOXP2 promoter study](https://doi.org/10.1016/j.mcpro.2026.101570) provides a public, experimentally grounded list of possible regulatory proteins, but does **not** test the Brain6 variant rs1476535 or demonstrate a sleep–ADHD mechanism. It is gene-level context, not independent genetic replication or an AlphaGenome benchmark.

## What was measured

The authors used dCas9–APEX2 proximity labeling at the active FOXP2 TSS1 promoter in living HEK293 cells, followed by quantitative mass spectrometry. They identified 373 proteins enriched against controls under their stated thresholds (Storey q < 0.05 and fold change > 1.2). Proximity does not by itself establish direct DNA binding or a regulatory effect; the study compared findings with HEK293 ENCODE ChIP-seq and other evidence.

The authors also examined previously generated [IRF2BP2 knockdown and ChIP-seq data, GEO GSE124636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124636). In HEK cells constitutively expressing glucocorticoid receptor, the reported FOXP2 transcript increase after IRF2BP2 knockdown was approximately 1.6-fold, and dexamethasone treatment decreased expression (n = 4). In A549 lung cancer cells, the knockdown comparison approached but did not meet the nominal 0.05 threshold (p = 0.0558; n = 3). The authors report that the control IRF2BP2 ChIP peaks they examined did not include the targeted promoter. These observations suggest cell-context-dependent regulatory hypotheses, not proof that IRF2BP2 directly binds the targeted promoter or acts through the selected SNP.

The proteomics dataset is publicly indexed as [MSV000098217](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD065063) / PXD065063; the paper lists the ENCODE experiments in its supplementary table. Availability does not establish that these data are suitable training labels for a model.

## What it means for the selected variant

rs1476535 lies in the FOXP2 region and has previously reported sleep–ADHD association, but this study does not identify the allele-specific effects of rs1476535 or establish it regulates FOXP2. The promoter measurements were in a non-neuronal cell line. A prior [FOXP2 enhancer study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5826363/) reported reporter activity for regulatory fragments, but no build-harmonized match of rs1476535 to an assayed fragment has been verified here. No specific brain-cell, stress-hormone, sleep or ADHD pathway is demonstrated by this comparison.

A practical follow-up is to establish genome builds and alleles, intersect the variant and its fine-mapped proxies with the actual experimental intervals, and check molecular QTL colocalization in relevant cells. If a plausible regulatory interval survives those checks, allele-specific reporter and endogenous editing experiments could measure FOXP2 and nearby transcripts alongside accessibility, with independent genetic replication. The protein list can inform a prespecified biological hypothesis; it cannot supply positive or negative variant labels for assessing AlphaGenome's ranking.

For the existing quantitative feasibility check, see [regional reporter readiness](../../prediction_rank_comparison/regional_reporter_readiness/RESULTS.md). For other FOXP2 context and the limitations of the genetic evidence, see the [gene and pathway review](../../pathway_hypotheses/GENE_PATHWAY_REVIEW.md).
