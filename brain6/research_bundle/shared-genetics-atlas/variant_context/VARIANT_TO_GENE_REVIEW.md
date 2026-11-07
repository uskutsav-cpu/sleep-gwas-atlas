# Brain6: selected AlphaGenome variants and gene-context audit

Analysis date: 27 October 2026. We selected, **within each of the six already selected insomnia–ADHD regions**, the variant with the largest absolute predicted adult frontal-cortex DNase change, largest AVI PHRED, largest insomnia ABF PIP and largest ADHD ABF PIP. Ties are resolved by input-row order; overlapping choices were deduplicated, giving **20 distinct variant–region rows**. All selection is exploratory and restricted to the 425-variant union of the two trait-specific credible sets; these are not genome-wide maxima. See `selected_variant_gene_context.tsv` for selection labels, source A1 alleles, both trait-specific estimates and credible-set membership, all coordinates and query links.

We queried Ensembl's GRCh38 REST gene-overlap endpoint for each ±100 kb interval and computed distance to each annotated gene body. All 20 queries succeeded; the endpoint reported release **116** at retrieval. Raw JSON and SHA-256 hashes are preserved. A gene can overlap a variant without being regulated by it, and no gene within ±100 kb is not evidence that none exists further away. Ensembl may label a locus with an identifier rather than a gene symbol. Where multiple genes overlap, their printed order conveys no priority. The endpoint returns LINC03051 with `protein_coding` biotype in this release; its name is not proof of a particular transcript class. Track the stable Ensembl IDs and verify the release before comparing annotation sources.

## Region representatives by absolute local DNase prediction

| Representative | GRCh38 REF>ALT | Local DNase log₂ ratio | Gene-body context in this lookup | Insomnia / ADHD credible-set membership |
|---|---|---:|---|---|
| rs62264764 | chr3:117920728 G>A | −0.2734 | Overlaps LSAMP and LINC03051 | yes / yes |
| rs10835373 | chr11:28622395 T>C | +0.0917 | No overlap; nearest annotated body LINC02758 | no / yes |
| rs79470265 | chr14:98182792 A>C | +0.0853 | Overlaps ENSG00000259097 | no / yes |
| rs2447828 | chr5:104611504 C>T | −0.0590 | Overlaps ENSG00000251574 | yes / no |
| rs2189015 | chr7:114508188 A>G | +0.0499 | Overlaps FOXP2 | yes / no |
| rs324886 | chr5:88600784 C>T | −0.0139 | Overlaps MIR9-2HG | no / yes |

The sign compares the model's ALT and REF local predicted DNase signal; it is not the direction of the GWAS effect. The score is log₂((sum ALT + 1)/(sum REF + 1)) in a 501-bp window. It is not a p-value, an assay measurement, or a causal-variant probability. Both the full 1 Mb input and the 2,001-bp extracted windows were used as documented in the [expanded pilot](../../../alphagenome_pilot/results/jbrowse_expanded/EXPANDED_REPORT.md). Different per-region sample sizes and selected credible sets make maxima unsuitable for cross-region significance comparisons. One-trait credible-set membership alone is not evidence that the same variant causes both traits.

## What the other selection methods highlight

- **Chromosome 3:** rs62264764 is the largest local DNase and AVI candidate in its region, but rs1456193 has the highest insomnia conditional PIP among selected rows (0.0615) and rs73163723 has the highest ADHD conditional PIP (0.0473). All overlap LSAMP and LINC03051, so this query does not resolve the gene. rs73163723 belongs to the ADHD set only. The published Zu et al. study already reports rs62264764 in an insomnia–ADHD PLACO locus with single-signal PP4=0.91. Brain6 chr3 PP4 falls from 0.9207 to 0.5373 under its lower p12 prior. Model-dependent agreement using related discovery inputs is not independent replication. [Zu et al.](https://www.nature.com/articles/s41398-026-04166-4)
- **Chromosome 5 (GRCh37 region 103.45–104.45 Mb):** rs2447828 has the largest local DNase change but is only an insomnia credible-set member. rs2431108 has the greatest AVI and insomnia PIP among selected rows; rs77960 has the highest ADHD PIP, belongs to both sets, and was reported by Zu et al. None of these facts demonstrates that the same base change regulates a gene or causes both phenotypes. Gene bodies returned include ENSG00000251574; Zu et al. use the label RP11-6N13.1 for rs77960. Confirm identifiers and release before treating names as equivalent.

  A later [published reporter and gene-mapping follow-up](crossdisorder_mpra/published_table/GENE_LINK_RESULTS.md) gives an ABC-predicted link from moderately lead-linked rs6874138 to NUDT12. RNA processing and circadian-transcript regulation motivate a specific follow-up, but the mapping is computational and the circadian experiment is in mouse liver; neither establishes a human sleep–ADHD mechanism.
- **Chromosome 7:** FOXP2 overlaps all four selected representatives. rs1476535 (largest insomnia PIP here) was reported in Zu et al.; rs9969232 has the largest ADHD PIP. rs2189015's local DNase maximum does not belong to the ADHD credible set. FOXP2 overlap generates a gene hypothesis, not variant-to-gene or gene-to-disease proof.
- **Chromosome 5 (GRCh37 region 87.19–88.19 Mb):** rs4916723 has ADHD conditional PIP 0.884 but insomnia PIP 0.000243 and is only an ADHD credible-set member. A large single-trait PIP does not establish a shared causal signal. Both selected variants overlap MIR9-2HG; the narrower ±100 kb query may omit more distant proposed genes from the published regional map.
- **Chromosomes 11 and 14:** no selected representative is in the insomnia credible set, even where the table reports an insomnia PIP value. Those values are model-conditional masses for the selected ADHD union, not complete insomnia credible sets. Low predicted functional change or a missing gene annotation is not evidence of no function.

## Pathways and cell types: what is measured or inferred

[Zu et al.](https://www.nature.com/articles/s41398-026-04166-4) report enrichment of mapped shared genes for synapse organization and growth-related functions and exploratory neuronal enrichment. [Xue et al.](https://pubmed.ncbi.nlm.nih.gov/41065713/) studied insomnia with 12 psychiatric disorders, reporting GABAergic-synapse enrichment and cortical GABAergic and glutamatergic cell types in an aggregate analysis. These are independent published analyses at pathway or gene-set level: they do **not** establish that rs62264764 acts through GABA signaling, a specific cortical neuron, or LSAMP. Brain6's 18-region positional sensitivity had no tissue or pathway at BH q<0.05; no cell-resolved enrichment was estimable in the reviewed panel. Frontal cortex in the pilot names the reference track, not a discovered cell type.

A sensible mechanistic follow-up is to compare LSAMP and FOXP2 with other regional genes using measured brain cis-eQTL and sQTL summary statistics, multi-signal trait–QTL colocalization, and enhancer-to-gene evidence; test predeclared cell classes with expression and accessibility reference panels and a matched background; then validate prioritized variant–gene links experimentally. Preserve negative and unavailable tests alongside positives. Frontal cortex is an initial hypothesis, and other brain regions and developmental stages should be considered where adequately powered data exist.

## Reproduction and limitations

`build_context.py` reads the frozen AlphaGenome evidence TSV in the sibling pilot, selects representatives and downloads or reuses the 20 raw JSON responses; the script assumes this original workspace location, so set its `base` if moved. SHA-256 hashes bind each stored response. The source TSV, model track definitions, genome-build liftover and prediction hashes are documented in the sibling pilot. The REST response is live and release-dependent; use the saved JSON for an exact reconstruction. The JBrowse service did not expose the model release. No new RNA gene scoring, QTL colocalization, independent cohort replication, or wet-lab experiment was performed here. Existing AlphaGenome terms and the user's non-commercial-research authorization apply to its outputs.

## New exact-variant measured QTL follow-up

See [GTEx v8 QTL review](gtex_v8/GTEX_QTL_REVIEW.md) for static significant-pair queries of all 20 representatives and dynamic exact rs62264764 tests of LSAMP and RP11-384F7.2 in cortex, frontal cortex and adrenal gland. This evidence adds tissue-specific observations but does not change any causal or Brain6 validation classification.

A [decision-ready rs62264764 follow-up plan](RS62264764_FOLLOWUP_PLAN.md) separates competing LINC03051, LSAMP and alternative-target hypotheses and specifies the evidence needed to distinguish them.
