# Nearby-gene annotations for the 15 exploratory candidates

I queried the Ensembl GRCh37 human gene-overlap endpoint for each candidate's position ±100,000 base pairs. The table `candidate_nearest_gene_annotations.csv` records the nearest annotated **gene body**: distance is zero when the position lies inside a gene's genomic span and otherwise the shortest base-pair distance to its boundary. All equally near gene bodies are retained, so overlapping genes can appear together. All 15 queries returned at least one annotation. Ensembl identifiers, gene biotypes, query URLs and response hashes are included for traceability. An annotation's gene-body span includes introns, and a zero distance is not evidence of a coding change.

| Candidate | Nearest annotated gene body or bodies | Distance (bp) |
|---|---|---:|
| rs113851554 | MEIS1 | 0 |
| rs8067974 | AC002117.1 (antisense) | 0 |
| rs1909650 | PCDH8P1 (pseudogene) | 15,348 |
| rs11152363 | TCF4 | 0 |
| rs1032524 | IPO9-AS1 (antisense); NAV1 | 0 |
| rs323509 | RP11-6N13.1 (lincRNA) | 0 |
| rs78206187 | AC009474.1 (pseudogene) | 25,582 |
| rs1942262 | TCF4 | 16,245 |
| rs4886140 | RPP40P2 (pseudogene) | 3,295 |
| rs12536335 | FOXP2 | 0 |
| rs35881094 | LINC01122 (lincRNA) | 0 |
| rs17550821 | LINC01122 (lincRNA) | 0 |
| rs7734045 | AC009014.3 (lincRNA) | 9,519 |
| rs11790060 | RP11-165J3.6 (antisense) | 3,529 |
| rs10865954 | KLHDC8B | 0 |

These are **location labels only**. They do not establish that a gene changes expression, causes insomnia or frailty, participates in a specific cell type, or carries a causal shared-risk variant. Names and gene boundaries depend on the Ensembl GRCh37 annotation and can differ across releases or approaches. This method includes all annotated gene types and calculates distance to the gene body—not a promoter or transcription start site. For example, rs10865954 is annotated here within KLHDC8B's gene body, whereas Song et al. highlight MIR191 in the broader region using their own mapping; neither label is established as the causal gene by proximity alone. [Song et al.](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full)

The annotations do not convert the candidate set into a valid genetic risk score. Candidates were selected from prior cross-trait findings, are not fully LD-pruned, and nine rsIDs have more than two alleles in the reference. Source-specific allele representations and QC remain important. In particular, rs12536335 overlaps FOXP2's genomic span, but its FI effect had opposite signs in the UKB and TwinGene cohort estimates after the project's direction adapter; the difference in those estimates does not prove a gene mechanism. See `candidate_twingene_sensitivity.md` and `prototype_variant_five_trait_evidence.csv`.

## Follow-up sequence

1. Verify the exact study alleles, genome build, phenotype, sample overlap and cohort-specific effects for each region; retain discordant and imprecise findings rather than choosing only favorable estimates.
2. With dense, QC-passing GWAS and appropriate reference LD, carry out prespecified locus-level shared-signal testing and fine-mapping. The existing project FI×sleep LAVA results failed their prespecified QC gates and cannot be promoted as shared-locus findings.
3. For regions with adequate evidence, compare credible variants to independent gene-mapping evidence—such as colocalized molecular QTLs and regulatory chromatin contacts—then use tissue and single-cell annotations as hypotheses to test. Neither a nearby gene nor expression in a cell proves mediation.
4. Develop any sparse score with frozen, LD-aware selection and training, then independently validate for a specified outcome and time horizon. The existing 15-row prototype was computed for a public reference panel only as a technical check; no outcome was tested and it is not a clinical prediction.
