# Musculoskeletal factor follow-up for the insomnia–frailty candidates

## What the public study adds

[Gong et al. (2026)](https://www.nature.com/articles/s41467-026-72164-7) used genomic structural equation modeling to construct three latent musculoskeletal factors. Their degenerative factor is a model of shared genetic liability across selected musculoskeletal phenotypes; it is not the clinical frailty index, Fried frailty score, HFRS, or a new individual-level predictive outcome. Its input studies include FinnGen and UK Biobank sources already used, or potentially overlapping sources, in this project. This is therefore **related-phenotype triangulation, not independent replication** of the sleep–frailty findings.

I streamed all 4,367,175 records from the public [MSKdg Zenodo file](https://zenodo.org/records/19012445) and matched the project's existing 15 candidates by GRCh37 rsID, position and allele pair. Twelve matched; three did not have an rsID row and are unobserved in this file, not failed associations. Two matched variants had p<5e−8 in the published factor analysis:

| Candidate, counted allele | Factor beta (SE) | P | Relationship to project evidence |
|---|---:|---:|---|
| rs12536335 G, FOXP2 region | −0.01557 (0.00135) | 7.87e−31 | G also has negative signed associations with MVP insomnia and FinnGen HFRS in the project's aligned results. |
| rs10865954 T, chromosome 3 | +0.01158 (0.00143) | 5.56e−16 | T has a positive insomnia weight; its HFRS estimate was positive but imprecise and did not reach nominal p<0.05. |

These comparisons are descriptive and postselected. Factor coefficients, insomnia odds ratios and HFRS coefficients have different units and cannot be added or used interchangeably. Consistent allele directions are not a colocalization test, do not prove mediation and do not validate the score. The three missing candidate rsIDs are rs8067974, rs1032524 and rs1942262; reasons for absence are not established. All 15 lookups and source alleles are in `stage6_gong_mskdg_15_candidate_lookup.csv`. The factor GWAS uses a model-derived effective sample size (1,131,013 in these rows), which is not a count of unique participants.

## FOXP2 region and cell evidence

The study's Supplementary Data 9 reports a HyPrColoc result for seven musculoskeletal phenotypes at its FOXP2-region locus, with posterior 0.9253 and candidate rs2894699 accounting for 0.328 of the posterior. That SNP also appears in the project's previously obtained MVP and HFRS credible sets. This coincidence is a **hypothesis for follow-up**, not a three-way sleep–frailty–musculoskeletal colocalization: neither insomnia nor frailty was included in the published seven-phenotype test, the disease datasets overlap or may overlap, and the posterior depends on model assumptions and priors. The published candidate probability does not prove the variant is causal, and locus assignment does not prove FOXP2 is the effector gene.

Supplementary Data 14 lists 13 selected RBM6–degenerative-factor expression colocalizations under immune-cell labels, with PP.H4 from 0.905 to 0.928. These are correlated selected rows, not 13 independent replications, and concern the musculoskeletal factor and expression—not insomnia or frailty. Their coloc model assumes at most one causal variant per trait in each region. They are not evidence that rs12536335 in the FOXP2 region acts through RBM6, which is a different genomic locus.

**Cell-data provenance needs clarification before using these as single-cell results.** Gong et al. describe an 89-donor CD4 T-cell single-cell eQTL source, but their selected table includes labels such as B cells, CD8 T cells, monocytes and M2. [DICE](https://dice-database.org/) contains multiple resources: its earlier [15-population DICE study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289654/) used purified immune-cell populations from 91 donors, while a later [single-cell DICE study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9035271/) profiled activated CD4 T cells from 89 donors. The paper and table do not by themselves resolve which eQTL file produced each listed row; we have not independently reconstructed them. These results should be described as **reported immune-cell-context evidence with unresolved source identity**, not as verified single-cell or cell-specific mediation. FOXP2 and SMAD5 do not appear among the selected rows in Supplementary Data 14; that is not evidence they lack an eQTL or a colocalization.

## Reproducibility and limits

`stage6_gong_mskdg_candidate_scan.py` streams the complete compressed source through HTTP range reads, scanning the entire ZIP member. Python's ZIP reader verified the member CRC; the complete ZIP MD5 was not independently verified. The full source was not stored locally. `stage6_gong_mskdg_15_candidate_receipt.json` records the source, row counts and matching. The project uses string-level alleles and positions; multiallelic representation and prior source harmonization remain limitations.

`stage6_gong2026_source_table_extract.py` extracts the paper's relevant source rows and records the downloaded supplement SHA-256. Tables 6, 7, 9 and 14 are selected reports; the reported RBM6 table is `stage6_gong2026_supplement_table14_selected.csv`. In particular, the latter contains significant results only and cannot establish negatives for absent genes. The HFRS and insomnia source lookups, their constraints and the current prediction-validation gap are documented in `stage5_validation_conclusions.md`.
