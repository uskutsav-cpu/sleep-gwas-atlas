# rs7734045: FI, FFS and insomnia follow-up

rs7734045 is one of the 15 post-selected project candidates. In the reviewed inputs, the C allele has positive point estimates for FI, FFS and all three insomnia sources. The unadjusted p-values are 0.00545 for FI meta-analysis, 0.015 for FFS, 1.59×10⁻⁶ for UKB insomnia, 0.0352 for FinnGen insomnia and 0.0184 for MVP insomnia. This is a useful pattern for prioritizing follow-up, **not** a corrected association or independent validation: it was selected after viewing related cross-trait results, source cohorts may overlap, phenotype scales differ, and the FFS effect-allele interpretation remains inferred.

The two FI component estimates are positive after the project's inferred TwinGene adapter: UKB β=0.00934 (SE 0.00367) and TwinGene β=0.02131 (SE 0.01599). The TwinGene estimate is imprecise (unadjusted p=0.183); a nonsignificant result does not show no effect. Their estimated difference is not nominally significant under an exploratory independent-estimate assumption (p=0.466); this does **not** establish equivalence, particularly with different FI measures and an unconfirmed TwinGene allele adapter. The FI meta-analysis contains both components and should not be counted as another independent result.

## Nearby genes, not proven effectors

At GRCh37 chr5:135,538,341, the Ensembl gene-body inventory within ±100 kb contains:

| Gene | Gene type | Distance to gene body |
|---|---|---:|
| AC009014.3 | Long noncoding RNA | 9,519 bp |
| TRPC7 | Protein coding | 10,658 bp |
| TRPC7-AS1 | Antisense RNA | 11,395 bp |
| SMAD5 | Protein coding | 13,906 bp |
| SMAD5-AS1 | Antisense RNA | 67,762 bp |

These distances come from [Ensembl GRCh37 gene annotations](https://grch37.rest.ensembl.org/overlap/region/human/5:135438341-135638341?feature=gene;content-type=application/json); annotations and transcript boundaries depend on database release. The nearest gene is AC009014.3, but naming it alone would obscure the nearby protein-coding alternatives. The complete identifiers and source hash are in `rs7734045_nearby_gene_inventory.csv`.

The large [Jansen et al. insomnia GWAS](https://research.vu.nl/files/83591711/Genome_wide_analysis_of_insomnia_in_1_331_010_individuals_identifies_new_risk_loci_and_functional_pathways.pdf) highlighted TRPC7 in a gene-based insomnia analysis. This is gene-level context, **not** proof that rs7734045 is the causal variant for TRPC7 or that TRPC7 mediates frailty. That GWAS included UK Biobank and is related to the insomnia discovery evidence in this project, so it should not be presented as an independent confirmation of this candidate. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/57113) describes TRPC7 as a protein-coding cation channel with predicted roles in calcium transport. Those general molecular functions do not demonstrate a sleep–frailty mechanism.

The existing allele-specific VEP query returned no transcript consequence for the A/C contrast at this site. This does not rule out distal regulation, unannotated effects, or a role for a nearby gene. A credible next step is to obtain dense, harmonized FI, FFS and insomnia data and test regional shared signals and gene expression QTL colocalization, including TRPC7, SMAD5 and AC009014.3 where measurements exist. The published significant-result tables screened in `rs8067974_foote_qtl_followup.md` did not list rs7734045 as a top SNP; because those tables were filtered, that is not evidence of no QTL effect. No single-cell assignment, causal gene or validated score weight has been established.
