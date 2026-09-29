# Sample-overlap assessment

This is a cohort-level source audit, not an estimate of the number of shared participants. Exact pairwise intersection counts are unavailable in the current source metadata and remain `UNKNOWN` in the matrix. Do not describe same-UK Biobank GWAS results as independent replication.

The registered primary Frailty Index source combines UK Biobank and Swedish TwinGene. The Fried Frailty Score source GCST90295968 is identified in the PGS Catalog as 386,565 European-ancestry participants from UK Biobank. Eleven of the twelve locked sleep trait sources are registered as UK Biobank-based; their participant-level overlap with the frailty estimates is expected at the cohort level, but exact intersection sizes remain unknown. The twelfth, sleep apnea, is FinnGen R9, so no common cohort is identified with the FI/Fried source descriptions; this is not proof of zero overlap.

For global LDSC, retain cross-trait intercept and sample-overlap diagnostics and report their limits. Treat UKB-based Fried results as construct comparison, not independent replication. Do not meta-analyze or call replication based only on different phenotype labels.

Evidence sources:

- Local source metadata: `config/traits.tsv` and `frailty_paper/config/target_traits.tsv`.
- Fried GWAS sample description: https://www.pgscatalog.org/score/PGS005229/.
- FI GWAS source paper and downloadable summary statistics: https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/.


## Phase 8 cohort evidence ledger

`frailty_analysis/manifests/cohort_overlap.tsv` records source-level cohort declarations for the frozen 12-trait sleep panel and the acquired/planned frailty and aging endpoints. Validate panel membership and source IDs with `make -C frailty_paper validate-overlap`. `YES` and `NO` are used only when supported by the local registered source descriptions; unresolved consortium membership is `UNKNOWN`. Rows preserve source acquisition state and evidence paths. The `exact_participant_intersection` field remains `UNKNOWN` throughout: cohort membership alone cannot quantify shared participants. In particular, FinnGen sleep apnea and UKB/TwinGene frailty have no common cohort identified in the registered descriptions, but zero participant overlap is not established.

The FI/Fried-by-sleep matrix above remains the pairwise primary-family summary. Same-UKB comparisons are not independent replication. No association has yet been assigned an exact or comparable independent replication status because eligible headline estimates and a validated external replication dataset are not yet available.

## Latent-factor and pneumonia cohort evidence

The Nature Genetics 2025 Supplementary Table 1 identifies UK Biobank in the source cohort of each of the 30 deficit GWAS retained for the latent frailty model, including an insomnia indicator GWAS from UK Biobank. This establishes UKB as a contributing cohort for the latent-factor inputs and shows cohort-level overlap risk with UKB-based sleep GWAS. It does not provide exact shared-participant counts; those remain `UNKNOWN` in `frailty_analysis/manifests/cohort_overlap.tsv`. Multi-cohort sources do not establish the presence or absence of FinnGen, 23andMe or CHARGE at the participant level beyond what each source explicitly documents.

The same article’s Supplementary Table 4 states that the pneumonia context GWAS is a fixed-effect meta-analysis of Pan-UK Biobank and FinnGen release 10. The ledger records UKB=YES and FinnGen=YES, while retaining exact participant intersection as `UNKNOWN`. Pneumonia is not treated as a frailty factor and its non-finite beta/SE rows are not filtered. Publisher workbook URL and checksum evidence are recorded in `frailty_paper/review/targeted_resource_audit.md`.
