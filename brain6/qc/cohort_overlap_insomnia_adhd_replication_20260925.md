# Insomnia–ADHD and FinnGen replication cohort-overlap audit

Audit date: 2026-09-25. Scope: public study descriptions for the locked insomnia and ADHD discovery inputs and the archived FinnGen R13 ADHD directional replication. This is a cohort-description audit; it does not use individual-level records or estimate participant overlap.

## Evidence reviewed

- The locked insomnia input is the UK Biobank-only release of Jansen et al. 2019 (N=386,533; 109,402 cases and 277,131 controls), as recorded in `manifests/gwas_master.tsv`.
- The locked Demontis et al. ADHD discovery description reports iPSYCH, deCODE, and ten European-ancestry PGC cohorts. The article identifies iPSYCH diagnoses through Danish psychiatric registers and describes the deCODE cohort separately. The ten PGC cohort identities were not recovered from the main-text source used in this audit.
- FinnGen describes its cohort as Finnish biobank samples linked to longitudinal Finnish health-register data. Its R13 F5 ADHD endpoint is based on register diagnoses (including ICD-10 F90.0 and ICD-9 3140 in the endpoint documentation), a distinct phenotype definition and ascertainment route from the ADHD discovery meta-analysis.
- FinnGen R13's archived ADHD replication estimate uses the same UK Biobank insomnia exposure as discovery. The outcome estimate is therefore directional/partial replication, not fully independent two-trait replication.

## Assessment

The reviewed sources do not establish a shared named cohort between the UK Biobank insomnia input, the Demontis ADHD discovery, and the FinnGen R13 outcome. This is limited negative evidence: the ten PGC cohorts were not enumerated in the main-text material reviewed, the study sources do not provide person-level linkage, and exact participant/control overlap was not checked. The `insomnia__adhd` overlap classification therefore remains `UNKNOWN`. FinnGen is recorded as a distinct national biobank/register outcome source, while retaining the Finnish-founder and phenotype-ascertainment caveats. No result estimate or replication classification is changed.

## Sources

- Demontis et al. 2023, primary article: <https://pmc.ncbi.nlm.nih.gov/articles/10914347/> (DOI: 10.1038/s41588-022-01285-8).
- Jansen et al. 2019, primary article: <https://www.nature.com/articles/s41588-018-0333-3> (DOI: 10.1038/s41588-018-0333-3).
- FinnGen cohort and data description: <https://www.finngen.fi/en/cohort_and_data>.
- FinnGen Risteys R11 F5_ADHD endpoint definition: <https://r11.risteys.finngen.fi/endpoints/F5_ADHD>.
- FinnGen DF13 access/results information: <https://www.finngen.fi/en/access_results>.

This note records the evidence boundary only. It does not certify no overlap, upgrade the replication level, or support a causal interpretation.
