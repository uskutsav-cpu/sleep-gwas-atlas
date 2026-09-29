# Latent-factor cohort/overlap source recheck

Date: 2026-09-23.

## Primary-source findings

The primary Foote et al. paper describes selecting publicly available GWAS for 30 frailty deficits from an initial 52 traits. Eligibility required at least 10,000 European-ancestry participants; consortium GWAS were prioritized, with Pan-UK Biobank summary statistics used when consortium data were unavailable. The paper states that Genomic SEM can model input GWAS with varying and unknown sample overlap, including mutually exclusive samples. The final seven latent-factor GWAS were created from this multivariate model.

Primary source: [Foote et al., Nature Genetics (2025)](https://www.nature.com/articles/s41588-025-02269-0), Methods “Phenotype selection” and Results “Multivariate genetic architecture of frailty”. The article gives aggregate inclusion criteria and method capability, not a complete accession-by-accession participant/cohort overlap matrix in the accessible main text.

## Relevance to this project

This confirms that the seven latent-factor GWAS are European-ancestry and that unknown overlap among their 30 input GWAS was intentionally handled by the authors’ Genomic SEM procedure. It does **not** establish whether particular source GWAS include UK Biobank, TwinGene, FinnGen, or other cohorts shared with any of the 12 sleep GWAS used here; nor does it quantify pairwise intersections with those sleep GWAS. Therefore the latent-factor sleep rg family remains secondary/sensitivity evidence and is not independent replication. The article’s tolerance of unknown overlap inside its own model is not evidence of zero overlap in downstream cross-study comparisons.

## Decision / follow-up

Retain exact cohort overlap as UNKNOWN. Use the primary article’s Supplementary Table 1 and source GWAS metadata to map each contributing frailty-deficit GWAS where possible, then compare against each sleep source. If any component source or intersection remains unidentified, preserve UNKNOWN rather than infer zero. Do not change estimates, thresholds, or significance classifications on this basis.
