# Long-sleep replacement source triage — 2026-09-27

This is a source-only feasibility comparison for a separately versioned LAVA
rescue. It does not admit a new trait, run association analysis, change the
canonical seven-trait family, or alter any frozen output.

The locked [Dashti 2019 source card](../../manifests/gwas_external_availability.tsv)
records 34,184 long-sleep cases (self-reported at least 9 hours) and 305,742
7–8-hour controls, all European ancestry. Its sample-size-equivalent balanced
case-control N, `4 * cases * controls / (cases + controls)`, is 122,985.41.
This arithmetic comparator is not a variant-specific N or a prediction of
local heritability.

| Public candidate | Long-sleep cases / controls | Definition | Sample-size-equivalent balanced N | Source-only decision |
|---|---:|---|---:|---|
| [SleepChart 2026](https://www.nature.com/articles/s41586-026-10524-5) | 25,049 / 300,420 | More than 8 h versus 6–8 h, UK Biobank European ancestry | 92,484.64 | Public summary statistics exist, but the phenotype differs and the nominal effective N is 25% lower than the locked source. No evidence yet that it repairs local-h² failures. Do not treat it as a stronger equivalent replacement. |
| [UKB/MVP 2023 meta-analysis](https://www.nature.com/articles/s41467-023-41249-y) | 15,962 long-sleep cases in its EUR meta-analysis | At least 10 h versus normal duration, UKB plus MVP | Not calculated; exact comparator control count not independently verified here | The phenotype is a more extreme tail and the case count is below the locked 34,184. The published result does not establish a stronger equivalent replacement. |

The SleepChart paper's [Methods and data availability](https://www.nature.com/articles/s41586-026-10524-5)
identify REGENIE GWAS, European-ancestry UK Biobank data, and public full
summary statistics at the SleepChart portal/Synapse. The UKB/MVP paper's
[data availability statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC10539313/)
points to dbGaP phs001672.v1.p1 and the Gelernter Lab. Neither source has
been downloaded, hashed, harmonized, or run through LAVA here. Sample count
alone cannot rule out better reference overlap, but it does not justify a
full-family rescue run; a prospective, outcome-blinded source and pilot
protocol would be needed first.

The dominant barrier remains 1,291 canonical `NOT_RUN` cells for long sleep
and 3,720 across all seven traits, against the unchanged maximum of 873.
The MDD2025 trait-only pilot cannot by itself clear the family bound. Its
existing outputs must be audited after `/Volumes/Extreme SSD` remounts.
