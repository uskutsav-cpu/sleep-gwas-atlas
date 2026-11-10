# Genome-wide audit of the TwinGene direction adapter

The expanded source check gives strong **within-source evidence** that the TwinGene beta column must be sign-reversed relative to the shared allele label to reproduce the FI meta-analysis. This is an inference from the released files, not a correction confirmed by the study authors. The original cohort README describes `EA` as the effect allele for both cohorts, which conflicts with this reconstruction.

I deterministically sampled rsIDs from the project's filtered FI meta-analysis using an rsID-only hash (approximately 1 in 64) before comparing effects. After matching rsID, GRCh37 chromosome/position and the unordered allele pair to the original UK Biobank/TwinGene cohort archive, and removing eligible duplicate source rsIDs, **97,901 unique matches** had finite UKB and TwinGene betas and positive standard errors. For each, I computed the usual inverse-variance weighted meta beta on the meta file's A1 under two hypotheses: the reported TwinGene beta is on the shared source EA, or its sign is reversed relative to that EA. The UKB effect was aligned to meta A1 in both calculations.

| Check in the 97,901 sampled matched rows | Treat TwinGene beta as reported | Reverse TwinGene beta |
|---|---:|---:|
| Mean absolute difference from meta beta | 0.002373 | 0.0000250 |
| Meta beta reconstructed within 0.0001 | 3,362 | 97,901 |
| Lower absolute reconstruction error | 852 | 97,049 |

The FI meta betas in the harmonized file are rounded to four decimals, so a residual around 0.000025 is consistent with ordinary rounding; this is not exact equality. To avoid ambiguous cases where the two predictions are nearly identical compared with rounding, I also examined the **91,143 rows where the two hypothesized meta betas differ by at least 0.0002**. The reversed convention fit better in all 91,143 and fell within 0.0001 in all of them; the unreversed convention did not fall within 0.0001 for any. This threshold was an explanatory sensitivity, not a formal independent-replication test. The variants are correlated and no inferential p-value is attached to the counts.

The sampled meta file has 6,291,635 rows and is a filtered, harmonized derivative, not every variant in the original GWAS. The original cohort archive has 16,557,800 rows. MHC, ambiguous and other previously filtered variants are not represented; matching and available-coefficient criteria further restrict the check. Counters for allele mismatches and invalid values are row-level diagnostics and can overlap with other records for the same rsID; do not subtract them as mutually exclusive categories. The genome-wide pattern substantiates the previously inferred adapter beyond the selected 15 candidates but **cannot determine why the file disagrees with its README**, authenticate every allele representation, rule out source processing or establish biological effects. It also does not validate the FFS effect convention.

## Implication for the leading candidate

For rs12536335, counted G, the adapted cohort estimates remain UKB FI −0.01889 and TwinGene FI +0.03922. With inverse-variance weighting they reconstruct the meta estimate as −0.015788, consistent with the published four-decimal −0.0158. The apparent sign disagreement is therefore not easily dismissed as an arbitrary direction choice: the alternative sign fails to reproduce the supplied meta estimate under the documented calculation. However, the FI instruments, samples and analyses differ; study-author confirmation of cohort coding and an independent, harmonized replication are still needed. Neither this audit nor a disagreement identifies a causal mechanism or proves an effect of removing insomnia items.

## Provenance and reproducibility

The audit reads the existing harmonized meta file (SHA-256 `241ed294ccc367b22787a0c051f8d377c020e16089e82920a555286f3148e6de`) and the [Atkins original cohort archive](https://figshare.com/articles/dataset/Genome-Wide_Association_Study_of_the_Frailty_Index_-_Atkins_et_al_2021/9204998). Its deterministic hash rule, match counts, thresholds and summary results are stored in `twingene_direction_genomewide_audit.json`. Script: `work/audit_twingene_adapter_genomewide.py`. No source files or effect estimates were changed. The findings remain conditional on the harmonized meta A1/BETA convention and the inverse-variance reconstruction.
