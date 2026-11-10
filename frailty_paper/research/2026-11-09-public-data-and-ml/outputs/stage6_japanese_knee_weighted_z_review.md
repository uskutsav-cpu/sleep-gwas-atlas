# Japanese knee-strength GWAS: exploratory partial-score screen

Using the checksum-verified Japanese GWAS, I combined the available candidate associations for each existing weight set. The source measures **higher** knee strength; I reversed the signed SNP statistics so positive Z in this table would correspond to the score's weighted direction aligning with **lower** strength. All four observed Z values were negative and statistically uncertain:

| Prototype | Markers used | Lower-strength weighted Z | Nominal two-sided p | Full vs partial genotype-score correlation in JPT reference |
|---|---:|---:|---:|---:|
| FI weights | 12/15 | −0.52 | 0.60 | 0.956 |
| FFS weights, inferred direction | 12/15 | −0.61 | 0.54 | 0.998 |
| Insomnia weights | 12/15 | −0.69 | 0.49 | 0.967 |
| Equal-weight count | 14/22 | −0.38 | 0.70 | 0.970 |

These postselected, approximate partial-score statistics do **not** establish an association in either direction, refute a small effect, validate a score, or rank FI versus FFS predictive accuracy. The four tests overlap and were not corrected for the broader analyses already conducted. The small GWAS (n=3,452), phenotype differences, ancestry transfer, missing markers and 104-person reference LD make uncertainty substantial. The data are summary associations for knee strength rather than FI, FFS or a joint sleep–frailty outcome. Full participant/relatedness non-overlap with all project sources is not independently established.

For each matched SNP, the source effect is aligned to the prototype's counted allele assuming the unharmonized source effect-allele coding is correct, then multiplied by −1 to represent lower strength. The paper's lead variant matches the source's signed effect, providing a control but not proving each project variant's strand or multiallelic representation. The FFS source weight directions remain inferred. With `z` the oriented marginal beta/SE and `a=weight × sqrt(2p(1−p))` using the source frequency, the statistic is `aᵀz / sqrt(aᵀRa)`; `R` is JPT reference dosage correlation within chromosomes, with cross-chromosome correlations zeroed. This approximates target Z covariance; it is not a regression of the actual score in target individuals. Reference genotype score correlations compare full and partial sums only in 104 1000 Genomes JPT samples, not outcomes or the target cohort. The reference VCF was accessed by indexed region without full-file checksum verification.

A subsequent [allele and frequency audit](stage6_japanese_knee_allele_frequency_audit.md) found reference allele-order matches for all 23 distinct included variants and no palindromic pairs, with the source and small JPT-panel effect frequencies differing by at most 0.0634. This strengthens the consistency check but does not prove study-specific allele semantics or prediction.

The [receipt](stage6_japanese_knee_weighted_z.json), [component and leave-one-out table](stage6_japanese_knee_weighted_z_components.csv), and `python3 outputs/stage6_japanese_knee_weighted_z.py` document the calculations. The [candidate screen](stage6_japanese_knee_candidate_review.md) documents source QC and missing variants. Proper individual score validation remains outstanding.
