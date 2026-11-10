# Reference-panel sensitivity of the FI–FFS summary contrast

The [15-marker All of Us insomnia contrast](stage6_aou_eur_fi_ffs_paired_contrast_review.md) used an external 503-person European reference for LD. The [subsequent omission analysis](stage6_aou_paired_contrast_locus_sensitivity.md) found that leaving out rs113851554 reduced the projected contrast. We resampled the 503 reference individuals 5,000 times, recalculating within-chromosome LD, while holding the target GWAS statistics, allele frequencies and score weights fixed. The same resampled individuals were used for both scenarios.

| Score scenario | Original reference contrast Z | Reference-bootstrap 2.5th–97.5th percentiles |
|---|---:|---:|
| All 15 markers | 2.82 | 2.69–2.97 |
| Omit rs113851554 | 1.39 | 1.28–1.51 |

Within this limited resampling exercise, the drop after omission persists. The percentiles **are not confidence intervals for the score effect, association, difference between FI and FFS, replication or prediction**: the target GWAS, selection process, sample overlap, score weights and allele conventions were held fixed. The reference may contain relatives or population structure, and bootstrapping it cannot correct systematic differences between the reference and All of Us. It cannot establish whether the modified score performs better, whether the variant is causal, or whether the FFS weight direction is correct. The marker omission remains post hoc and may leave linked variants in the model.

The script reproduces the original reference correlation matrix and previous full and omitted Z values before resampling. See [results](stage6_aou_paired_reference_bootstrap.csv), [hash and seed receipt](stage6_aou_paired_reference_bootstrap_receipt.json), and `work/aou_complete/eur_paired_reference_bootstrap.py`; run with `/opt/miniconda3/bin/python3` from the project root. No individual target-cohort data were used.
