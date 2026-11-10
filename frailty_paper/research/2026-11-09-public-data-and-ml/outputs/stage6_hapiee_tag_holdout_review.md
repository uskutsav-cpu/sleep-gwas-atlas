# Reference-tag sensitivity: non-Finnish EUR selection, Finnish holdout

The earlier GSA v3 position-based tag search selected and evaluated tags in the same 503-person 1000 Genomes European reference, potentially overstating how well selected tags generalize. As a sensitivity analysis, `stage6_hapiee_tag_holdout.py` selected the best tag within 100 kb using 404 non-Finnish European individuals and evaluated the same tag in 99 Finnish individuals. All 33 unique candidates had an evaluable reference tag and holdout r².

| Candidate set | Number | Training r² ≥ 0.8 | Holdout r² ≥ 0.8 |
|---|---:|---:|---:|
| 15-variant prototype | 15 | 7 | 7 |
| 22-variant prototype | 22 | 6 | 7 |

These are marginal counts, not necessarily the same variants. For example, rs10865954's chosen tag had training r² 0.905 and holdout r² 0.778. Full per-variant results appear in `stage6_hapiee_tag_holdout.csv`. Thresholds are descriptive; no claim of reliable imputation follows. The holdout is small and Finnish LD can differ from Czech LD. The broader 1000 Genomes European data had already been examined in the earlier exercise, so this is a secondary sensitivity analysis rather than a wholly independent study.

The GSA v3 manufacturer manifest represents candidate probe *positions*, without verifying assay alleles, strand, actual cohort platform, observed calls, sample-level quality, or target imputation INFO. In particular, the HELIAD dataset uses a different GSA version. Tag r² is neither a calibrated estimate of imputation accuracy nor evidence of predictive score validation. Obtain actual cohort genotype dosages, variant-level QC and phenotypes for those purposes.

The reference data are from the [1000 Genomes Project](https://www.internationalgenome.org/data); the script uses the public UCSC hg19 mirror recorded in the adjacent receipt. Reproduce with the local manifest and panel listed in the script and network access to UCSC.
