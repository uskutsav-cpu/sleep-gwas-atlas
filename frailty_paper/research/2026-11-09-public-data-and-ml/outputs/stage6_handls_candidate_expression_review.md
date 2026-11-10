# Candidate genes in a public middle-aged frailty blood-expression dataset

The public [GEO series GSE129534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129534), from the HANDLS cohort, contains PBMC RNA-seq count files for 8 frail and 8 non-frail participants, ages 45–49. The authors describe frailty using the FRAIL scale in their [follow-up paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9361278/). This is not a participant-level genotype dataset or a study of insomnia, and does not validate the project scores. It is bulk PBMC, not cell-resolved expression.

To check whether four previously prioritized genes have measurable counts, we calculated counts per million using each distributed file's sum of gene counts as denominator. This denominator is the sum of the released table, not independently verified sequencing depth; library preparation and filtering are not fully reconstructed. This post hoc descriptive screen is not a differential-expression analysis.

| Gene | Non-frail samples with nonzero counts | Frail samples with nonzero counts | Median CPM, non-frail / frail |
|---|---:|---:|---:|
| FOXP2 | 0/8 | 0/8 | 0 / 0 |
| SMAD5 | 8/8 | 8/8 | 166.3 / 185.5 |
| RBM6 | 7/8 | 7/8 | 28.0 / 28.5 |
| MEIS1 | 2/8 | 1/8 | 0 / 0 |

SMAD5 and RBM6 are measurable in most files. Zero detected FOXP2 counts and sparse MEIS1 counts concern this assay and set of samples only; they do not imply those genes are absent in other tissues or have no role in sleep or frailty. Group medians are not adjusted for race, sex, cellular composition, batch, sequencing quality or age, and the very small, selected sample does not support an inference that expression differs or is equal by frailty. Gene expression is not a genetic association, a variant-to-gene link, or evidence of causation.

Reproducibility: [aggregate results](stage6_handls_candidate_expression.csv) and `work/handls_candidate_expression.py`. The 16 distributed count files were extracted from the [GEO supplementary archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129534/suppl/GSE129534_RAW.tar), SHA-256 `3cb5d9dd11aa5956137703fa7a7404921d0d88ce4f598cd00851b96c29d88c5a`; `work/handls/matrix.gz` contains source metadata. Each per-file released count sum ranged from 123,441 to 533,134. Gene symbols were matched exactly; no transcript-level ambiguity or annotation version check was performed.
