# Health ABC array subset: genotype-only reference sensitivity

The Health ABC array inventory matched 8 of the 22 markers in the equal-weight count. Using the existing score alleles and public 1000 Genomes Phase 3 reference genotypes, we compared that 8-marker subset with the complete 22-marker count. This is **not** a Health ABC analysis or validation against insomnia, frailty or future outcomes.

| Reference group | N | Correlation of 8-marker subset with complete count | SD of omitted 14-marker count |
|---|---:|---:|---:|
| EUR | 503 | 0.656 | 2.303 |
| AFR | 661 | 0.568 | 1.866 |

These correlations are moderate and demonstrate that the typed subset does not reproduce the full count in these reference samples. Correlations are invariant to offsets and do not establish similar absolute scores, prediction, or effects on outcomes. The 1000 Genomes groups are heterogeneous and are not substitutes for Health ABC's ancestry distribution or self-identified racial groups. Exact imputed and QC-passing marker coverage in Health ABC is unknown; the 8-marker subset is defined by a public platform annotation rather than its actual sample genotypes. No score should be shortened silently.

The [receipt](healthabc_partial_count_reference.json) records the exact 8 markers, input hashes and full results. The [script](../work/healthabc_partial_count_reference.py) uses one exact biallelic GRCh37 position, rsID and allele-pair match per score marker and complete diploid reference genotypes. Remote indexed VCFs were not fully checksum-verified, and project allele coding and multiallelic effects in the original GWAS remain separate uncertainties. Only aggregate results are saved.

## Separate 15-marker weighted prototypes

The array inventory matches 5 of the 15 markers. In the existing 503-person EUR reference genotype matrix, the partial and complete weighted sums correlate **0.535** for FI weights, **0.675** for inferred FFS weights and **0.540** for insomnia weights. These use the same five markers but different weights, not three independent validations; the FFS allele direction remains inferred. The omitted ten markers' weighted sums retain substantial variation. They are not Health ABC estimates and cannot establish predictive performance or justify ranking the weights. [Receipt and input hashes](healthabc_partial15_reference.json); [reproduction script](../work/healthabc_partial15_reference.py). Actual cohort imputation and QC may produce different coverage.
