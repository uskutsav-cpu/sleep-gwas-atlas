# Regional LD among the 22 shared-direction markers

The 22-marker count includes several markers in the same genomic neighborhoods. Among the 231 possible marker pairs, five are on the same chromosome and within 1 megabase (Mb) under the candidate coordinates. All five show some correlation in the 503 European-labeled 1000 Genomes reference samples; their estimated phased-haplotype r² ranges from 0.031 to 0.160. The other 226 pairs were not evaluated by this regional analysis, so this is **not** a complete pairwise LD screen.

| Marker and counted allele | Marker and counted allele | Distance (bp) | Phased r | r² |
|---|---|---:|---:|---:|
| rs1909650 G | rs9536347 T | 147,884 | +0.400 | 0.160 |
| rs12536335 A | rs12671330 C | 290,884 | +0.325 | 0.106 |
| rs12536335 A | rs17312861 G | 82,351 | +0.310 | 0.096 |
| rs7734045 C | rs7731706 G | 32,646 | +0.380 | 0.144 |
| rs12671330 C | rs17312861 G | 208,533 | +0.177 | 0.031 |

The three chromosome 7 markers are all labeled near or within FOXP2 in the project's positional mapping. Two of them, rs12536335 and rs12671330, are also among the five candidates nominally associated in the same direction across all eight inspected GWAS results. These observations do not establish three independent signals or FOXP2 as the causal gene. The chromosome 13 and chromosome 5 pairs also share regional LD despite having different nearest-gene labels. A gene label is not a locus boundary.

Positive r refers to the alleles counted by the shared-direction score: on the reference haplotypes, these alleles occur together more often than expected from their marginal frequencies. The underlying receipt uses the HFRS effect allele in some cases, so r was sign-flipped where needed to match the score alleles; r² is unchanged. This correlation does not mean the alleles cause either outcome. An unweighted sum counts each allele separately; correlated markers can therefore contribute overlapping information. Conversely, modest pairwise r² does not prove independent effects, exclude multiple causal variants, or justify removing one marker without an outcome-based evaluation.

## How this informs the next analysis

Keep the 22-marker count as the frozen exploratory prototype rather than silently pruning it. In a separate development and independent evaluation design, compare it with a prespecified LD-aware alternative and relevant genome-wide score benchmarks. Conditional analysis and fine-mapping with dense, ancestry-matched LD are needed to distinguish signals; colocalization with suitable molecular data is needed before making a shared-gene claim. The appropriate LD and score performance must be evaluated in the intended target population. No phenotype data or predictive evaluation are present here.

## Provenance and scope

The [source receipt](song124_22_nominal_reference_ld.json) records the candidate and reference metadata, all five pair results and exact allele matching. It uses 1,006 phased haplotypes from 503 European-labeled individuals in 1000 Genomes phase 3. The [score definition](unvalidated_22_shared_direction_score.csv) supplies the counted alleles used above. The original analysis assumed GRCh37 coordinates, and matching a biallelic VCF record does not resolve all source-study multiallelic representations. Pooled reference estimates have sampling uncertainty, and different cohorts can have different LD. The table does not measure trait correlation, variance in disease explained, or score accuracy.
