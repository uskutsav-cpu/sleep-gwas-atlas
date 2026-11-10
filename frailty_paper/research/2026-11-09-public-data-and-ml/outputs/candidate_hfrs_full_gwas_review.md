# All 15 candidates in the complete FinnGen hospital-frailty GWAS

**Correction:** an earlier search note said the full HFRS GWAS had not been located. The repository already contains the complete [GWAS Catalog file GCST90707291](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90707001-GCST90708000/GCST90707291/) for Mak et al.'s FinnGen R12 hospital frailty risk score (HFRS). I scanned all 21,294,560 variant records and reproduced its catalog metadata MD5, `38d9434f50a4493def9d2c9d1cc1b4b3`. The source metadata identifies GRCh38, 500,737 participants, SAIGE, and effect/other allele fields. This is the **full** HFRS phenotype; it is not the dementia-excluded summary file.

All 15 candidate rsIDs had an exact match at their independently mapped GRCh38 position and the same pair of allele strings. Betas were reversed when the source effect allele was the candidate's comparator; the full result, source alleles, p-values and match counts are in `candidate_finngen_hfrs_full_gwas_lookup.csv`. Nine of the rsIDs are multiallelic in the broader Ensembl reference, so matching the two source strings does not independently settle every original GWAS coding detail.

Five candidates had nominal HFRS p < 0.05, and all five had the same estimated direction as the project FI meta-analysis for the candidate's counted allele:

| Variant; counted allele | Nearby gene-body label | HFRS beta per counted allele | SE | HFRS p |
|---|---|---:|---:|---:|
| rs1909650 G | PCDH8P1 | +0.0452 | 0.0127 | 0.000386 |
| rs323509 A | RP11-6N13.1 | +0.0398 | 0.0139 | 0.00425 |
| rs4886140 A | RPP40P2 | −0.0394 | 0.0129 | 0.00226 |
| rs12536335 G | FOXP2 | −0.0728 | 0.0126 | 7.85×10⁻⁹ |
| rs7734045 C | AC009014.3 | +0.0491 | 0.0137 | 0.000325 |

The source beta units are on the HFRS outcome scale, not odds ratios or absolute probabilities; these coefficients cannot be added to FI, FFS or insomnia weights. Four p-values (all except rs323509) are below the descriptive 0.05/15 threshold, but these candidates were selected after cross-trait discovery, input data may overlap, and this is not a prespecified independent multiple-testing analysis. Do not label them confirmed shared causal variants. For example, rs1909650 and rs4886140 lack nominal FFS evidence in the reviewed source, while rs12536335 and rs323509 have discordant UKB/TwinGene FI component directions under an inferred adapter. Nearby genes are positional labels, not established causal genes.

The HFRS and FinnGen insomnia results use overlapping FinnGen releases, so the new HFRS associations and FinnGen insomnia results are **not independent cross-trait replications**. The HFRS is a different construct from FI and FFS and includes broad codes that can capture drowsiness or vital exhaustion (`hfrs_sleep_overlap_audit.md`). Participants and study effects may differ across UKB and Finnish sources. The ten candidates without nominal HFRS significance must not be treated as no-effect or failed replication on p-values alone.

## Existing genome-wide result in the project

The repo also contains an exploratory LDSC result for UKB insomnia versus this FinnGen HFRS: genetic correlation **0.439**, SE 0.026, approximate 95% CI 0.388–0.490, with 1,111,759 regression SNPs. The stored run passed its project h² precision gate and reports an unconstrained cross-trait intercept of −0.00037 (SE 0.00603). The separate HFRS univariate LDSC report has intercept 1.191 (SE 0.012), which warrants investigation of model assumptions and residual confounding; a precision-gate pass alone does not establish unbiased estimates. Cohorts are described as UKB versus FinnGen, but exact participant overlap is unknown. This result is not a genetic correlation for FI versus FFS, not a local shared-signal result, not a causal effect and not independent replication of exactly the original frailty phenotypes. It does not bypass the failed FI×sleep local-correlation gates. Project source: `nature_communications_continuation_v6/repro_inputs/clinical_validation_v5_1/tables/NEW_VALIDATION_FAMILY_RESULTS.tsv`, `insomnia_x_HFRS.exact.json`. The new lookup script is `work/hfrs_15_full_lookup.py`.
