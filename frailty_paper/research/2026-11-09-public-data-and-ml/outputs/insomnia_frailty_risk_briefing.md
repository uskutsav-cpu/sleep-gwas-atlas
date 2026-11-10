# Insomnia and frailty: risk comparison and candidate evidence

Research briefing for review • 27 October 2026

## Key findings

1. **FI and FFS have substantial genome-wide genetic correlation in the project's existing data.** Auxiliary LDSC estimates were 0.757 (SE 0.021; 95% CI 0.715–0.799) and 0.755 (SE 0.020; 95% CI 0.716–0.794) on two matched SNP sets. They use the same two frailty GWAS, so they are sensitivity checks, not independent replications. These are genetic correlations, not correlations between measured scores in the same people or percentages of shared individual risk.
2. **Insomnia correlates more strongly with this FI source than with this FFS source** in each of three source-specific comparisons. Why remains unknown. The UK Biobank component of the FI contains an insomnia item and a fatigue item; FFS includes exhaustion. Score content, participant overlap, age, ascertainment and other differences are plausible explanations that this analysis cannot separate.
3. **There is no validated 15-marker shared insomnia–frailty risk score.** A research-only prototype and 15-marker evidence table have been prepared. The markers were chosen from a published post-selected cross-trait list and their weights have not been evaluated against outcomes in an independent cohort. They do not yield an individual's calibrated risk.

## FI versus FFS

| Feature | Frailty Index (FI) | Fried Frailty Score (FFS) |
|---|---|---|
| Construct | Accumulation of multiple health deficits | Five physical domains: weakness, slow walking, exhaustion, low activity and weight loss |
| Source used here | Atkins et al. FI meta-GWAS: 164,610 UK Biobank + 10,616 TwinGene participants, European descent | Ye et al. quantitative 0–5 FFS GWAS: 386,565 UK Biobank participants, European descent |
| Content concern | UKB FI includes an insomnia and a tiredness item; TwinGene's FI is not identical | Exhaustion overlaps conceptually with fatigue; detailed original field/cutoff crosswalk remains unresolved |
| Study overlap | FI and FFS both include UK Biobank; exact participant intersection unknown | Same limitation |

Sources: [Atkins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/), [Ye et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/). The genetic-correlation evidence is reproduced from the project's saved native estimates; the FI–FFS point estimate and 200-block jackknife SE were recalculated from those outputs, not refit from raw GWAS. See `direct_fi_ffs_genetic_correlation.csv` and `fi_ffs_risk_followup.md` for checksums and input details.

The UK Biobank 49-item FI includes both insomnia and tiredness. A related published 44-item SALT index underlying the Swedish component lists neither explicitly, but the Atkins GWAS-specific item table remains unverified. I found a discrepancy between a public FI derivation script and a later published table for one tiredness response category, so the exact original Atkins coding needs verification before reconstructing that GWAS phenotype. Ye’s FFS includes exhaustion, but its exact original item implementation also remains to be verified. See `fi_ffs_component_crosswalk.md`.

## How similarly do FI and Fried classify people?

In a separate Toledo cohort of 1,563 older adults with both measures, a reconstruction of the published eight-category table gives 103 classified as frail by both, 125 by FI only, 18 by Fried only, and 1,317 by neither. Among those frail by either measure, 41.9% were frail by both; overall agreement was 90.9% because most were non-frail by both. These are threshold-based classification results in a separate cohort, not a continuous score correlation or a genetic correlation. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC13264346/); see `fi_ffs_same_cohort_comparison.md` for the derivation and mortality-model limitations.

A Spanish study’s Figure 2 reports a **Spearman correlation of 0.542** between its 27-item FI score and count of Fried phenotype criteria. It is a small cross-sectional study with different inclusion criteria and frailty definitions, not an estimate for the project's GWAS participants. [Article and figure](https://link.springer.com/article/10.1007/s11357-025-01931-z/figures/2).

A separate Canadian cohort of 28,527 adults reported a **phi correlation of 0.38 between binary FI and Fried frailty classifications**, with 92.2% overall agreement. Their 93-item FI and clinical cutoffs differ from our GWAS definitions. This is a same-participant categorical association, not the correlation of continuous scores or the genetic correlation. [Pasquet et al.](https://link.springer.com/article/10.1007/s40520-024-02706-w); see `fi_ffs_same_cohort_comparison.md`.

## Insomnia correlations with each measure

| Insomnia source | Insomnia–FI r_g | Insomnia–FFS r_g | FFS minus FI | 95% CI for difference |
|---|---:|---:|---:|---:|
| UK Biobank symptom measure | 0.640 | 0.483 | −0.156 | −0.213 to −0.100 |
| FinnGen clinical insomnia | 0.426 | 0.245 | −0.181 | −0.260 to −0.102 |
| MVP European clinical insomnia | 0.549 | 0.398 | −0.151 | −0.220 to −0.082 |

The difference uses matched SNP support and paired jackknife covariance. All three pass the existing conservative 396-test Bonferroni overlay. The same FI and FFS GWAS are reused, and the comparison has not tested cross-source differences. The original UKB symptom contrast had paired-uncertainty validation; the two clinical contrasts additionally passed the saved three-trait matrix checks. These are descriptive comparisons of specified data, not independent replication of frailty, evidence of causality, or proof that removing an FI item would eliminate the difference.

## Fifteen-marker risk evidence

The 15 rows come from Song et al.'s published cross-trait list, taking the 30 lowest reported CPASSOC p-values and retaining those available in the three project inputs after existing QC. This is a **post-result, availability-driven research set**. It is not LD-pruned or an independently validated score. Nine of the 15 rsIDs have more than two alleles listed in the Ensembl reference. The FI and FFS inputs overlap, exact overlap with insomnia sources remains unknown, and the raw non-HapMap3 FinnGen/MVP markers have not undergone full project QC. Source p-values below are **unadjusted**. A p-value describes a statistical association, not a person's chance of disease or a causal effect.

| Variant | Counted allele | FI p | FFS p | UKB insomnia p | FinnGen insomnia p | MVP insomnia p | All five p<0.05? |
|---|---|---:|---:|---:|---:|---:|---|
| rs113851554 | T | 0.199 | 0.19 | 2.35e−21 | 1.16e−10 | 4.22e−15 | no |
| rs8067974 | G | 7.77e−7 | 0.027 | 1.52e−11 | 0.223 | 0.261 | no |
| rs1909650 | G | 6.32e−7 | 0.097 | 9.02e−10 | 0.102 | 1.58e−5 | no |
| rs11152363 | A | 0.0871 | 0.032 | 2.18e−7 | 0.018 | 0.144 | no |
| rs1032524 | T | 3.53e−4 | 1.00e−4 | 2.56e−6 | 0.408 | 0.161 | no |
| rs323509 | A | 1.80e−6 | 8.00e−10 | 6.83e−9 | 1.58e−4 | 0.0509 | no |
| rs78206187 | G | 6.94e−4 | 0.75 | 6.30e−8 | 4.05e−4 | 3.19e−4 | no |
| rs1942262 | A | 0.00594 | 0.003 | 2.57e−6 | 0.443 | 0.574 | no |
| rs4886140 | A | 7.89e−5 | 0.45 | 1.05e−7 | 0.166 | 0.0154 | no |
| rs12536335 | G | 2.63e−6 | 0.0017 | 7.65e−7 | 0.00625 | 4.57e−8 | yes |
| rs35881094 | G | 0.602 | 0.2 | 4.37e−5 | 0.254 | 0.133 | no |
| rs17550821 | A | 0.886 | 0.072 | 2.13e−8 | 0.247 | 0.413 | no |
| rs7734045 | C | 0.00545 | 0.015 | 1.59e−6 | 0.0352 | 0.0184 | yes |
| rs11790060 | C | 0.00106 | 6.00e−5 | 2.20e−5 | 0.162 | 0.0627 | no |
| rs10865954 | T | 0.00181 | 7.80e−7 | 2.65e−6 | 0.0741 | 0.333 | no |

Two variants—rs12536335 and rs7734045—show p<0.05 in all five sources with the same aligned sign. This is a **descriptive observation**, not a multiple-testing-adjusted discovery: the markers were selected after earlier results, the data are not fully independent, and LD and allele representations need validation. None of the 15 satisfies p<0.05/15 in all five sources; that numerical benchmark is also not a valid post-selection corrected claim. Non-significance is not evidence of no effect.

A cohort-specific follow-up adds a caution about rs12536335: after applying the project’s independently audited TwinGene allele-direction adapter, the counted-G FI coefficient is negative in the UKB component and positive in TwinGene. Their phenotypes and samples differ, and this post-selected observation does not explain the difference or validate a shared variant. See `candidate_twingene_sensitivity.md`.

A cautionary example is rs113851554: despite strong insomnia p-values in all three sleep sources, its FI p is 0.1987 and its FFS p is 0.19 in the project inputs. Its sleep evidence must not be presented as a confirmed shared frailty association.

The 124 original Song candidates and the full 15-marker estimates, SEs and approximate intervals are preserved separately. Song's cross-trait approach does not prove that each listed variant independently affects both traits or that a shared causal variant exists. [Song et al.](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full). The project's own FI×sleep LAVA run completed all 29,940 slots, but all 12 prespecified failure gates failed; those results cannot be used to establish or exclude local sharing.

The 15 candidates now have GRCh37 nearest-gene-body annotations in `candidate_nearest_gene_annotations.csv`; the annotation method and cautions are in `candidate_nearby_gene_notes.md`. For example, rs12536335 overlaps the FOXP2 gene body, but proximity does not establish FOXP2 as the effector gene or resolve the opposite FI cohort signs.

## Score status and what the published scores show

The `unvalidated_15_variant_score_prototype.csv` gives allele-aligned FI, FFS and insomnia coefficients for three separate raw weighted sums. A later technical check computed the three sums in 503 public 1000 Genomes reference individuals without outcome data; it has no validated scaling, calibration, discrimination, missing-genotype handling, LD pruning, or independent outcome evaluation. See `reference_15_score_diagnostic.md`. The three phenotype-specific weights are on different scales and must not be added into a single score. The prototype is not for clinical or individual risk use.

Published genome-wide scores exist: [PGS005228](https://www.pgscatalog.org/score/PGS005228/) was derived from FI and [PGS005229](https://www.pgscatalog.org/score/PGS005229/) from FFS. In the published study these were evaluated against FI in older-adult cohorts; reported performance varies by cohort. Their whole-score results cannot be assigned to the 15-marker subset and do not establish a joint insomnia–frailty score. [Foote et al.](https://www.nature.com/articles/s41588-025-02269-0)

A later primary-supplement audit found that several FI-outcome performance values in PGS Catalog do not agree with the paper’s Supplementary Table 44; the reason is unresolved. The evaluated ELSA and PISA FI definitions include sleep-related items; the listed LBC1936 FI items contain no explicit sleep or fatigue label under a keyword check. Table 44 also reports FI- and FFS-derived score associations with categorical physical frailty in LBC1936. See `fi_ffs_pgs_published_validation_audit.md` before quoting or comparing performance metrics.

A file audit found that all 2,379,429 rows in the original FI PGS file repeat the effect allele in the `other_allele` field. This does not prove incorrect weights or invalidate the published study, but the two columns cannot be treated as a distinct allele pair. In a deterministic sample, 2,322 of 2,323 score rows matched a project FI variant at the same position, and every matched counted allele agreed with one allele in that FI input. This supports the counted allele within that filtered sample only; it is not a validation of all source rows or predictive performance. The original score was not changed.

A published sparse insomnia-only benchmark, [PGS003859](https://www.pgscatalog.org/score/PGS003859/), contains 129 variants and was associated with self-reported insomnia symptoms in a UK Biobank evaluation (OR 1.08 per score SD, 95% CI 1.07–1.09). It is not a shared frailty score. An exact-rsID lookup finds one overlap with our 15 candidates, rs113851554; the score file lacks coordinates, so this does not establish complete locus or allele harmonization. See `published_insomnia_score_comparison.md` for discovery provenance, limitations, and a proposed validation design.

## What would resolve the main uncertainty

* **Why FI and FFS differ:** obtain appropriately authorized, same-participant FI and FFS data; reconstruct source-authentic or clearly named alternative scores; compare full FI with FI excluding the insomnia item and then insomnia plus fatigue, using covariance-aware estimates and a fixed analysis plan. Exact source scoring and adequate power must be verified first.
* **A genuine shared score:** define the outcome and time horizon, select and estimate weights with an LD-aware, prespecified procedure, then evaluate calibration and incremental performance against existing FI, FFS and insomnia scores in independent genotyped cohorts. No such evaluation is available here.
* **Candidate loci:** retain the 15 and 124 tables as hypotheses. A separately justified, QC-passing analysis with dense inputs, ancestry-matched LD, overlap review and fine-mapping/colocalization is required before claiming a shared causal variant. The failed LAVA gate remains in force.

## Supporting files

* `prototype_variant_five_trait_evidence.csv` — 15 variants, five source effects and intervals where SEs exist.
* `prototype_five_source_descriptive_triage.csv` — transparent numerical screen; not adjusted significance.
* `published_insomnia_fi_124_cross_trait_candidates.csv` — complete transcribed published candidate list.
* `unvalidated_15_variant_score_prototype.csv` — exploratory coefficients and allele warnings, not a clinical score.
* `published_fi_ffs_pgs_source_audit.json` and `fi_pgs_allele_sample_summary.json` — published-score provenance and allele audit.
* `fi_ffs_risk_followup.md` — full methods, source limitations, example calculations and links.

## Full FI score counted-allele check

A later full-file audit compared all 2,379,429 published FI score records with the project's harmonized FI GWAS. **2,377,040 matched by rsID and GRCh37 position, and the score's counted allele matched a GWAS allele in every matched record.** The remaining 2,389 were absent from that filtered GWAS input, including 2,359 in its excluded MHC interval; the other 30 have unresolved absence reasons. This is a within-source allele-string check, not validation of the score's weights or predictive performance. The duplicate `other_allele` field in the original score remains unresolved. See `fi_pgs_full_counted_allele_audit.json` for the complete counts and source checksums.

For the FFS score, a parallel full-file audit found that 2,373,832 records matched the article-linked FFS input by rsID and coordinate, with the counted allele matching its reported `A1`. The remaining 388,036 have no score-side coordinates; they matched the source by rsID and allele only, which is weaker evidence. These source checks do not validate predictive performance. See `ffs_pgs_full_counted_allele_audit.json` and `ffs_pgs_missing_coordinate_rsid_audit.json`.

A subsequent weight-orientation audit found that published score weights were positively correlated with the allele-aligned marginal effects in their respective underlying discovery GWAS (FI r=0.898 over 2,377,040 matched rows; FFS r=0.730 over 2,373,832 matched rows). This is a within-source diagnostic, not an independent prediction test or a comparison of which score performs better. LD-aware SBayesR weights can differ from marginal effects. See `fi_pgs_weight_orientation_audit.json` and `ffs_pgs_weight_orientation_audit.json`.

## Additional full HFRS lookup

The project already contains the complete checksum-matched FinnGen HFRS GWAS, GCST90707291. All 15 candidate allele pairs matched by rsID and mapped position; five had nominal HFRS p < 0.05 with the same estimated direction as FI, subject to selection and cohort-overlap limits. The stored exploratory UKB-insomnia versus FinnGen-HFRS LDSC estimate is rg=0.439 (SE 0.026); its HFRS univariate intercept warrants further review. See `candidate_hfrs_full_gwas_review.md`. An earlier statement that the full HFRS GWAS could not be located was corrected.
