# Executed robustness diagnostics

Audit date: 2026-10-08. All analyses below are retrospective source-free diagnostics. They were not prespecified confirmatory experiments. The original 1,200-pair BH and 217-candidate Bonferroni families remain unchanged.

`tables/sensitivity_results.tsv` contains 281 executed diagnostic rows: three normal-P/BH precision summaries, two conservative printing-bound summaries, one stricter discovery Bonferroni summary, three measurement strata, eighteen phenotype-domain strata, two definition strata, twelve sleep-specific intercept summaries,205 pair-level covariance-grid calculations, five covariance-grid aggregate summaries,28 profile-similarity comparisons and two negative-outcome summaries. The generator derives these values from original tables and retained replication logs.

Precision checks preserve all 603 discovery hits and all 23 replication successes. The stricter exploratory discovery Bonferroni cutoff, 0.05/1200, retains 337 comparisons. Changing correction methods after observing outcomes would be a new analysis and is not presented as an original decision.

Measurement stratification keeps the same global BH labels:

| Sleep measurement | Tests | BH-significant |
|---|---:|---:|
| Self-report/questionnaire,8 traits | 800 | 474 |
| Accelerometer,3 traits | 300 | 61 |
| Clinical register sleep-apnea endpoint,1 trait | 100 | 68 |

Twenty-two of 23 replication successes involve self-report traits and one uses register-defined sleep apnea; none uses any of the three accelerometer traits. These are descriptive differences in findings and measurement coverage, confounded by sample size, heritability, phenotype and replication availability. No formal test of superiority between measurement groups is justified. The audit does not infer absence of actigraphy associations from smaller significant counts.

The58 source-available candidates comprise 37 exact-definition and 21 comparable-definition pairs, including QC-ineligible sources. They yield 18 and 5 successes respectively under the same 0.05/217 threshold. These proportions are not a randomized or source-matched comparison. Comparable definitions remain visibly labeled in every output; excluding them retrospectively is not the primary analysis.

For the eight replicated external phenotypes, all 28 Pearson similarities between their vectors of 12 discovery sleep-rg estimates range from 0.7752 to 0.9851. These high profile similarities motivate interpreting overlapping gastrointestinal, smoking/pulmonary, pain and musculoskeletal manifestations cautiously. They are **not** external-phenotype genetic correlations: correlated sleep measurements, common cohorts, correlated estimation errors and selection of positive results can contribute. No independence-adjusted cluster P, effective number of diseases, conditional relationship, causal mediator or replacement correction denominator was computed. A valid external-phenotype genetic/sampling covariance matrix is needed to distinguish clusters statistically.

Sleep-specific cross-trait-intercept summaries report nominal intercept P diagnostics and largest absolute values. They quantify recorded diagnostics without converting them into a pass/fail rule or proof of independence. In particular, a small intercept does not exclude shared sleep input, structure or correlated errors. Covariance-grid calculations are detailed in report05.

Smoking/adiposity conditional genetic-correlation or mediation analyses are blocked by missing eligible multivariate genetic and sampling covariance inputs. No unadjusted correlation was relabeled as a direct effect. A scientifically valid follow-up would freeze traits and source families, establish allele/build/ancestry eligibility and estimate joint covariance before fitting the selected model.

Source-free numerical robustness is `PASS_WITH_QUALIFICATIONS`; native sample-overlap correction, conditional smoking/adiposity analysis and formal phenotype-cluster covariance analysis are `BLOCKED_EXTERNAL_INPUTS`.
