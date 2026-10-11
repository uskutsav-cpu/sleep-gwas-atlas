# Independent scientific source admission: one exploratory FinnGen–MVP transport analysis

Reviewed 11 October 2026. This review generated **zero genetic-correlation fits**. Historical V4/V7 results, thresholds and source files were not changed. The input contract and bounded source diagnostic are companion JSON files in this directory.

## Decision

**Conditional GO for one separately frozen exploratory analysis:** FinnGen R13 registry insomnia (`F5_INSOMNIA`) × MVP genetically inferred European ancestry abdominal pain (`GCST90479148`, phecode 785). Existing verified derivatives can be used under explicit CI95/log-Wald, effective-N and HapMap3 fallback assumptions. No fundamental invalidity of these assumptions was established by the reviewed evidence. Their remaining uncertainty prevents treating this experiment as confirmatory source admission or fully independent replication.

This experiment asks whether a selected historical sleep–abdominal-pain association is compatible with a different sleep-side source and ascertainment construct. FinnGen supplies a new sleep source; the MVP outcome was used previously and remains identical. Therefore this is **source-and-construct transport robustness**, with reused outcome data and historical selection. A favorable result cannot establish a new discovery, an independent outcome confirmation, causality, or a statistically significant difference between clinical and questionnaire sleep constructs.

The V4 frozen feasibility protocol authorized preprocessing and one heritability diagnostic only. Its original source admission requirements remain unresolved. A new, explicitly qualified exploratory protocol must be frozen before this new fit; its permissive assumptions must not be retrospectively presented as passing V4's stricter requirements.

## Source-specific schema contract

| Contract field | FinnGen insomnia | MVP abdominal pain |
| --- | --- | --- |
| Source | R13 `F5_INSOMNIA`, generation 1777989563097164 | `GCST90479148`; PMID 39024449; GIA EUR, phecode 785 |
| Raw coordinates | GRCh38 | GRCh38, 1-based; public metadata describes unharmonised source |
| Effect orientation | ALT is effect allele; supplied `beta` is log OR | `effect_allele` is effect allele; supplied `odds_ratio` |
| Existing signed Z | `beta/sebeta`, aligned to retained HM3 allele orientation | `log(OR)/SE_logOR`, with `SE_logOR=(log(ci_upper)-log(ci_lower))/(2*1.959963984540054)` |
| Z qualification | Supplied effect SE; source P is retained as diagnostic evidence | Assumes the exported interval is a 95% normal/Wald interval on log-OR scale; no supplied SE available in the bounded sample |
| N convention | Assumed constant `4/(1/51643+1/446273)=185146.70377332723`; serialized 185146.703773 | Existing constant effective N from 89,262 cases and 299,053 controls: 274973.347782; raw total 388,315 |
| N qualification | Endpoint case/control totals, not verified SNP-wise analysed N | Bounded raw sample has those counts consistently; whole-body SNP-wise N equivalence not newly established |
| Imputation/QC route | No INFO column; HM3 fallback, MAF >0.01, nonambiguous alleles, coordinate/allele mapping and duplicate exclusions | Supplied `r2>0.9`, MAF >0.01, nonambiguous HM3 alleles; historical pipeline unchanged |
| Source derivative | 1,156,359 SNPs; current compressed SHA matches receipt | 826,029 SNPs; current compressed SHA matches receipt; existing independent replay matched original ordered content/QC |
| LD reference | Existing 1000 Genomes EUR LD scores and weights; Finnish LD mismatch remains a limitation | Same pinned reference; broad European matching does not prove absence of mismatch |

FinnGen's documented file schema identifies ALT, log-OR beta, effect SE and GRCh38 coordinates; neither INFO nor N is a source column. Batch annotation is separate and does not certify endpoint-specific INFO for the retained derivative. R13 used regenie 3.3 with relatedness handling, approximate Firth testing and demographic/genotyping covariates. These facts support the signed-effect route but do not establish SNP-wise N or elimination of residual confounding. [FinnGen data schema](https://finngen.gitbook.io/documentation/data-description), [R13 GWAS methods](https://finngen.gitbook.io/documentation/methods/phewas).

The registry endpoint includes F51.0/G47.0 diagnoses and excludes sleep-disorder controls. The inspected endpoint does not impose a minimum-event rule. Registry diagnosis is not interchangeable with questionnaire symptoms, chronicity, an identical observation horizon, or standardized clinical adjudication. Official endpoint information supplies genotype-QC counts 51,643/446,273; this differs from the R13 study-wide total and is the relevant count convention. [R13 insomnia endpoint](https://r13.risteys.finngen.fi/endpoints/F5_INSOMNIA).

MVP's primary study identifies SAIGE-GPU mixed-model association analysis and GIA-defined ancestry strata. It does not, in the material inspected, establish this export's precise confidence-interval construction. The exact public source metadata confirms build, EUR total N, accession and MD5; a GWAS-SSF label is not proof that an exported interval is 95% Wald. [MVP source study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/), [exact accession metadata](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479148/GCST90479148.tsv.gz-meta.yaml), [GWAS-SSF specification](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.0.2.tex).

## Bounded new diagnostic and its proper interpretation

`MVP_bounded_CI_schema_diagnostic.json` records a first-1000 eligible-row inspection from the first 6,966 raw rows. Eligibility used finite OR/CI, r2>0.9, MAF>0.01 and nonambiguous alleles; it was not an exact HM3 sample, not a genomic random sample, and did not select favorable P-values. All 1,000 supplied SE values were `#NA`. Their counts were consistent and their intervals nearly log-symmetric (maximum log-center discrepancy 0.000715882449).

27/1,000 rows had >10% relative disagreement between supplied P and normal P reconstructed from existing CI-derived Z. This is **a bounded schema/calibration observation, not a validated source defect**. Mixed-model score/SPA P-values need not equal effect/Wald P-values; rounding is another possible contributor. Exact attribution was not established. Do not impose equality as an invented validity criterion, replace the historical Z route after viewing a result, or describe CI95 semantics as verified. The confidence-level and construction assumption must accompany every result from this route.

The LDSC authors recommend HapMap3 filtering when INFO is absent and permit assumed constant N when SNP-wise N is unavailable. Thus these absences do not universally invalidate LDSC; they warrant explicit qualification here and do not satisfy a frozen requirement to verify INFO>0.9. The existing European reference is the standard route, with Finnish LD limitations retained. [Official LDSC input guidance](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

Constant-N scale assumptions limit interpretation of absolute observed h². Normalization can remove simple constant scale changes in rg, but cannot cure SNP-wise N variation, test-statistic miscalibration, ascertainment or reference mismatch. The fit must not be labelled population/liability heritability, and clinical diagnostic rates must not be treated as population prevalence.

## Existing executed evidence: no repeats needed

FinnGen's existing native h² is 0.05353575912413689 (SE 0.0040705198925973796; Z 13.15206915497371), nominal interval 0.04555768673629218–0.0615138315119816 under the assumed N. Intercept is 1.077872498550233 (SE 0.009015015588198625), mean χ² 1.2705273128003078, λGC 1.2387232869346392. The original necessary Z≥4/intercept≤1.2 conditions pass; they are not sufficient source-validity or pair-power tests. Existing MVP h² is 0.054710911106158965 (SE 0.003506309835863687), intercept 1.132466986433778 (SE 0.010561499362995507).

The actual input paths, current one-time compressed input hashes, pinned LDSC code hashes, interpreter, all-22-chromosome reference presence, and **actual prior h² commands** are in `qualified_native_input_contract.json`. No raw GWAS was redownloaded, no million-row calculation repeated and no new execution controller created.

Exact evidence paths, all under `/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/`:

- `new_source_feasibility/finngen_R13_F5_INSOMNIA/finngen_R13_insomnia_acquisition_verified_v4_2.json`
- `new_source_feasibility/finngen_R13_F5_INSOMNIA/pipeline_replay_v3_8/derived/insomnia.preprocessing.json`
- `new_source_feasibility/finngen_R13_F5_INSOMNIA/pipeline_replay_v3_8/receipts_v4/preprocessing.execution_receipt_v3_8.json`
- `new_source_feasibility/finngen_R13_F5_INSOMNIA/observed_h2_v3_9/fits/insomnia_observed_h2.full_precision.json`
- `new_source_feasibility/finngen_R13_F5_INSOMNIA/observed_h2_v3_9/receipts_v4/observed_h2.execution_receipt_v3_9.json`
- `validation_pipeline_replay_v4/receipts_v4/gwas_catalog_GCST90479148.source.json`, `.adapter.json`, `.comparison.json`
- `native/original_runner_support_v4_2/native/validation_reproduction_v1/validation_h2_gwas_catalog_GCST90479148.full_precision.json` and `.execution_receipt.json`.

## Exact remaining essential gates for this experiment

1. Freeze one pair, one source release, the existing Z/N/QC assumptions, stock allele alignment, unconstrained h² and genetic-covariance intercepts, 200-block jackknife, output location and complete-result reporting before execution. Declare the new finite family of one as exploratory; historical results informed its selection.
2. Run only the new two-trait rg, using the verified inputs and pinned existing pipeline. Preserve resource guards and exclude the separate IBD/Crohn/UC source orientation issue from this experiment.
3. Review the actual estimator result: intersection count/allele handling, finite estimate and SE, nominal interval, intercepts, warnings, native jackknife evidence and successful completion. A numerical failure or nonpositive h² making rg uninterpretable is a no-go for interpretation; do not rescue a favorable point estimate by suppressing warnings or changing methods.
4. Report regardless of sign or significance, as qualified source/construct transport. Independent reviewers must assess the actual estimate, uncertainty and source qualifications. Do not repeat completed h² or rerun unrelated historical analyses.

Leave overlap intercepts unconstrained because participant intersections are not quantified. LDSC can model overlap-related intercept terms; this does not prove cohort disjointness or supply independence from historical selection. Do not use CI overlap or a point-estimate change as a calibrated clinical-vs-questionnaire difference test; shared-outcome covariance remains uncalibrated. [LDSC FAQ](https://github.com/bulik/ldsc/wiki/FAQ), [cross-trait LDSC methods](https://www.nature.com/articles/ng.3406).

For this selected one-pair experiment, show nominal exploratory P and optionally the stricter historical sensitivity cutoff 0.05/217=0.0002304147465437788. The latter does not by itself define or control the full source/method/selection family, remove winner's curse, create conditional post-selection validity, or make reused MVP independent. For a truly fixed full family with valid marginal P-values, Bonferroni can control unconditional FWER despite dependent tests; that is a different statement from asserting this selected new test has a proven confirmatory family.

Precision cannot be promised from the passing univariate diagnostics. An analytic normal approximation at two-sided 0.05/217 and 80% power would require SE≤0.06630 to detect a true |rg|=0.30. With SE=0.10, the corresponding illustrative detectable magnitude is about 0.4525 and a nominal 95% interval has half-width 0.196. These are illustrative prospective requirements, **not estimates of actual pair power** or forecasts based on the historical selected effect.

## Alternatives and publication consequences

No inspected alternative closes the independent two-trait gap now. The All of Us public browser is a lead, while full Hail tables require controlled-tier access; current official guidance flags corrections to v8 missing nonsignificant rows and deprecation of some v7 tables. No exact public, full insomnia source with usable schema, counts and adequate precision was verified here. Different cohort names also do not establish participant disjointness. Restricted access, forms, investigator contact and bulk retrieval were not pursued. [All of Us current All-by-All documentation](https://support.researchallofus.org/hc/en-us/articles/27049847988884-Overview-of-the-All-by-All-tables-available-on-the-All-of-Us-Researcher-Workbench).

A 2025 preprint already studied clinical insomnia versus ordinary sleep variation using UKB, MGB, MVP, FinnGen and All of Us, with Estonian Biobank replication. It is relevant prior art, rather than an automatically admissible independent source: its meta-analysis contains historical cohorts and exact independent public inputs were not established by this review. The paper's existence narrows novelty claims around clinical-vs-questionnaire sleep genetics. [Kunorozva et al. primary preprint](https://doi.org/10.1101/2025.11.21.25340447).

This proposed single analysis can add genuinely new **project evidence** about transport to a new sleep-side source once executed; it does not guarantee novel field-level discovery or acceptance at Cell Genomics. The present review adds a bounded source-schema diagnostic and an evidence-backed admission decision, not a new biological association. Essential next work is the one qualified fit and its actual scientific review, followed by figures, use-case interpretation, machine-readable qualified release and manuscript materials. Additional controller development, repeat hash campaigns, repeated completed h² and new exhaustive audits are optional and unwarranted absent a concrete defect. Exact CI export documentation, independently sampled validation or population-matched LD can strengthen later claims, but should not silently convert this qualified transport experiment into replication.
