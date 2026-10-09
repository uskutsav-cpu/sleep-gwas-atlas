# Adversarial biological interpretation review v1

Client review date: 2026-10-08, America/Chicago. Phase 10 independent automated reviewer. This is a technical claim audit, not manuscript text or human specialist approval. No original, frozen result, SSD file or source GWAS was modified. No new MVP insomnia association results were inspected; no GWAS, local, QTL, causal or clinical analysis was run. The accompanying standard-library script inventories existing evidence and independently recalculates descriptive profile similarities; successful assertions are not biological validation.

**Disposition: qualified global common-variant sharing is the strongest supported biological interpretation. No effector gene, tissue, pathway, shared causal variant, causal direction, mediation, treatment benefit or clinical decision value is established. The new-primary-question NO-GO remains scientifically appropriate.** Conditional clinical-insomnia transport validation could improve generalization after its frozen gates, but would not remove these boundaries.

## Evidence inventory and current versus historical status

| Evidence checked independently | Observed inventory | Interpretation boundary |
|---|---:|---|
| Historical extension claim ledger | 1,200 rows, all `OBSERVED_GLOBAL_GENETIC_CORRELATION` | Every row permits global association only |
| Historical threshold-positive external checks | 23 pairs, 7 sleep traits, 8 outcome labels | All 23 reuse their discovery sleep GWAS |
| Positive-pair measurement composition | Insomnia 7, short sleep 6, long sleep 5, snoring 2, apnea 1, napping 1, reported duration 1 | 22 questionnaire-based pairs, 1 registry-based pair, 0 device-based pairs |
| Outcome definition labels | 18 `EXACT`, 5 `COMPARABLE_WITH_DOCUMENTED_DIFFERENCES` | Source-curated concepts; complete ascertainment equality not certified |
| Local, pleiotropic, colocalized-gene and mechanism fields | Each contains `NA_BLOCKED_UPSTREAM` for all 23 ranked pairs | Unavailable experiments; no zero-locus or no-mechanism conclusion |
| Eight-outcome sleep-profile similarities | 28 pairwise Pearson values, 0.7752–0.9851; independent recomputation agrees | Correlations of 12 sleep-rg point estimates; not outcome-outcome rg or validated clusters |
| New native historical insomnia–BMI control | rg ratio 0.17229925855319342 | Reproduction control; no mediation or clinical validation |
| New independent clinical-insomnia validation | Source and protocol review only in this audit | No new positive, null or eligible-rg result adjudicated here |

Inputs and exact hashes are in `biological_interpretation_v1.json`; structured claim dispositions are in `biological_interpretation_v1.tsv`. This reviewer accepts the parent/reproduction reviewers' byte recovery and frozen-family arithmetic as separate findings. Exact recovery of 18 artifacts and one native pair cannot establish complete native scientific reproduction or clinical meaning.

Historical `discovery_extension/final_report.md` retains the headline “independent replication” and the stronger independent-EUR wording in its claim boundary. Those claims are superseded by source reuse evidence and the unified frozen protocol. Preserve the historical file, but do not propagate that wording to unified interpretations. The prior submission package correctly qualifies outcome-side validation. Its missing-input statements describe the earlier detached-SSD state, not this recovery. Clinical and biological restrictions persist even when original inputs become accessible.

## Attempts to falsify clinically meaningful conclusions

### B1 — Strong symptom–disease sharing does not identify a sleep mechanism

Insomnia–abdominal pain has discovery rg 0.5591 and historical outcome-side rg 0.4098; insomnia–reflux has 0.5216 and 0.3640. These are substantial common-variant association estimates under the named sources and model. They cannot distinguish direct sleep-related biology from broad symptom liability, genetic influences on smoking/adiposity, correlated illness liabilities, healthcare capture or selection. These are competing explanations, not explanations demonstrated by this audit. Replacing the external outcome leaves the same UKB sleep phenotype in both estimates and cannot establish independent clinical-insomnia biology.

Inherited influences on a behavioral or disease determinant can contribute to both GWAS effect profiles; that does not require the measured sleep trait to cause the outcome. Conversely, a correlation alone cannot exclude a causal relationship. Clinical labels, large rg and small P do not decide between these possibilities. The [primary cross-trait LDSC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4797329/) defines rg as covariance standardized by SNP heritabilities and models sample overlap in the intercept under its assumptions. It supplies a global sharing estimand, not treatment efficacy, individual risk or variant-level localization. A small intercept cannot certify independence of reused data or exclude selection and ascertainment effects.

### B2 — Smoking/adiposity explanations remain unresolved

The positive set includes three tobacco-use-disorder relationships and three COPD relationships, alongside gastrointestinal and musculoskeletal outcomes. No joint smoking/adiposity conditional model was executed for this family. The historical insomnia–BMI pilot confirms a known input relationship computationally; it does not quantify mediation of any extension relationship.

Do not describe any extension association as independent of smoking, independent of adiposity, a direct sleep effect, or a proportion mediated. A valid residual genetic-covariance model could ask a narrower statistical question, provided exact smoking and adiposity constructs, a complete joint S/V system, stable residual variances, uncertainty propagation and independent validation are admitted. Tobacco use disorder, smoking dependency, initiation and heaviness are not interchangeable covariates. Attenuation after conditional covariance is not causal mediation; covariate-adjusted GWAS are a different estimand and require their own selection assumptions. Rank B is a feasibility lead, not an admitted mechanism.

### B3 — The observed breadth does not establish eight independent biological processes

The eight external positive labels occur 4/5/1/2/3/3/4/1 times for abdominal pain/reflux/diaphragmatic hernia/noninfectious gastroenteritis/tobacco use disorder/COPD/disc disorders/cholelithiasis. Repeated labels, clinically correlated outcomes, correlated sleep traits, shared source samples and selection can produce similar association profiles. The 28 similarities cannot identify independent disease components, a common pathway, an effective number of diseases or a clinical syndrome. No outcome-outcome genetic/sampling covariance model was fitted for these eight outcomes.

The archived 42-trait odd/even-chromosome factor exercise failed its held-out residual criteria; it is not a validated latent biological phenotype. Its chromosome split assesses stability within the same study sources, not independent population replication. Do not import its exploratory factors, smooth matrices or discard indicators after outcomes to manufacture a disease-specific component. The separate statistical review's rg uncertainty export issue also prevents casual reuse of the archived fixed-diagonal uncertainty for new contrasts.

### B4 — Objective-versus-self-report conclusions are untested

The extension has 474/800 questionnaire, 61/300 device and 68/100 registry BH positives, with zero device traits among the 23 historical external positives. These counts do not test a genetic difference or prove that subjective sleep is more clinically relevant. Estimator precision, source availability, sampling, heritability, trait selection and measurement definitions vary. Within the device source, efficiency and duration share an estimated-sleep component; midpoint uses window boundaries. They are not three independent biological replications.

The device study used UKB participants, up to seven days of recording 2.8–9.7 years after baseline, and inactivity-based estimates; it acknowledges selection and difficulty distinguishing quiet wake from sleep. These source facts undermine a pure-modality interpretation and do not invalidate every device association. [Jones 2019 primary methods and limitations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6451011/). The present package has no covariance-aware clinically selected disease contrast or eligible external device-GWAS validation. A valid contrast would still concern the exact measurements, sampling and time windows; it would not isolate modality as a biological intervention.

### B5 — Clinical constructs and generalization remain limited

The original insomnia phenotype is frequent UKB complaints, not a diagnosis requiring chronicity, adequate sleep opportunity and daytime impairment. Ordinal sleepiness and napping do not supply hypersomnia, narcolepsy or severity; respondent-reported snoring is not diagnosed apnea. Duration tails against shared 7–8-hour controls do not establish a clinically optimal interval or a causal U-shaped disease curve. Signed interpretation differs across increased morning preference, later midpoint, longer duration, greater efficiency and symptom liability.

FinnGen R9 G6_SLEEPAPNO is the exact coded endpoint, not PSG severity, treatment response or certified exclusively obstructive apnea. The sleep-phenotyping reviewer documents an unresolved Finnish ICD-9 3472A narcolepsy/cataplexy label in the endpoint's official code table. This prevents certifying pure OSA; it does not establish an erroneous fraction of cases, since code counts overlap and coding history is unresolved. No case membership should be repaired from aggregate counts.

MVP PheCode 327.4 is a clinically recorded diagnosis construct related to frequent complaints. It is not an exact replacement. The [authors' primary resource](https://www.med.upenn.edu/syspharmatt/assets/user-content/documents/science.adj1182.pdf) acknowledges automated clinical-code phenotyping limitations. A future positive MVP×FinnGen test could support cohort-distinct transport across related definitions, with exact individual intersection unknown. It could still reflect healthcare utilization, coding, symptom burden or source composition. A negative test without suitable precision could reflect power or construct differences. Neither result establishes clinical equivalence, sleep treatment benefit or symptom-to-disorder progression.

### B6 — There is no locus-to-gene or clinical-actionability bridge

Source-ready tissue resources, annotations, workflow code, a nearest-gene column, and prior sleep-GWAS genes do not constitute a verified disease-pair molecular result. All 23 ranked local/molecular fields are blocked; the unified protocol admits no locus-level primary signal. There are no verified credible sets, matched signed LD analyses, molecular-trait colocalizations or multipronged effector-gene evidence for these relationships. Do not import a locus or pathway from a source paper as a new sleep–disease mechanism. Even an eventual colocalization would require its model/QC and would not by itself establish an effector gene, mediation or causality.

No eligible patient-level or longitudinal clinical analysis, predictive model, calibration, discrimination, absolute-risk estimate, severity endpoint, intervention outcome or treatment-response analysis was executed here. Therefore clinical risk stratification, screening, prevention and therapeutic targeting are untested. The frozen |rg|=0.15 power-planning target is a statistical design target, not a validated clinically important effect. Liability h2 conventions and nonphysical h2 values identified by the statistical reviewer cannot be presented as literal disease-risk fractions. A mechanism figure would currently lack traceable evidence.

## Human specialist questions and interpretation stops

| Reviewer | Exact question | Required disposition until resolved |
|---|---|---|
| Sleep clinician/phenotyping expert | Do field 1200 complaints and MVP PheCode 327.4 identify overlapping symptom/severity domains under their chronicity, opportunity, impairment and exclusion rules? | Keep transport wording; reject exact clinical replication |
| Finnish coding expert and sleep clinician | What is the code-history meaning of R9 3472/3472A, its overlap with G47.3, and the obstructive/central subtype coverage? | Retain coded sleep-apnea endpoint; stop pure-OSA, severity and treatment claims |
| Clinical epidemiologist | For each of the 13 historical outcome sources, what full code maps, observation time, encounter thresholds and case/control exclusions make a concept match defensible? | Qualify `EXACT` as curated concept-level match; stop identical-ascertainment claims |
| Statistical geneticist | Can exact smoking/adiposity constructs be jointly modeled with each sleep/disease pair using calibrated covariance and stable residual variances? | No independent-of-smoking/BMI or mediation claims |
| Clinical epidemiologist/statistical geneticist | How do UKB participation, device participation and VA/Finnish healthcare selection limit generalization? Are demographics, observation opportunity and comorbidities comparable enough for the frozen transport estimand? | No unconditional generalization across populations or clinical severity |
| Statistical geneticist/sleep researcher | Is a commensurate continuous-duration contrast identifiable with shared-estimator covariance and adequate power, and which source differences remain inseparable? | No modality superiority or distinct-mechanism claim from significance counts |
| Statistical geneticist/clinical disease expert | Are repeated outcome indicators clinically coherent and separable by a prespecified, stable, independently tested model? | No independent disease count, pathway cluster or validated latent phenotype |
| Molecular geneticist | Is any newly admitted shared locus supported by dense statistics, appropriate signed LD, credible-set QC, relevant molecular tissue and independent gene-prioritization evidence? | No gene/tissue/pathway or shared-causal-variant claim; no mechanism figure |
| Causal-inference specialist | What eligible design could separate horizontal/vertical sharing, source selection and reverse direction, and is it actually executed? | No causal direction, mediated fraction or intervention benefit |
| Clinical investigator | Which eligible clinical cohort and prespecified endpoint quantify prediction or treatment utility beyond established clinical factors? | Clinical utility remains untested; rg magnitude alone does not authorize action |
| Investigators/literature specialist | Does a source-eligible new result materially change prior conclusions after phenotype/source equivalence and recent competitors are checked? | No first-ever or novel mechanism claim; retain current NO-GO |

## Completion judgment

No important positive global estimate was numerically disproved by this biological review. Stronger interpretations fail because the required identifying design, independent sleep evidence or molecular/clinical bridge is absent or unadmitted. This distinguishes unsupported mechanisms from demonstrated absence of biology.

Recovery and reproduction support reliability of a numerical resource; they do not establish a distinct, clinically meaningful SLEEP advance. Historical 396- and 1,200-test families may be integrated descriptively with their separate families and source restrictions. The current evidence is not cleared for human authoring around a new mechanism, causal finding or clinically actionable discovery. A qualified resource handoff still requires the exact human reviews above and the parent's remaining native/source checks. Any future result must receive an updated review; this report does not prejudge the pending transport experiment.
