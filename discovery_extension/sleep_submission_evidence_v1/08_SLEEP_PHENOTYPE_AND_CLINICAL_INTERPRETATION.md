# Sleep phenotyping and clinical interpretation audit

The 12-trait locked panel contains seven self-report questionnaire measures, one participant questionnaire item about a partner/close-contact complaint of snoring, three wrist-actigraphy measures and one hospital-registry sleep-apnea phenotype. All exact definitions, source releases and sample sizes are in `tables/sleep_trait_metadata.tsv`.

Chronotype is a self-reported morning-versus-evening preference, not objectively recorded circadian phase. Insomnia is frequent trouble falling asleep or waking during the night in UKB field 1200, not clinical insomnia diagnosis. Short sleep is <7 hours and long sleep ≥9 hours against 7–8-hour controls. Duration tails share controls and are not independent sleep mechanisms. Daytime sleepiness and napping are questionnaire phenotypes. Snoring is a participant report of others' complaints, not confirmed sleep apnea.

Actigraphy measures are rank-normalized SPT-window sleep efficiency, duration and midpoint in approximately 84,810–85,449 UKB participants, with a single shared source archive. They estimate wearable-derived behavior and do not constitute polysomnography. They have substantially smaller N than the questionnaire GWAS; direct comparisons of significance conflate measurement and power.

FinnGen R9 sleep apnea is hospital-discharge/cause-of-death ICD coding with specified control exclusions. The locked panel reports 38,998 cases and 336,659 controls. Old `config/traits.tsv` values conflict and must not replace the locked metadata. Diagnosis registries may miss undiagnosed cases and do not provide apnea severity or treatment response.

Of the 23 qualifying validation pairs, 20 use self-report sleep traits, two use reported snoring and one uses registry sleep apnea. None uses one of the three actigraphy sleep traits. This concentrates support on self-reported behavior and symptoms. It is not evidence that objective sleep has no relationship: the validation family was selected, source availability is incomplete, and actigraphy has less power.

Effect conventions are retained per source: higher insomnia/short-sleep/long-sleep liability, morning-person liability, greater sleepiness/napping, greater snoring/apnea liability, higher efficiency, longer duration or later midpoint. A negative rg does not have a uniform clinical meaning across these constructs. Binary effective N is calculated as `4/(1/cases+1/controls)` and never substituted for original total N without labeling.

Genetic correlation describes covariance of common-variant effects under LDSC assumptions. It does not establish an individual's risk, causal direction, shared causal variant, treatment target or benefit of changing sleep. Questionnaire and registry associations do not provide clinical screening accuracy. Shared smoking/adiposity-related architecture remains a plausible confounder of interpretation; it is not quantified by an unsupported conditioning procedure.

Native sleep heritability replay is blocked by absent summary statistics and LD reference files. A human sleep researcher should check construct comparability and effect coding before manuscript interpretation.

Bibliographic repair in the new metadata namespace: five sleep-trait rows carry newly verified DOI fields for [Wang 2019 daytime sleepiness](https://pubmed.ncbi.nlm.nih.gov/31409809/), [Dashti 2021 napping](https://www.nature.com/articles/s41467-020-20585-3), and [Jones 2019 actigraphy](https://www.nature.com/articles/s41467-019-09576-1). Historical `UNRESOLVED` DOI strings remain unchanged in the original and `original_DOI`/`doi` columns. This is a source-citation correction, not a scientific-result change. Receipt: `sources/provenance_sleep_bibliography_corrections.json`.
