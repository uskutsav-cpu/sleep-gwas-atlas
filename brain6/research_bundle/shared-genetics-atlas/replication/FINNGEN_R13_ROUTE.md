# Candidate replication route: FinnGen R13 ADHD

The [FinnGen R13 endpoint F5_ADHD](https://r13.risteys.finngen.fi/endpoints/F5_ADHD) is called “Disturbance of activity and attention,” a registry-defined phenotype. Risteys displays 5,559 FinnGen individuals for that endpoint; this descriptive count is **not verified here as the GWAS case count or effective sample size**. The related endpoint F5_HYPERKIN is labeled as excluding ADHD in an earlier release and must not be substituted based on its name. Registry phenotypes can differ from the Psychiatric Genomics Consortium ADHD study's clinical cohorts and case definitions.

The [FinnGen access page](https://www.finngen.fi/en/access_results) states that R13 summary-statistic users should fill an online request form and receive download instructions by email; the browser is separately offered for exploring results. R13 effects and frequencies are reported for the GRCh38 alternate allele. We have existing R13 insomnia summary data, but did not locate an F5_ADHD file with verified authorization and provenance in the project paths reviewed. A public-browser search for rs62264764 did not produce a usable result during this check; this is an observation failure, not evidence of no association. No request form was submitted and no ADHD summary statistics were downloaded.

If the appropriately authorized F5_ADHD source becomes available, a useful **independent-cohort candidate consistency analysis** would:

1. Verify release, exact endpoint definition, cases, controls, source model and accession; document whether FinnGen participants overlap any PGC discovery ADHD cohort, and whether the proposed sleep cohort overlaps the discovery sleep GWAS. Unknown overlap prevents calling the result independent.
2. Freeze the previously identified Brain6 regions and their discovery-selected lead or credible-set rules before examining the new results. Validate hg19-to-hg38 mappings, reference and alternate alleles, multiallelic representation, frequency, imputation quality, and source-specific effect direction. Use the same pre-specified allele and region rule across traits.
3. Report each independent sleep and ADHD marginal beta, SE, p and sample size, and a conservative max-p test with correction for the number of frozen locus–pair tests. A null or weak result with this endpoint may reflect limited power or phenotype differences. If the two new FinnGen traits use overlapping participants, the max-p intersection–union test remains conservative when each marginal p is valid, but the two estimates are not independent of one another.
4. Match population-specific LD and run multi-signal colocalization only when variant coverage and QC are adequate. Replication of two marginal associations is not proof that the same causal variant drives both, and exact source-cohort independence must be established separately.

Neither an ADHD association nor a sleep association has been estimated in F5_ADHD for the Brain6 candidates in this check. This route does not change Brain6's existing `UNSUPPORTED` local-evidence class or rescue the failed LAVA gate.

A separate [illustrative rs62264764 power calculation](rs62264764_power/RESULTS.md) shows how an ADHD replication can be underpowered under plausible assumed effect sizes and a six-test correction. The available endpoint count and hypothetical control count are not verified GWAS sample sizes.

A separate [childhood ADHD-symptom lookup](childhood_symptoms/RESULTS.md) found same-direction nominal rs77960 evidence in a symptom GWAS whose authors report no participant overlap with the diagnostic ADHD source. This is related-phenotype consistency on the ADHD side, not independent two-trait sleep–ADHD replication.

A [2026 family-based ADHD source assessment](familial_adhd_2026/RESULTS.md) found useful polygenic background but no full variant-level results for the Brain6 leads in the checked article and supplement. It does not establish independent two-trait replication.

The public GWAS browser phenotype table subsequently displayed **5,559 cases and 489,493 controls** for F5_ADHD; see the [updated power sensitivity](rs62264764_power/finngen_browser_counts/RESULTS.md). The variant association page did not load, and participant overlap remains unverified.

A subsequent [archived munged-file check](finngen_adhd_archived_munged/RESULTS.md) found a file labeled F5_ADHD, but could not verify its raw source or processing receipt. It contains only one of the six leads, with a weak derived association; it does not establish clinical or two-trait replication. The earlier statement that no verified F5_ADHD source was located remains correct.
