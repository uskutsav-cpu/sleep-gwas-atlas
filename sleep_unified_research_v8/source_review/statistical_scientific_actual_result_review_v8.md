# Independent review of the actual FinnGen insomnia–MVP abdominal-pain result

Reviewed 11 October 2026 against the separately frozen `protocols/exploratory_registry_insomnia_abdominal_pain.md`. Scope: the one completed new native fit and its six small output files. No estimator was rerun, no historical computation repeated, and no between-source contrast or between-job block covariance was estimated.

## Decision and scientific result

**Retain as qualified exploratory evidence of positive source/construct transport.** Under the frozen source and LDSC assumptions, FinnGen R13 registry-coded insomnia and MVP GIA European abdominal pain have an estimated common-variant genetic correlation of **0.39912075217856924**, SE **0.050695938374874**, nominal 95% interval **0.2997585388013542–0.49848296555578436**, Z **7.872834885257436**, two-sided nominal P **3.4669408448782324×10⁻¹⁵**. The allele-validated native intersection contains **815,347 SNPs**.

The interval width is 0.19872442675443014. This result provides considerably more information than a favorable sign alone: sampling uncertainty under the specified model supports positive sharing and excludes zero. It does not account for all uncertainty in source schema, ascertainment, LD reference or historical selection. Genetic correlation describes shared common-variant association structure; it does not show that insomnia causes abdominal pain or identify a molecular mechanism.

This is genuinely new **project evidence**: the FinnGen sleep-side source had not previously been paired with this MVP outcome. The result extends a selected symptom/complaint relationship to registry-coded insomnia under this qualified analysis. MVP abdominal-pain data are reused from the historical analysis, and those historical results informed selection. Therefore it is not fully independent two-trait replication, unbiased replication of a selected effect, a new field-level discovery, or an independently replicated outcome. Cohort, recruitment and phenotype ascertainment change together, so a pure clinical-versus-questionnaire measurement effect is not identified.

## Independent statistical reproduction from actual outputs

The companion JSON records 80-digit Decimal arithmetic from the actual captured totals and all three exported 200-row deletion files. The ratio was independently reconstructed as `gencov/sqrt(hsq1*hsq2)`. For each native deletion, the reconstructed ratio was `gencov_b/sqrt(hsq1_b*hsq2_b)`; SE was computed as `sqrt((B−1)/B * sum((ratio_b−mean(ratio_b))²))`. This is equivalent to the stock ratio-jackknife pseudovalue SE. The pinned `RG` and `RatioJackknife` implementations were inspected to distinguish the reported full-data ratio from its bias-corrected jackknife estimate.

| Quantity | Independently reconstructed | Native agreement |
| --- | ---: | --- |
| Full-data rg | 0.3991207521785693 | Difference 5.55×10⁻¹⁷ |
| Centered 200-delete ratio SE | 0.05069593837487401 | Difference 6.94×10⁻¹⁸ |
| Bias-corrected jackknife rg | 0.39826198290040615 | Difference 7.16×10⁻¹⁵ |
| Z | 7.872834885257436 | Same displayed binary64 value |
| Two-sided normal P via independent `erfc` | 3.466940844878248×10⁻¹⁵ | Relative difference approximately 4.55×10⁻¹⁵ |
| Nominal 95% interval | 0.2997585388013542–0.49848296555578436 | Uses full-data ratio and stock SE |

The publication estimate should use the stock reported full-data ratio 0.39912075217856924, not silently substitute the bias-corrected value. The small difference between these two estimates is expected from ratio jackknifing. No alternative estimator was selected after seeing the result.

All 200 deletion ratios were finite, between 0.3906100120867539 and 0.4102742736158114, with mean 0.3991250676020776. Every underlying h² and gencov deletion was positive and finite. The largest absolute change from the full-data rg was 0.011153521437242155; the largest individual squared centered deviation contributed 4.812432256612888% of their total. These descriptive checks reveal no obvious near-zero-denominator, sign-switching or single-deletion instability. They are not a newly invented acceptance threshold, leave-one-chromosome sensitivity, proof that no locus matters, or protection against systematic bias. The three deletion vectors belong to this native pair; their indices were not aligned to other jobs.

The separately reconstructed h²1, h²2 and gencov centered deletion SEs also match their actual native values within relative 10⁻¹²/absolute 10⁻¹⁵. These are bounded checks of the new output, not renewed million-row source audits.

## Pair-specific diagnostics and their limits

| Actual pair quantity | Estimate | SE | Nominal 95% interval |
| --- | ---: | ---: | --- |
| FinnGen observed h² under assumed N | 0.05202247810759592 | 0.004556401352741753 | 0.043092095557112506–0.06095286065807934 |
| MVP observed h² under assumed N | 0.05494454705381611 | 0.0034988122034778028 | 0.04808700114633039–0.061802092961301826 |
| FinnGen h² intercept | 1.084990267297004 | 0.010948801293037238 | 1.0635310110887655–1.1064495235052425 |
| MVP h² intercept | 1.1316006496545958 | 0.011207934743943389 | 1.1096335012153917–1.1535677980938 |
| Genetic covariance under assumed N | 0.021338413475760064 | 0.002789676715447691 | 0.0158707475849726–0.02680607936654753 |
| Genetic-covariance intercept | 0.014914689572365932 | 0.00786843457320914 | −0.0005071588058337684–0.030336537950565633 |

Pair h² Z values are 11.41745 and 15.70377, comfortably positive; both intercepts are below the historical necessary 1.2 bound. Passing these conditions does not certify source validity. Both h² intercept intervals exceed one. From actual outputs, independently reproduced `(intercept−1)/(meanχ²−1)` values are 0.28850609404716093 (FinnGen) and 0.28400515070860866 (MVP). These describe intercept-associated inflation under the fitted model, not fractions of causal confounding or a proof that the remaining signal is unbiased.

The estimated cross-trait intercept was unconstrained as frozen. Its interval includes zero and small positive values; it cannot establish absence of participant overlap, shared ascertainment, correlated bias or relatives. The model allows an overlap-associated intercept contribution, but a null-compatible estimate is not a cohort-disjointness certificate. [Official LDSC overlap/intercept guidance](https://github.com/bulik/ldsc/wiki/FAQ), [original cross-trait LDSC methods](https://www.nature.com/articles/ng.3406).

Pair-specific h² need not equal the earlier standalone h² because the pair uses a smaller intersection. These modest changes are not evidence of a failed reproduction or statistically significant source heterogeneity. No shared-outcome cross-source covariance is calibrated here; comparison of rg point estimates or confidence intervals cannot supply the missing formal difference test.

## Assumptions that remain attached to the result

The exact source contract and its evidence are documented in `source_admission_review_v8.md` and `qualified_native_input_contract.json`; this result does not retroactively clear the original V4 confirmatory requirements.

- FinnGen effective N is an assumed constant 185146.703773, derived from endpoint counts 51,643/446,273; no SNP-wise N or endpoint INFO column was verified. The existing HM3 fallback is an explicit source route, rather than a fabricated INFO>0.9 claim. The official LDSC workflow permits HM3 fallback and assumed constant N, with limitations. [LDSC source guidance](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).
- MVP effective N is the existing constant 274973.347782. Its Z was reconstructed from OR and exported CIs under the assumption of a 95% approximately log-Wald interval because supplied SE is unavailable in the bounded inspected sample. CI confidence level and exact construction remain unverified. The supplied SAIGE P need not equal Wald P; the earlier bounded disagreement diagnostic does not establish an error and was not used to alter inputs.
- The common European LD reference is an approximation for the Finnish source. Good numerical precision, positive deletion values and clean logs cannot eliminate systematic LD mismatch.
- Registry insomnia and MVP abdominal-pain coding do not certify identical chronicity, observation horizon, clinical adjudication or causal disease entities. Source and ascertainment variation are part of the transport question.
- Historical MVP outcome reuse and selection remain. All uncertainty quoted above is marginal model-based uncertainty conditional on the adopted input specification; it does not include a correction for the research selection process.

The nominal P is smaller than the historical 0.05/217=0.0002304147465437788 sensitivity cutoff. That observation demonstrates strength under a stricter descriptive threshold. It is not a newly executed 217-pair confirmatory family, a complete correction for source/method selection, removal of winner's curse or proof of independence.

## Native completion and bounded integrity evidence

The actual receipt is `native_transport/receipts/one_pair_execution.json` (the initially supplied `terminalreceipts` location did not exist). It records `WORKER_COMPLETE_VERIFIED`, return code 0, null stop reason, unchanged plan, completed process-group teardown and peak observed worker RSS 763,084,800 bytes. The native log reports the 815,347-SNP valid-allele intersection and successful completion. Neither the native log nor worker stdout contains warning/error/NaN/Inf/out-of-bounds lines. Six small current result hashes match this execution receipt. These facts support numerical completion, not biological validity.

Primary actual-result evidence:

- `native_transport/fits/finn_registry_insomnia_MVP_abdominal_pain.full_precision.json`
- `native_transport/fits/finn_registry_insomnia_MVP_abdominal_pain.log`
- The same output prefix followed by `insomnia.sumstats.gz_gwas_catalog_GCST90479148.sumstats.gz.hsq1.delete`, `.hsq2.delete` and `.gencov.delete`
- `native_transport/logs/one_pair_execution.stdout.log`
- Companion `source_review/statistical_scientific_actual_result_review_v8.json`, containing every independent calculation and exact paths.

## Manuscript use and next work

Recommended sentence: “In a separately specified exploratory source-transport analysis, FinnGen registry-coded insomnia showed positive common-variant genetic correlation with MVP abdominal pain (rg=0.399, SE=0.051, nominal 95% CI 0.300–0.498; P=3.47×10⁻¹⁵; 815,347 allele-aligned SNPs). This analysis used a new sleep source but reused the historically selected MVP outcome and retained source-specific effective-N, CI-reconstruction and LD-reference qualifications.”

This result strengthens the resource's demonstrated use case and supports a qualified figure panel showing clinical-source transport. It does not close mechanism, locus-level colocalization, independent two-trait validation or disease-specificity gaps, and cannot guarantee Cell Genomics acceptance. The related clinical-versus-habitual-sleep multi-biobank preprint means that broad clinical-sleep framing itself is not a new contribution. The defensible contribution must be established against actual resource prior art and the integrated atlas's added harmonization, uncertainty and reusable output. [Relevant clinical-sleep prior art](https://doi.org/10.1101/2025.11.21.25340447).

Essential next work is to incorporate this one result and qualifications into the release table, research figure, source-aware methods and scientific interpretation; then independent scientific reviewers should assess the integrated evidence and manuscript claims. No repeat fit, controller version, broad unchanged-source audit, source retry or calibrated cross-source difference claim is justified by this clean completed result.
