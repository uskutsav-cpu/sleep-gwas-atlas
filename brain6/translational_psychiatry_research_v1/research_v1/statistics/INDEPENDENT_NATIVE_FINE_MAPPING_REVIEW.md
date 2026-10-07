# Independent automated numerical review of native candidate C fine-mapping

The new native results support a **source/reference-qualified shared-signal model at the previously reported chr5 insomnia–ADHD region**. Separate numerical review found no inconsistency in the saved fits, credible-set coverage/purity, variant PIPs, ELBO histories or colocalization calculations. This is materially stronger methodological support than the historical single-signal ABF alone. It is not independent cohort replication, a newly discovered genomic region, a uniquely identified causal variant or a causal gene/mechanism.

This review reads existing native RDS and exact signed LD without refitting. It uses independently written R and Python calculations, not a second human statistical review. All newly read raw GWAS inputs, native per-SNP BFs and fitted RDS remain local under ignored `work/`; the public review contains only derived checks and source hashes.

## Source and method contracts checked

The fixed `FINEMAP_PRE_FIT_PROTOCOL.md` SHA256 matches its freeze receipt and the fine-mapping input manifest. The actual shared LD, two local trait-input hashes and complete chr5 LD source hash match their recorded receipts. The shared 2185×2185 signed matrix is exactly the allele/key-ordered submatrix of the full2694-variant source LD: maximum difference **0**. The reference manifests contain **503** EUR samples, with503 recorded at every reference variant. This corroborates matrix ordering/identity; it is not an independent replay of raw genotype parsing or a guarantee of source-matched LD.

The native fitted methods are pinned susieR0.14.2/coloc5.2.3. Four fixed model conditions produce **8 trait fits**: L10 with median effective N, L1 and L5 with the same N, and L10 with n=∞. The **12 posterior conditions** are4 signal-pair comparisons ×3 p12 settings, not12 additional fine-mapping fits.

The original source-N contract remains an approximation. ADHD effective N is constructed from per-variant case/control counts; insomnia effective N uses literal per-variant analyzed N and the published overall case fraction. Variant-specific case fraction is unavailable. The source data are binary association statistics and consortium/meta-analysis contributions; a single median effective N does not recover a common-sample Gaussian regression likelihood. The n=∞ sensitivity is also an approximation. Similarity across these arms checks model sensitivity; it does not validate exact N semantics.

The current candidate C inputs supply frequency and INFO for both traits. Harmonization requires finite fields, INFO≥0.9, source MAF≥0.01, and an absolute source/reference ALT-frequency difference≤0.15. Insomnia uses the source MAF field as A1 frequency; ADHD uses FRQ_U_186843, the control frequency. Missing frequency/INFO documentation in some historical inputs does not apply to this current experiment. Operational `estimate_s_rss` estimates are about0.01409 for insomnia and0.00917 for ADHD, with zero markers meeting the frozen combined allele-switch rule. These are useful compatibility diagnostics, not full empirical calibration or proof of correct allele identity. The source-frequency interpretations, binary effective-N approximation, finite external503-person LD, and source ancestry/cohort heterogeneity remain limitations.

## Saved-fit probability and convergence checks

All8 saved fits converge and have finite alpha probabilities in[0,1], with rows summing to1 within5.56e-15. A separate stable calculation of `PIP=1-exp(sum(log(1-alpha)))` reproduces the full exported variant PIPs within **5.60e-16** maximum absolute error. All probabilities are in[0,1]. Overall PIPs are not within-signal posterior probabilities and are not summed to assess a particular credible set.

The saved ELBO histories match their exported tables within **4.95e-10** absolute error, consistent with text roundtrip rounding. The minimum ELBO step across all arms is **0**; none drops beyond the frozen1e-6 tolerance. Convergence after2–6 iterations is the software's variational stopping criterion, not proof of a global optimum, exact likelihood specification or posterior calibration. No additional refinement/multistart fits were run by this reviewer; the fixed L sensitivity arms are retained as executed.

Each trait/model has **one qualifying credible set**. The L10 model does not establish10 causal effects, and two credible sets in the primary analysis means one per trait, not two independently resolved shared signals.

| Primary fit | Credible-set size | Independently summed within-signal coverage | Minimum absolute LD purity | Maximum variant PIP |
|---|---:|---:|---:|---:|
| Insomnia | 70 | 0.95112159 | 0.64644916 | 0.13730068 |
| ADHD | 8 | 0.95828605 | 0.72791690 | 0.29472186 |

Every saved credible-set membership is consistent with the minimal highest-alpha set crossing95% within-signal probability. Exact signed LD yields matching minimum/mean/median absolute purity, with maximum discrepancy **6.78e-15**. The smallest conditional coverage across all8 arms is **0.95000777**; minimum purity is **0.64644916**. Coverage is posterior mass under the fitted model, not a verified95% frequentist coverage guarantee.

All8 primary ADHD credible-set markers occur in the70-marker insomnia set. Both historical lead rs2431108 and published rs77960 occur in both primary sets. Their strong reference LD and similar model support prevent a distinct-causal-variant claim from their different names. No single variant approaches posterior certainty; this region remains unresolved at variant level.

## Independent colocalization checks

All12 H0–H4 vectors are finite, bounded and sum to1. An independently implemented bounded **explicit i≠j BF configuration sum** reproduces every vector from the saved native signal BFs within **3.45e-15** maximum absolute error. Independently reweighting only H4 by the p12 ratio and renormalizing reproduces the prior-sensitive posteriors within **1.89e-15**. The reported `H4/(H3+H4)` values also match their direct ratios.

Primary H4 is **0.99156823** at p12=1e-5 and **0.92162951** at1e-6. All fixed L/N arms remain above0.8 at the lower prior. This supports sharing **conditional on the source/reference, fitted signal model and priors**. High H4 is neither a frequentist family-significant P value nor proof that one specific SNP causally mediates both traits. The cohort sources are the original discovery insomnia and ADHD studies; independence count is **zero**.

The remaining interpretations to withhold are:

* A first discovery of the chr5 region: Zu2026 already reported rs77960 and strong single-signal coloc there.
* A new distinct causal signal: the same prior lead lies inside both credible sets and tags the historical lead strongly.
* Multiple distinct shared effects: each trait has only one qualifying credible set in every tested arm; weak/outside-window effects cannot be excluded.
* A causal gene, cell type, molecular mechanism or therapeutic target: none follows from trait–trait H4; the larger historical functional window supplies no newly admitted strong-QTL candidate C evidence.
* Independent replication, clinical prediction or symptom relevance: no new psychiatric cohort was used.
* Restoration of the frozen LAVA family or previously failed32LD gate: this is a distinct PSD-compatible RSS experiment, not a change to those historical decisions.

## Artifacts and reproducibility

`native_fine_mapping_RDS_checks.tsv` contains8 saved-fit/CS reviews. `native_coloc_configuration_and_prior_checks.tsv` contains12 independent posterior checks. `native_fine_mapping_review_summary.json` summarizes the results, and `native_fine_mapping_review_input_manifest.json` records exact authoritative input hashes with local-only restrictions. No full fit was repeated.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/homebrew/bin/Rscript brain6/translational_psychiatry_research_v1/research_v1/statistics/review_native_fine_mapping.R /absolute/path/to/sleep-gwas-atlas
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/review_native_fine_mapping.py
```

These local integration commands require actual fitted objects, source inputs and LD. Missing inputs fail; they are not replaced with simulated data or silently skipped. The final source-free checks assess only committed derived evidence, and actual local-data integration tests are separately recorded. Independent human statistical-genetics review remains required before scientific promotion beyond this known-region model evidence.
