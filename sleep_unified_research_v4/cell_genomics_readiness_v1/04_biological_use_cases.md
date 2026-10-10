# Demonstrated biological use cases

The atlas supports source-aware retrieval and prioritization of published GWAS genetic-sharing estimates. It does not identify causal effects, clinical equivalence, molecular mechanisms or differences between correlated estimates.

1. **Retrieve a well-estimated relationship with its actual phenotype definition.** The core table contains insomnia–frailty at r_g = 0.640506 (SE 0.023535; marginal 95% CI 0.594378 to 0.686634). The paired construct record specifies frequent UK Biobank sleep-initiation/maintenance complaints, not a clinical insomnia diagnosis. A researcher can locate the result and inspect the source before interpreting its scope.
2. **Find a disease-sharing signal in a broad, separately corrected screen.** The extension table includes insomnia–abdominal pain at r_g = 0.559063 (SE 0.035898). It belongs to the fixed 1,200-pair family; it is not pooled with the 396 core tests. The result supports shared common-variant signal for the source-defined traits, not causation or mechanism.
3. **Distinguish a qualified external-outcome check from full replication.** The 217-candidate summary reports 23 qualified outcome-side positives, 18 directionally concordant below the frozen threshold, 17 QC-ineligible, and 159 unavailable. The 41 estimated rows reuse discovery sleep GWAS, so the table cannot establish independent two-trait replication.
4. **Inspect whether a result is sensitive to the fitted specification.** Device-derived sleep duration–HDL changes from r_g = −0.049705 (95% CI −0.095451 to −0.003959) at baseline to −0.028438 (−0.064299 to 0.007423) under the prespecified sensitivity. The latter marginal interval includes zero. The difference itself has no calibrated uncertainty or test; this is a useful flag for follow-up, not evidence that the effects differ.

These examples are selected to show both positive and limiting use. They reuse existing results and are not new analyses.
