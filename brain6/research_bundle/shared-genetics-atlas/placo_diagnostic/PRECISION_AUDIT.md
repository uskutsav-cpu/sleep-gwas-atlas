# Source-statistic precision check

The same deterministic hash-selected 2,191 markers used in the control diagnostic were matched by rsID and aligned allele pair to the filtered R13 insomnia source, the raw Atkins FI file, and the Ye FFS source. All 2,191 matched in each source; one duplicate FI rsID was encountered after its first record and recorded. The FFS A1 interpretation is the project's previously inferred convention, not newly author-confirmed. The complete matched rows and reproducible script are in `precision_audit_sample.tsv` and `precision_audit.py`; full input source paths and audit counts are in the script and `precision_audit_summary.json`.

For all three sources, the signed Z in the scan input matched the normal two-sided Z derived from the **reported p-value and beta sign** to within 0.0005 in every sampled row. This validates the rounded, p-derived representation for these sampled markers; it does not verify the other 1.136 million or establish which of a source's statistics is best for cross-trait inference.

| Source | Median absolute difference: scan Z vs beta/SE | 95th percentile | Count with difference > 0.01 of 2,191 |
|---|---:|---:|---:|
| R13 insomnia | 0.00025 | 0.00048 | 0 |
| FI | 0.00640 | 0.0220 | 659 |
| FFS | 0.02157 | 0.0632 | 1,644 |

These differences can reflect source precision, test statistics or model conventions, among other possibilities; this check alone does not determine the cause or prove that either representation is wrong. For example, at rs17124028 the FFS p-derived aligned Z was +1.227 while aligned beta/SE was +1.340. Near a decision threshold, using p-derived versus Wald Z may change a cross-trait p-value or ordering. FI beta and SE in the source are printed to limited precision, and the FFS source uses BOLT-LMM; exact statistical definitions should be checked against source documentation before replacing either input. Run and compare a full-precision, source-justified sensitivity analysis, especially for candidate hits. The scan's approximate results are not publishable confirmations until that discrepancy, method calibration, coverage and independent validation are resolved.

The [BOLT-LMM manual](https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/BOLT-LMM_manual.html) distinguishes beta and SE from an infinitesimal-model approximation and an association p-value for its non-infinitesimal model. Ye et al. used BOLT-LMM, so distinct model outputs are a plausible explanation for FFS discrepancies. The public FUMA-style file labels its field only `P`; this review has not verified the exact upstream column mapping or BOLT configuration, and cannot say that this is the definitive explanation for the observed differences. Do not replace the reported p-values with Wald p-values solely because they differ. [Ye et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/)
