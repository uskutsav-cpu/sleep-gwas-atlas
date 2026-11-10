# Exploratory genetic correlations among frailty measures

I ran two new LDSC comparisons with the project's existing processed sources: Atkins FI meta-analysis versus FinnGen R12 HFRS and Ye FFS versus the same HFRS. Both used the existing European reference LD scores, unconstrained intercepts and 200 jackknife blocks. They are exploratory analyses using pre-existing data, not independent remeasurement of the frailty phenotypes.

| Pair | Genetic correlation | SE | Approximate 95% CI | Regression SNPs |
|---|---:|---:|---:|---:|
| FI–HFRS | 0.640 | 0.0246 | 0.591–0.688 | 1,150,066 |
| FFS–HFRS | 0.553 | 0.0234 | 0.507–0.599 | 1,143,326 |

Both point estimates are positive. FI has the larger point estimate, but **the two results use different SNP support and share the same HFRS data**. Those separate runs do not themselves supply a covariance-aware comparison. A subsequent paired analysis on common SNP support is reported in `hfrs_fi_ffs_paired_contrast.md`; comparing these separate intervals or p-values alone does not answer that question. The project's earlier FI–FFS result around 0.75 concerns a different pair and matched supports; the FI–FFS correlation is not itself compared with either HFRS correlation here.

The HFRS single-trait intercept in the paired runs is ~1.189, and the HFRS input is Finnish while the reference LD is general European. Its per-variant sample sizes and INFO scores were not available, so preprocessing used the total N and could not screen INFO. This reference/model concern is material even though HFRS h²/SE exceeds the project's ≥4 precision threshold. A near-zero fitted cross-trait intercept cannot prove absence of bias or participant overlap. The UKB-based FI/FFS sources and FinnGen HFRS are distinct cohorts by description, but exact intersection has not been verified. HFRS is a diagnosis-weighted measure, not FI or FFS, and includes symptom categories related to drowsiness or exhaustion. None of the estimates proves a shared causal variant, disease mechanism or predictive risk.

## Reproducibility and checks

The original processed inputs were read unchanged. For both runs I independently recomputed rg as the saved genetic covariance divided by the square root of the two saved SNP-heritability estimates, and recomputed its SE from all 200 delete-block estimates; both matched the saved output to numerical precision. The 95% CIs use a normal approximation. `hfrs_fi_ffs_new_exploratory_ldsc.csv` includes exact values, intercepts, sample support and hashes; `FI_HFRS_exploratory.exact.json` and `FFS_HFRS_exploratory.exact.json` preserve the observed native fitted results and deletes. Full command logs are in the task's `work/` directory. The existing HFRS–insomnia sensitivity and intercept audit are documented in `hfrs_insomnia_ldsc_intercept_audit.md`. A suitably matched Finnish LD reference, review of the HFRS intercept, and independent replication are needed before treating the subsequent common-support paired contrast as conclusive.
