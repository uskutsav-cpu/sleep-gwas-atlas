# Audit of paired FI–physical-frailty genetic-correlation comparisons

The comparisons use **paired genomic-block jackknife uncertainty**, which accounts for covariance between the two estimated genetic correlations under the LDSC model and jackknife assumptions. Recomputing the standard errors from the saved delete-block estimates reproduced the reported individual SEs and the paired-difference SEs for all six comparisons below.

| Insomnia source | Frailty pair | Physical minus FI r_g | Paired SE | SE if covariance were set to zero, illustration |
|---|---|---:|---:|---:|
| FinnGen | Atkins / Ye | −0.185 | 0.040 | 0.053 |
| FinnGen | Mak FI / FP | −0.164 | 0.026 | 0.045 |
| MVP | Atkins / Ye | −0.155 | 0.036 | 0.053 |
| MVP | Mak FI / FP | −0.131 | 0.026 | 0.044 |
| UK Biobank | Atkins / Ye | −0.157 | 0.028 | 0.036 |
| UK Biobank | Mak FI / FP | −0.148 | 0.021 | 0.033 |

The final column is a calculation showing what happens if covariance is set to zero; it is not a recommended alternative analysis. The paired SE is smaller in these examples because the estimated jackknife covariance is positive. The covariance itself does not quantify participant overlap or eliminate bias due to sample overlap, phenotype definition, population structure or LD-reference mismatch. The estimates are post hoc, and their uncertainty is not adjusted for the wider set of analyses explored.

Each row compares correlations calculated on exactly the same ordered SNP support and 200 matching genomic blocks within its insomnia source. Support differs between insomnia sources. These values are from the four-frailty-source matched-support analyses; they can differ slightly from the original-pair table in [the main assessment](mentor_risk_assessment_updated.md), which uses a different intersection of source files. Neither set is a new independent sample.

The [machine-readable audit](insomnia_frailty_paired_uncertainty_audit.json) records each estimate, covariance-derived correlation and SHA-256 hashes of its exact LDSC inputs. The reproducible script is `work/audit_paired_ldsc_covariance.py`. Full source-specific assumptions and diagnostic outputs are in the [FinnGen](mak2023_insomnia_fi_fp_comparison.md), [MVP](mak2023_mvp_insomnia_fi_fp_comparison.md) and [UKB](mak2023_ukb_insomnia_fi_fp_comparison.md) reports. No participant-level observations were used.
