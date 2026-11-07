# UK Biobank and FinnGen hospital frailty: insomnia-source comparison

A subsequent [dementia-weight exclusion sensitivity](dementia_sensitivity/RESULTS.md) finds similar genome-wide insomnia correlations for original and modified UKB HFRS on the same markers. This is not a sleep-item-excluded FI or FFS.

The authors' [public HFRS release](https://doi.org/10.6084/m9.figshare.31259023.v1) includes full UK Biobank and FinnGen summary statistics from [Mak et al. 2025](https://www.nature.com/articles/s43587-025-00925-y). Using the newly available UK Biobank file, we compared both HFRS sources with questionnaire and registry insomnia on **1,020,621 common markers**. This is exploratory and uses the original cohorts, not independent two-trait replication.

| HFRS source | Questionnaire insomnia r_g (SE) | Registry insomnia r_g (SE) |
| --- | ---: | ---: |
| UK Biobank | 0.5288 (0.0334) | 0.3755 (0.0416) |
| FinnGen | 0.4363 (0.0253) | 0.6132 (0.0326) |

All four individual correlations meet Bonferroni correction over these four tests; exact p-values are in [results.json](results.json). For a separate family of four paired comparisons, UKB minus FinnGen HFRS gives **+0.0925** (paired SE 0.0355, adjusted p=0.0366) with questionnaire insomnia and **−0.2377** (SE 0.0426, adjusted p=9.87×10⁻⁸) with registry insomnia. Questionnaire minus registry insomnia gives **+0.1533** (SE 0.0477, adjusted p=0.00525) with UKB HFRS and **−0.1769** (SE 0.0397, adjusted p=3.31×10⁻⁵) with FinnGen HFRS. Contrasts use the rounded r_g values printed by LDSC and a 200-block paired jackknife. These small families were planned after prior analyses, so the correction is not project-wide.

The pattern is stronger when each HFRS source is paired with insomnia from the same named cohort. **This does not establish that participant overlap caused the pattern.** The HFRS sources also differ in hospital ascertainment, ancestry, phenotype distributions, GWAS models and reference fit. Exact participant overlap and relatedness have not been audited. FI includes sleep and tiredness, FFS includes exhaustion, and neither was modified in this analysis. HFRS is a separate diagnosis-based score. These genome-wide results do not demonstrate causation, a shared locus, or rescue the failed local-correlation gates.

The UKB HFRS–questionnaire insomnia cross-trait LDSC intercept is 0.0672 (SE 0.0062); FinnGen HFRS–registry insomnia is 0.1795 (0.0078). The two cross-cohort intercepts are near zero in these fits. An intercept alone is neither a count of shared participants nor proof of overlap or its absence. Finnish HFRS has an elevated univariate intercept **1.1866 (0.0139)** with the general-European LD reference; ancestry mismatch, sample-size conventions, ascertainment and model fit remain unresolved.

## Important source discrepancy

The authors' article and release description state **407,463** UKB participants, while **all 9,812,854 rows** in the released UKB summary file report N=**407,643**. The reason is unknown. We used the per-row file value in this exploratory run and did not silently change it; this discrepancy prevents treating the source's sample metadata as verified. Registry insomnia and Finnish HFRS also use study-level N rather than independently verified per-variant effective N. Results may be affected by model and source assumptions beyond this numerical discrepancy.

## Reproducibility and quality checks

The complete 280,558,723-byte UKB gzip was streamed without retaining the full archive, verified against its publisher MD5 and read to gzip completion. Source SHA-256 is `975a6c9b0e8db48326b5c02d711c2f9c16180ff9595eef52d69ec8a73c4b5aff`. Among the previously selected 1,020,637 markers, 1,020,631 matched the source by rsID, GRCh37 reference position and allele pair; 10 failed the specified INFO≥0.9 or MAF≥0.01 rule, leaving 1,020,621. No matching IDs were duplicated. For all 1,020,621 retained rows, normal p-values computed from BETA/SE were within 0.02 on the log10 scale of the source p. This checks numeric consistency, not GWAS calibration or source-level allele semantics. A1 follows the conventional [fastGWA effect-allele definition](https://yanglab.westlake.edu.cn/software/gcta/); the entire source's build and processing were not independently authenticated beyond this reference crosswalk.

All four LDSC runs retained the same markers; shared-trait heritability delete vectors matched exactly across each paired contrast, and reconstructed jackknife SEs matched logs to printed precision. Files, selection thresholds, checksums, N counts and code are retained in the [protocol](PROTOCOL.md), [stream receipt](stream_receipt.json), [matching receipt](prepare_receipt.json), [source streaming code](stream.py), [matching code](prepare.py), [run script](run.sh), [paired calculation](summarize.py) and complete per-run logs. Filtered inputs remain on the SSD. The 16-marker difference from the earlier [five-outcome comparison](../../component_exclusion/physical_function_ldsc/insomnia_source_sensitivity/frailty_extension/RESULTS.md) means exact numerical results should be read with their own support and correction family.
