# Snoring liability-scale heritability source recheck

Date: 2026-09-26 (UTC)

## Source and prevalence choice

The locked sleep panel defines snoring as UK Biobank field 1210, the Campos et al. habitual-snoring GWAS (PMID 32060260; DOI [10.1038/s41467-020-14625-1](https://doi.org/10.1038/s41467-020-14625-1)). The registered source reports 408,317 participants, 152,302 cases, and 256,015 controls. The analyzed case fraction is 152,302/408,317 = 0.3729994. The frozen panel records K=0.37 and identifies this as the analyzed UK Biobank European sample prevalence. The paper's Table 3 footnote states that liability-scale estimates assume equal population and sample prevalence. This is a cohort-frequency assumption, not an independently measured population prevalence; that limitation should remain explicit.

The GWAS source, phenotype, and case/control counts match the current munged summary-statistics input. No alternative prevalence or phenotype was introduced. The source article reports a rounded liability-scale estimate near 9.9%; the rerun below is a project-specific LDSC result and is not represented as an exact reproduction of that published estimate.

## Reanalysis

The original observed-scale h2 log is preserved in repository history at `frailty_paper/results/frailty_v1/logs_sleep_panel/h2_snoring.log` in commit `8aaacd898859178273c82a7df45fb0fb6dadba92` (SHA-256 `eeaef8932fed8f120234fb2f9b03f5fd4af9deba7b67ed07bf082554a73f476d`). A separate LDSC run used the same munged input, EUR LD scores, and HapMap3-filtered SNP set as the original result, with `--samp-prev 0.3729994097723093 --pop-prev 0.37`. Its captured run log is `results/frailty_v1/logs_snoring_liability/h2_snoring.log`; the retained log has only trailing whitespace normalized. Checksums for the munged input, 44 chr1–22 LDSC reference files, prevalence config/lock, LDSC script, runtime, and output log are recorded in `analysis/snoring_liability_h2_run_manifest_2026-09-26.json`. The run timestamp is 2026-09-26 04:32:13 UTC (2026-09-25 23:32:13 CDT). The current panel log and downstream table now use the liability-scale estimate.

- Observed-scale historical result: h2=0.0647 (SE=0.0025), intercept=1.0309 (SE=0.0091).
- Liability-scale result under the frozen prevalence convention: h2=0.1054 (SE=0.0041), Z=25.71, intercept=1.0309 (SE=0.0091), LDSC QC PASS.
- No association statistics, sample, LD reference, pairwise rg, frozen analysis plan, LAVA inputs, LAVA lock, or QC thresholds changed. Since the scale conversion is multiplicative, the h2 Z statistic and intercept gate are unchanged.

The aggregate h2 master and Supplementary Table S4 were rebuilt from the current panel logs. The prior observed-scale estimate is retained for traceability, not as an alternate primary estimate. Exact participant overlap with frailty GWAS remains unknown.
