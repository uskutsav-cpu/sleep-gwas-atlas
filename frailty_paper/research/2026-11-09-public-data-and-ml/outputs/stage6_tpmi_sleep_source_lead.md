# TPMI sleep-disorder GWAS: source lead

The [Taiwan Precision Medicine Initiative PheWeb](https://pheweb.ibms.sinica.edu.tw/pheno/327) displays a sleep-disorders phenotype (phecode 327) with 45,836 cases and 287,785 controls, and links to SAIGE and PLINK2 summary-statistic downloads. Its category is neurological; it is a broader diagnostic grouping, not an insomnia-only phenotype or measured joint insomnia–frailty outcome. The phenotype list includes related subtypes such as sleep apnea and sleep-related movement disorders. Source genotype build, case-code composition, effect-allele convention, variant coverage, cohort overlap and permissions for analytical reuse must be verified before querying candidates or interpreting results. The site reports temporary download restrictions and requests prioritizing downloads one at a time; respect those limits.

This is a promising separate-population summary-level follow-up, not a confirmed independent replication. Individual genotype and phenotype data, FI/FFS outcomes, exact source/sample independence and score calibration have not been established. No summary file was downloaded or candidate association computed in this screen. The [about page](https://pheweb.ibms.sinica.edu.tw/about) describes the Taiwan Precision Medicine Initiative. The search page's fuzzy matching should not be treated as an exhaustive absence test for other phenotypes.

The source [About page](https://pheweb.ibms.sinica.edu.tw/about) specifies GRCh38 and asks users of the data to acknowledge TPMI and cite Chen et al. (2025), DOI 10.1038/s41586-025-09350-y, and Yang et al. (2025), DOI 10.1038/s41586-025-09680-x. The browser-mediated SAIGE file request did not yield an accessible completed file within 60 seconds. A local HTTPS request encountered certificate verification failure; no certificate checks were disabled. This outcome does not establish that the data are unavailable, and no candidate-level results or score test can be reported from it. The site asks users to download one file at a time and retry after a few days if limited; do not assume a timeout proves a request finished or start repeated parallel downloads.

A subsequent browser search successfully displayed an individual candidate association; see [rs8067974 lookup](stage6_tpmi_rs8067974_lookup.md). It is a single rounded, exploratory estimate in a broader sleep-disorders phenotype, not score validation.

A manual browser lookup for a second candidate, rs12536335, found a rounded sleep-disorders estimate close to zero with substantial uncertainty. See [FOXP2-region lookup](stage6_tpmi_foxp2_lookup.md); its effect allele remains unconfirmed and the outcome differs from insomnia.

A third manual lookup, [rs10865954](stage6_tpmi_rbm6_lookup.md), returned a rounded sleep-disorder estimate with p = 0.084. Its project-coordinate lift matched, but effect-allele semantics and multi-allelic representation remain unresolved. These three postselected single-marker checks are not a score-level test.

Two further rounded candidate estimates, rs7734045 and rs1032524, appear in the [two-marker lookup](stage6_tpmi_two_candidate_lookup.md). These remain exploratory single-marker findings with unconfirmed effect-allele orientation.

A further postselected lookup, [rs1909650](stage6_tpmi_olfm4_lookup.md), returned nominal p = 0.026 for broad sleep disorders; effect-allele direction remains unconfirmed and no multiple-testing-adjusted or score-validation claim is made.

A further manual screen for [rs1942262](stage6_tpmi_tcf4_lookup.md) returned a rounded p = 0.14; its allele direction remains unverified. It is not a whole-score validation.

Two further chromosome 2 variants were matched to their project coordinates and alleles in the [TPMI lookup](stage6_tpmi_chr2_lookup.md). Both have statistically unclear rounded sleep-disorder estimates; source effect alleles and full-score validation remain unresolved.

A further browser result for [rs11790060](stage6_tpmi_chr9_lookup.md) was nominally associated with broad sleep disorders (p = 0.011), but allele direction and the multi-allelic contrast remain unverified and this postselected lookup does not validate a score. Two other rsID searches showed no matches; that is not evidence of variant absence.
