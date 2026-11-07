# Unsent source-metadata request: physical frailty component GWAS

Proposed recipient: Jiajin Chen, `jiajinchen@xmu.edu.cn`, as listed in the [published article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12972468/). Draft only; no message has been sent. Sender name and any affiliation require the user's choice.

Subject: Technical and provenance details for Zenodo record 14011550

Dear Dr. Chen,

I am auditing the five physical-frailty component summary-statistic files deposited at Zenodo record 14011550 for a reproducible genetic analysis. The public record and `0_Readme.txt` identify UK Biobank and Regenie v3.5, but do not specify several details needed before allele harmonization. Could you please clarify, for each component file:

If a public methods supplement, repository README or analysis code already documents these details, a link would be very helpful; the analysis is prioritizing openly accessible data without account registration. Full-file checks have confirmed the depositor's checksums and GRCh37 agreement for matched reference rows, but these do not establish the study design.

- the genome assembly and chromosome/position conventions;
- ancestry definition, sample exclusions, and variant-level QC/imputation filters;
- the source cohort(s) and exact generating study or analysis;
- phenotype coding, case/control definitions, and whether the effect is logistic-regression log-odds, linear-model beta, or another scale;
- whether `BETA` and `A1FREQ` refer to `ALLELE1` in these processed files, and how multi-allelic records and repeated rsIDs were represented;
- covariates and Regenie command/model settings;
- per-variant sample-size semantics and any missingness; and
- participant overlap among component files and with other published UK Biobank GWAS.

The README cites the heart-failure study, which reports 49,530 UK Biobank participants and refers to an earlier physical-frailty GWAS, while the five deposited files report total sample sizes around 397,000–407,000. Please confirm whether the cited paper generated these summary statistics or cited them as an input, and identify the earlier source if applicable.

Thank you,

[Name]
