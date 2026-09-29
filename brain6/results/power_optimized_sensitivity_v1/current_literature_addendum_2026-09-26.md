# Power-screen literature addendum — 2026-09-26

This dated source recheck follows the 2026-09-25 literature note and does not modify the frozen outcome-blinded candidate screen or its promotion rule.

## Portas et al. 2026 device-measured sleep GWAS

The peer-reviewed *Nature Communications* paper reports four accelerometer-derived traits (night-time sleep duration, sleep efficiency, and inferred REM/NREM duration) in 80,013 European-ancestry UK Biobank participants, with over 9.8 million common variants. Its GWAS Catalog accessions are GCST90824104–GCST90824107 (GCP001581). The article’s GRCh37-to-GRCh38 statement refers to a GTEx colocalization conversion, so it does not establish the build of an association file.

The paper’s data statement says its full summary statistics are deposited in the GWAS Catalog. On 2026-09-26, Catalog study pages for GCST90824104, GCST90824105, and GCST90824107 still display “Full Summary Statistics: Not available” and incorrectly retain the “not yet associated with a peer-reviewed publication” label despite the article being published on 2026-04-01. The GCST90824106 page was not retrievable in this check. No full GWAS file was downloaded, so variant schema, effect fields, exact file-level terms, build, MD5/SHA256, and reference overlap remain unverified. Catalog documentation says full files and their metadata are hosted on FTP when available; publication deposition statements alone are not a verified download.

**Decision: exclude from the long-sleep power screen.** Night-time accelerometer duration is not the self-reported >=9-hour binary tail, and N=80,013 is below both canonical long sleep (N=339,926) and the harmonized continuous-duration candidate (N=446,118). No material power gain is expected, independent of downstream associations. Do not acquire or harmonize this as a long-sleep replacement. Reconsider only for a separately specified device-sleep question or once a usable full-statistics file can be directly verified.

Sources: [Portas et al., Nature Communications (2026)](https://www.nature.com/articles/s41467-026-71252-y), [Catalog GCST90824104](https://www.ebi.ac.uk/gwas/studies/GCST90824104), [Catalog GCST90824105](https://www.ebi.ac.uk/gwas/studies/GCST90824105), [Catalog GCST90824107](https://www.ebi.ac.uk/gwas/studies/GCST90824107), and [GWAS Catalog summary-statistics documentation](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/).
