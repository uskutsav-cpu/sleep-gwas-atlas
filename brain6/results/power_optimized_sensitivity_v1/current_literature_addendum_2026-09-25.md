# Power-screen literature addendum — 2026-09-25

This append-only addendum records one 2026 publication surfaced in a fresh
source check. It does not alter the hash-bound outcome-blinded candidate screen
or its promotion rule.

## Portas et al. 2026 device-measured sleep GWAS

Portas et al. report four accelerometer-derived sleep traits in 80,013
European-ancestry UK Biobank participants: night-time sleep duration, sleep
efficiency, and inferred REM/NREM duration. The study considered over 9.8
million common variants. Its GWAS Catalog accessions are GCST90824104 through
GCST90824107 (group GCP001581). The paper reports GRCh37-to-GRCh38 conversion
for a GTEx colocalization analysis; that does not establish the build of each
released GWAS file. Catalog metadata currently labels GCST90824104 as not yet
associated with a peer-reviewed publication and “Full Summary Statistics: Not
available,” while showing an FTP link. No file was downloaded or admitted, and
its trait/file schema, effect fields, checksum, and exact sample overlap were
not verified.

**Decision: exclude as a long-sleep replacement or power-optimized sensitivity.**
The objectively measured nighttime-duration phenotype is not the self-reported
≥9-hour binary tail. Its N=80,013 is also materially below the current long-sleep
GWAS N=339,926 and the related continuous-duration GWAS N=446,118. It shares the
UK Biobank source with the existing analyses. Summary-statistic availability
does not overcome the phenotype and power criteria. This dataset may be
reconsidered only for a separately defined device-sleep question; it is not
eligible for the present replacement screen.

Sources: [Portas et al., Nature Communications (2026)](https://www.nature.com/articles/s41467-026-71252-y), [GWAS Catalog GCST90824104](https://www.ebi.ac.uk/gwas/studies/GCST90824104), [GWAS Catalog summary-statistics access](https://www.ebi.ac.uk/gwas/downloads/summary-statistics).
