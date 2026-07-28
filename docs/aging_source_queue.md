# Aging-source verification queue

This is a discovery log for selected aging traits. A row here is **not** a
source-registry entry and does not make a trait downloadable or `CURATED`.
It records only sources whose primary release page and metadata were checked
on 28 July 2026; the registry is updated only after a local download passes
its own integrity and header checks.

| Trait | Official release | Verified facts | Still required before curation |
| --- | --- | --- | --- |
| `longevity` | [GCST008598](https://www.ebi.ac.uk/gwas/studies/GCST008598), [`Results_90th_percentile.txt.gz`](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008598/Results_90th_percentile.txt.gz) | Downloaded 28 July 2026: gzip integrity passed; 141,843,112 bytes; SHA-256 `99353dfd78f9cb86b35f6cea4fb86ccebed781064d03a8495f24ee56da7c6ac7`. README specifies GRCh37, effect/non-effect alleles, log-OR beta, and per-SNP effective N. | **Blocked:** the actual file’s SNP IDs are only `chromosome:position`; obtain a documented GRCh37 rsID mapping before standardization/LDSC. |
| `parental_lifespan` | [GCST009890](https://www.ebi.ac.uk/gwas/studies/GCST009890), [`lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz`](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009890/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz) | README specifies GRCh37, effect/reference alleles, effect frequency, log-hazard-protection beta, SE, P, and INFO. The catalog's MD5 is `9209bf29b5cf281bba35934dcf708538`. | Download, check the supplied MD5 plus local SHA-256/size, validate header/sentinel, and decide the defensible per-variant N rule. |
| `healthspan` | [GCST007406](https://www.ebi.ac.uk/gwas/studies/GCST007406), [`GCST007406_buildGRCh37.tsv`](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007406/GCST007406_buildGRCh37.tsv) | Catalog metadata names a European UK Biobank subset of 300,447, GRCh37/hg19 coordinates, and beta/SE/P/effect-allele fields. Its supplied MD5 is `08f9d5861b99d0b0c00d87822eaf995c`. | Download, check supplied MD5 plus local SHA-256/size, validate header/sentinel, and gzip-standardize the plain TSV reproducibly. |

The releases above are compatible in principle with the atlas's EUR/GRCh37
requirements, but each retains its source-specific phenotype definition and
sample-size semantics. They must not be treated as interchangeable ageing
measures or pushed through the generic ZIP materializer.
