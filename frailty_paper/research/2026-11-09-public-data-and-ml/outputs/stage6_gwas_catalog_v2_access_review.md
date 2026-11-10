# Updated public-source audit using GWAS Catalog API v2

The previous availability checks encountered HTTP 410 from the Catalog's retired API. The Catalog now directs users to [REST API v2](https://www.ebi.ac.uk/gwas/docs/news/rest-api-v2-migration-guide/). On 28 October 2026, publication-filtered v2 study queries gave the following results:

| Publication | PMID | Result from this query |
|---|---:|---|
| Edelson et al., periodic limb movements | 35670608 | Two curated entries, GCST90270120 and GCST90270121; both say full summary statistics unavailable in the Catalog. |
| Tseng et al., UKB electronic FI | 37523034 | Zero study records. |
| Kwak et al., Korean frailty | 41733842 | Zero study records. |
| Borhan et al., Canadian frailty phenotype | 42020437 | Zero study records. |
| Mak et al., metabolomics and modified FI | 37184129 | Zero study records. |

These are query results, not proof that the studies or complete statistics do not exist in another repository, under another publication record or later release. Zero records are distinct from an explicitly marked `false`, and neither is a negative scientific result. The previously identified Mak 41-item FI excludes eight mental-wellbeing items, including insomnia and fatigue, but this query did not find an available file; see `mak2023_sleep_excluded_fi_lead.md`. A public full signed GWAS and item coding remain necessary before the planned matched LDSC analysis can run. An exact query for the literal disease trait `Frailty index` returned two unrelated existing Catalog entries, illustrating that the API is working; it was not an exhaustive search over all frailty terminology.

`stage6_gwas_catalog_v2_public_access_audit.py` makes the five publication queries and saves the raw JSON, page counts, URLs and SHA-256 receipts in `stage6_gwas_catalog_v2_public_access_audit.json` and per-PMID JSON files. All five fit within the requested page size. The API queries the Catalog's curated data and is not a comprehensive catalogue of external repositories; no participant data were retrieved. Relevant source statements and independent cohort constraints remain in `independent_frailty_gwas_feasibility.md`, `korean_frailty_gwas_validation_feasibility.md`, and `independent_outcome_validation_cohort_feasibility.md`.
