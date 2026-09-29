# HLMA donor and library supplement audit

**Audit date:** 2026-09-24 UTC
**Study:** Multimodal cell atlas of the ageing human skeletal muscle (PMID 38649488; DOI 10.1038/s41586-024-07348-6)
**Scope:** integrity and descriptive structure of PMC Supplementary Tables 1–2 only. No HLMA processed matrices or other omics objects were acquired or analyzed.

## Acquisition and identity checks

The files were retrieved from the documented public PMC Article Datasets HTTPS distribution, using the PMC article-version metadata JSON for `PMC11062927.1`. The saved metadata identifies the version as open access under CC BY and lists both media objects with MD5 checksums. The JSON is recorded at `frailty_paper/manifests/HLMA_PMC_Article_Datasets_PMC11062927.1.json` (SHA-256 `9f02f18b33152ced30f5b109f15b136434e2e07eb6af9067fcfc0e65756f35d8`). The article-level license is not assumed to apply to distinct objects on the HLMA processed-data portal.

| File | Bytes | PMC MD5 | SHA-256 | XLSX package | Worksheet |
|---|---:|---|---|---|---|
| `HLMA_Supplementary_Table_1_Donor_Information.xlsx` | 12,117 | `0d2b1e477c72a8afa422cf1500a41757` | `55a69c768a6862ea5e1317c81a96f75dee7733b4ae406df8b22a67d630fe0773` | Valid | `Supplementary Table 1` |
| `HLMA_Supplementary_Table_2_Library_Information.xlsx` | 22,313 | `58fd305b32267f318abe37ee0987ae26` | `f5c4f3f0719802b68d9c187176003d6edea9490528267c477eb11d09aa1da33a` | Valid | `Supplementary Table 2` |

Both files matched the MD5 published in the PMC metadata and passed ZIP/XLSX package integrity checks. Their external-storage paths and SHA-256 values are registered in `frailty_paper/manifests/all_acquired_resources.tsv`.

## Descriptive findings

Supplementary Table 1 contains 31 donor records and 12 columns. Aggregate counts are 17 male and 14 female records; 18 European-cohort and 13 Asian Chinese-cohort records; the listed ages span 15–99 years. Individual-level dates and measurements are not reproduced here.

Supplementary Table 2 contains three assay sections and 190 library rows:

| Assay | Library rows | Distinct donor codes | Distinct samples | Cell/nucleus count sum | Missing count values |
|---|---:|---:|---:|---:|---:|
| scRNA-seq | 38 | 7 | 11 | 79,649 | 0 |
| snRNA-seq | 106 | 22 | 27 | 212,774 | 0 |
| snATAC-seq | 46 | 16 | 18 | 95,021 | 0 |

All three aggregate cell/nucleus totals agree with the primary article. Donor codes in the snRNA-seq and snATAC-seq sections match Table 1. The scRNA-seq section includes code `YM5`, which has no matching donor row in Table 1. This discrepancy is unresolved; no alias or identity is inferred.

## Interpretation and remaining limits

The supplements resolve donor descriptions and modality/library/sample metadata. They do **not** establish per-file CNP accession mapping for the HLMA processed-data portal, the coordinate build of a selected snATAC peak object, or the portal-specific reuse terms. The PMC article-version CC BY declaration is not extended to separately hosted portal objects. No locus overlap, expression analysis, cell-type enrichment, or mechanistic claim follows from this metadata audit. Any future use of portal objects remains gated on eligible project loci and verification of the exact object, build and reuse terms.
