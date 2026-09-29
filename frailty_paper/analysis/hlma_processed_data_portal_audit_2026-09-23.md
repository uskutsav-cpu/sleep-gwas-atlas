# HLMA processed-data portal audit

Accessed: 2026-09-23 13:31 UTC. Source was the live official [HLMA download page](https://db.genomics.cn/cell/hlma/download), read in the Codex in-app browser. The visible page lists five processed-data categories: snRNA-seq count matrices (10 files), snATAC-seq fragments (46 files), snATAC-seq peak matrices (6 files), an scRNA-seq count-matrix archive (1 file), and processed sc/snRNA H5AD objects (6 files). The targetable count/peak/H5AD objects (23 entries, including displayed sizes) are transcribed to frailty_paper/manifests/hlma_processed_file_inventory.tsv. The 46 fragment files were not inventoried individually because the matrix/peak/H5AD products provide a more targeted processed-data starting point and no eligible project locus currently justifies bulk acquisition.

The official portal currently exposes category-specific file names and display sizes but the accessible listing does not expose a per-file CNP accession crosswalk or checksum. Four project accessions are registered at the study level: CNP0004394, CNP0004395, CNP0004494 and CNP0004495. The Nature article reports 31 biopsied participants and a QC multimodal dataset of 22 people, with 212,774 snRNA-seq nuclei, 79,649 scRNA-seq cells and 95,021 snATAC-seq nuclei; it points to the HLMA portal for processed data. The article’s CC BY 4.0 status does not by itself verify reuse terms for repository-hosted objects. Portal data-sharing/terms links were visible, but their page text was not retrievable in the available web reader, so dataset-specific redistribution/reuse terms remain unverified.

No file was downloaded. No checksum, per-object genome build, or per-file project accession is claimed. In particular, the snATAC peak-coordinate build must be read from the selected object metadata before any variant overlap. Because no eligible frailty locus is frozen, object-level acquisition and cell-context analyses remain deferred.

## Evidence sources

- HLMA processed-data listing: https://db.genomics.cn/cell/hlma/download
- HLMA sample portal: https://db.genomics.cn/cell/hlma/sample
- Primary article: https://www.nature.com/articles/s41586-024-07348-6

## Publisher supplement check — 2026-09-23

The primary paper’s PMC full text exposes Supplementary Table 1 (donor information; 11.8 KB XLSX) and Supplementary Table 2 (library information; 21.8 KB XLSX). These tables are relevant to donor-level phenotype and library-to-modality/sample mapping. Ordinary link opening through PMC returned a “Checking your browser before accessing pmc.ncbi.nlm.nih.gov” reCAPTCHA interstitial for both files; no spreadsheet bytes were obtained or inspected. Nature’s direct supplement route previously redirected the web reader to its identity-protection page. Therefore donor/library fields, sample-to-library links, and any accession crosswalk remain unverified. This is a retrieval challenge, not evidence that the supplementary data are restricted. Do not bypass it; retry only through an ordinary permitted publisher route.

Source: [PMC full text of the primary Nature paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11062927/), Supplementary Information lines listing table names, descriptions, and sizes.


## PMC Article Datasets supplement acquisition addendum — 2026-09-24

The two cited workbooks were subsequently acquired through the documented public PMC Article Datasets AWS route, using the official article-version metadata JSON for `PMC11062927.1` (saved at `frailty_paper/manifests/HLMA_PMC_Article_Datasets_PMC11062927.1.json`, SHA-256 `9f02f18b33152ced30f5b109f15b136434e2e07eb6af9067fcfc0e65756f35d8`). Its metadata identifies the article as open access under `CC BY` and lists both workbook media entries with MD5 values. The article-version metadata is not a license determination for the separate HLMA processed-object portal.

Both XLSX files were retrieved over public HTTPS from `pmc-oa-opendata.s3.amazonaws.com`, passed ZIP package integrity checks, and matched PMC’s published MD5: Table 1, 12,117 bytes, MD5 `0d2b1e477c72a8afa422cf1500a41757`, SHA-256 `55a69c768a6862ea5e1317c81a96f75dee7733b4ae406df8b22a67d630fe0773`; Table 2, 22,313 bytes, MD5 `58fd305b32267f318abe37ee0987ae26`, SHA-256 `f5c4f3f0719802b68d9c187176003d6edea9490528267c477eb11d09aa1da33a`. They are stored under the configured external frailty-data root at `frailty_paper/data/single_cell/hlma_supplementary/`.

Table 1 has one worksheet, 31 donor records and 12 columns; aggregate counts are 17 male, 14 female, 18 European-cohort and 13 Asian Chinese-cohort records, with listed ages 15–99. Table 2 has one worksheet and three modality sections: scRNA-seq, 38 library rows / 7 donor codes / 11 samples / 79,649 total cells; snRNA-seq, 106 rows / 22 donors / 27 samples / 212,774 nuclei; snATAC-seq, 46 rows / 16 donors / 18 samples / 95,021 nuclei. All library count cells were populated. The cell/nucleus totals agree with the primary article. Donor code `YM5` appears among scRNA-seq library records but has no matching row in Table 1; the donor crosswalk is unresolved and no alias is inferred.

These supplements resolve descriptive donor and library metadata only. They do not resolve per-file CNP accession mapping, processed-object coordinate build, or reuse terms for data hosted by the separate HLMA portal; no processed matrices or omics analyses were undertaken.


## NGDC OMIX archive-link crosswalk — 2026-09-24

A live read of the official HLMA download page exposed OMIX archive IDs on the download links for all ten listed snRNA count archives. Nine Europe-labeled P-sample links map to OMIX006022 file IDs; the `snRNA-China.tar.gz` link maps to OMIX004308-06. Official release records identify their file rows and mark data accessibility “Open-access”: [OMIX006022](https://ngdc.cncb.ac.cn/omix/release/OMIX006022) and [OMIX004308](https://ngdc.cncb.ac.cn/omix/release/OMIX004308).

This is a link-level crosswalk, not byte-level file identity. The portal names one archive `P27_snRNA_count_matrix.tar.gz` but links it to OMIX006022-04, whose NGDC record names P17. The `snRNA-China` row links to OMIX004308-06, whose NGDC record is generically titled “snRNA-seq count matrix” and displays a different size (portal 5.7G; release 4.76 GB). Both discrepancies remain unresolved. The other 13 objects in the 23-row targeted inventory do not have a verified OMIX crosswalk in this audit. Global CNP accessions remain study-level only; no per-file CNP mapping, coordinate build or reuse license follows from the OMIX links. No archive was downloaded. See `frailty_paper/analysis/hlma_portal_omix_crosswalk_2026-09-24.md` and the updated `frailty_paper/manifests/hlma_processed_file_inventory.tsv`.
