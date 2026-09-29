# HLMA portal-to-NGDC OMIX crosswalk recheck

**Checked:** 2026-09-24 UTC
**Scope:** public listing metadata and downloadable archive IDs only. No archive or processed matrix was downloaded.

The official live [HLMA download page](https://db.genomics.cn/cell/hlma/download) now exposes CNGBdb/NGDC OMIX archive URLs in the link targets for its snRNA-seq count-matrix rows. The direct page links identify nine European-labeled P-sample archives under `OMIX006022` and a `snRNA-China.tar.gz` link under `OMIX004308`. Official release records are [OMIX006022](https://ngdc.cncb.ac.cn/omix/release/OMIX006022) and [OMIX004308](https://ngdc.cncb.ac.cn/omix/release/OMIX004308). Both report open-access data availability under BioProject `PRJCA017356`.

| HLMA portal filename | Portal link file ID | Official OMIX file title and displayed release size | Finding |
|---|---|---|---|
| `P29_snRNA_count_matrix.tar.gz` | `OMIX006022-09` | P29; 43.23 MB | Sample label agrees; listed sizes approximately agree. |
| `P27_snRNA_count_matrix.tar.gz` (25.3M row) | `OMIX006022-08` | P27; 25.33 MB | Sample label agrees; listed sizes approximately agree. |
| `P5_snRNA_count_matrix.tar.gz` | `OMIX006022-07` | P5; 277.81 MB | Sample label agrees; listed sizes approximately agree. |
| `P3_snRNA_count_matrix.tar.gz` | `OMIX006022-06` | P3; 236.04 MB | Sample label agrees; listed sizes approximately agree. |
| `P23_snRNA_count_matrix.tar.gz` | `OMIX006022-05` | P23; 86.84 MB | Sample label agrees; listed sizes approximately agree. |
| `P27_snRNA_count_matrix.tar.gz` (45.1M row) | `OMIX006022-04` | P17; 45.14 MB | **Unresolved label conflict:** the portal filename says P27, whereas the official OMIX file record says P17. The size agrees, but it does not resolve which label is correct. |
| `P13_snRNA_count_matrix.tar.gz` | `OMIX006022-03` | P13; 25.71 MB | Sample label agrees; listed sizes approximately agree. |
| `P26_snRNA_count_matrix.tar.gz` | `OMIX006022-02` | P26; 17.13 MB | Sample label agrees; listed sizes approximately agree. |
| `P21_snRNA_count_matrix.tar.gz` | `OMIX006022-01` | P21; 144.81 MB | Sample label agrees; listed sizes approximately agree. |
| `snRNA-China.tar.gz` | `OMIX004308-06` | Generic “snRNA-seq count matrix”; 4.76 GB | The portal link resolves to this OMIX file ID, but NGDC does not give the China filename in its record and the displayed sizes differ (5.7G vs 4.76 GB). Do not assume archive equivalence without a byte-level check. |

This resolves the **linked OMIX file IDs** for these ten snRNA portal rows; it does not establish file contents because no archives were downloaded. One P27/P17 sample-label conflict and one China archive size/title mismatch remain. The six peak matrices, the scRNA count archive, and six H5AD products (13 further rows in the 23-row targeted inventory) still have no verified OMIX file crosswalk in this audit. The 46 fragment listings were not expanded into an individual accession inventory. The global CNP accession list remains a study-level list, not a per-file CNP mapping. No coordinate build or reuse license is resolved by the OMIX links or their open-access availability field.

## Second pass: public API listing and release-record comparison (2026-09-24)

The portal's public client bundle calls `https://db.genomics.cn/cdcp_hlma/api/download/get_download` with no query parameters to populate its download table. A normal HTTP GET to that same endpoint returned HTTP 200 (`application/json`, 13,050 bytes). The exact response is preserved as `../manifests/hlma_portal_download_metadata_2026-09-24.json` (SHA-256 `bcacfd5f78af7745470ef49d0259ef5db1245fa25ec616e34f026ddf39c1c94d`). It contains five sections and 69 listing rows: 10 snRNA count, 46 snATAC-fragment, 6 peak-matrix, 1 scRNA-count, and 6 H5AD rows.

The page API supplies `file_path` for all 69 rows. Exact filename-and-size matching links all 23 targeted processed-object rows to a portal-reported archive path and maps the 46 fragment listings to ten portal-reported OMIX paths. The reproducible outputs are `../manifests/hlma_processed_file_inventory.tsv`, `../manifests/hlma_fragment_archive_crosswalk.tsv`, and builder `../scripts/46_build_hlma_portal_archive_crosswalk.py`. These are URL/identifier crosswalks only: no archive was fetched, and the link does not establish that the named portal file is present as an archive member.

Cross-checking the linked release pages exposes further inconsistencies that remain open:

| Portal listing row | Portal `file_path` ID | Official release-page description | Finding |
|---|---|---|---|
| `YM2_ST_snATAC_seq_1.fragments.tsv.gz` | `OMIX006021-05` | ATAC-seq fragments for P5 | Sample-code conflict; do not infer that YM2 aliases P5. |
| `P5_TA_snATAC_seq_1.fragments.tsv.gz` | `OMIX006021-03` | ATAC-seq fragments for P29 | Sample-label conflict. |
| `P29_GM_snATAC_seq.fragments.tsv.gz` | `OMIX006021-01` | ATAC-seq fragments for P26 | Sample-label conflict. |
| `P26_ST_snATAC_seq_1.fragments.tsv.gz` | `OMIX006021-06` | ATAC-seq fragments for P23 | Sample-label conflict. |
| `P27_snRNA_count_matrix.tar.gz` (45.1M row) | `OMIX006022-04` | snRNA count matrix for P17 | Existing P27/P17 conflict remains. |
| `Myofiber_sn_RNA.h5ad` (560M) | `OMIX004308-05` | snRNA H5AD for Myonucleus; 6.01 GB | Portal filename, cell label, and size do not agree with the release description. |
| `snRNA-China.tar.gz` (5.7G) | `OMIX004308-06` | Generic snRNA count matrix; 4.76 GB | Existing title/size discrepancy remains; archive bytes are needed to assess identity. |

The four P-sample fragment-ID disagreements, along with portal-row counts that differ from the corresponding release sample counts for `OMIX004305-02`, `-03`, `-04`, and `OMIX006021-06`, are recorded row-by-row in `../manifests/hlma_fragment_archive_crosswalk.tsv`. The count difference is metadata only and does not establish a missing or duplicated fragment. Official OMIX release descriptions also characterize `OMIX004308-01` as one combined Vascular/Immune/MuSC H5AD submission, while the portal links three object rows to that archive path; archive-member identities remain unverified. The FAPs/Stromal and All-CellType/Skeletal-muscle labels/sizes are broadly compatible but do not prove file identity.

The `per_file_accession_mapping` inventory field remains `UNKNOWN` for all 23 objects: OMIX archive identifiers are not CNP BioProject accessions. Per-file CNP mapping, object-level coordinate build, archive byte identity, donor identity for conflicting rows, and object-specific reuse terms remain unresolved. The API's `public_type=public` field and the release pages' `Open-access` label are access metadata, not reuse licenses.

The snapshot and builder can be regenerated/checked with `make -C frailty_paper audit-hlma-portal-crosswalk`. Official sources: [HLMA download listing](https://db.genomics.cn/cell/hlma/download), [HLMA portal API](https://db.genomics.cn/cdcp_hlma/api/download/get_download), [OMIX004305](https://ngdc.cncb.ac.cn/omix/release/OMIX004305), [OMIX006021](https://ngdc.cncb.ac.cn/omix/release/OMIX006021), [OMIX004308](https://ngdc.cncb.ac.cn/omix/release/OMIX004308), [OMIX006022](https://ngdc.cncb.ac.cn/omix/release/OMIX006022).
