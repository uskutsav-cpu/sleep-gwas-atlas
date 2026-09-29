# HFRS 109-code crosswalk to FinnGen DF13 — 2026-09-27

## Question

Can the official public DF13 endpoint catalog identify standard FinnGen endpoints containing the ICD-10 categories used by Mak et al.'s custom HFRS score? This is a source-discovery and definition crosswalk only; it does not qualify any standard endpoint as the paper's HFRS GWAS.

## Sources and integrity

- Mak et al., Nature Aging 2025 article and its official supplementary workbook: <https://www.nature.com/articles/s43587-025-00925-y>
- Direct supplementary workbook: <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43587-025-00925-y/MediaObjects/43587_2025_925_MOESM1_ESM.xlsx>
  - bytes: 1,059,282
  - SHA-256: `123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625`
  - Sheet `ST18` supplies 109 ICD-10 code categories, weights and descriptions.
- FinnGen's official clinical-endpoint list: <https://www.finngen.fi/en/researchers/clinical-endpoints>
- Linked public DF13 endpoint workbook: <https://www.finngen.fi/sites/default/files/inline-files/FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx>
  - bytes: 864,406
  - SHA-256: `0004624f727b96b4a56c6fac75ff182eb8bc3daf78cdf690bdad276b31fca6c7`
  - The worksheet's stored dimension is `A1:A1`; the scan therefore reads sheet XML rows directly and does not trust the dimension.

## Method and result

The reproducible scanner `../scripts/crosswalk_hfrs_codes_to_finngen_catalog.py` extracts the 109 code/weight/description rows from `ST18` and searches for exact code tokens in FinnGen DF13's `OUTPAT_ICD`, `HD_ICD_10_ATC`, `HD_ICD_10`, and `COD_ICD_10` fields. It does not expand ICD hierarchies, interpret bracket expressions/ranges, or infer endpoint equivalence.

Of 109 HFRS code rows, 99 had at least one literal exact-token occurrence in a public endpoint definition, and 10 did not. The latter codes are `I67`, `W18`, `M81`, `E16`, `N28`, `G31`, `H91`, `E05`, `M80`, and `M48`. A zero literal match does not establish that a code or descendant is absent; it may be encoded using an unexpanded pattern or hierarchy. The row-level crosswalk is `hfrs_df13_code_crosswalk_2026-09-27.tsv` (109 rows; SHA-256 `b6eae124e0f73aad65a0f5a0a1584822f9a15106852c3bac8904ae318b507e1e`).

## Interpretation and eligibility

The crosswalk confirms that many component ICD-10 categories appear in standard FinnGen endpoint definitions. Those endpoints can group codes differently, use different inclusion/exclusion logic, represent different case/control sets, and lack the paper's person-level sum of 109 weights across the specified age and follow-up window. Component endpoint summary statistics cannot be summed to reconstruct a GWAS of a weighted continuous score. Thus this check does not establish that the custom FinnGen HFRS or HFRS-without-dementia GWAS is present in DF13, nor make a component endpoint a valid substitute. No proxy was selected and no HFRS analysis was run.

The official supplementary workbook includes thresholded GWAS results and the score definition, but not the full FinnGen HFRS summary statistics. In Table 18, the FinnGen and UKB `N_individual` columns are 519,200 and 407,643, respectively; the article reports GWAS sample sizes of 500,737 and 407,463. These counts differ by 18,463 and 180. The reason is not established here, so Table 18 denominators are not treated as GWAS Ns. The exact custom full-statistics object, endpoint/source identifier, and checksum remain unresolved. Bulk FinnGen access still follows the official registration-and-email route documented in `hfrs_official_route_recheck_2026-09-27.md`; no form was submitted.
