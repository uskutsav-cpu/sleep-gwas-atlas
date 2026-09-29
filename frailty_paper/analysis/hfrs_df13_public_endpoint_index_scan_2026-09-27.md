# FinnGen DF13 public endpoint-index scan — 2026-09-27

## Scope and source

To extend the HFRS access-route review, the official FinnGen clinical-endpoint page was checked for its linked DF13 endpoint workbook. The downloaded workbook was saved temporarily outside the repository as `/private/tmp/frailty_hfrs_df13_endpoints.xlsx`.

- Official endpoint page: <https://www.finngen.fi/en/researchers/clinical-endpoints>
- Workbook URL: <https://www.finngen.fi/sites/default/files/inline-files/FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx>
- Download size: 864,406 bytes
- SHA-256: `0004624f727b96b4a56c6fac75ff182eb8bc3daf78cdf690bdad276b31fca6c7`
- Workbook sheet: `Sheet 1`

## Scan and result

The workbook's stored worksheet dimension incorrectly reports `A1:A1`. To avoid trusting that dimension, the full worksheet XML was iterated with an explicit upper bound of 100,000 rows and 20 columns. This reached 5,043 non-empty rows. A case-insensitive scan of every populated cell for `frailty`, `HFRS`, and `hospital frailty` found zero matches.

## Interpretation

This is evidence only about literal labels in this particular public endpoint-index workbook. It does not establish that the custom HFRS phenotype is absent from FinnGen, from an opaque-named object, or from a separately delivered custom analysis. Nor does it make a standard endpoint a valid substitute for the paper's weighted 109-code HFRS. The exact custom HFRS and HFRS-without-dementia full-statistics objects, endpoint codes, and checksums remain unidentified. No form was submitted and no protected or bulk files were accessed.

## Reproducibility note

The downloaded workbook is temporary and is not part of the analysis inputs. The source URL, byte count, checksum, scan terms, explicit scan bounds, row count, result, and interpretation limits are recorded here. The preceding route assessment remains in `hfrs_official_route_recheck_2026-09-27.md`.

The standard-library scanner is `frailty_paper/scripts/scan_finngen_endpoint_index_xlsx.py`; rerunning it against the downloaded workbook reports 5,043 worksheet rows, 5,043 non-empty rows, and zero matching rows.
