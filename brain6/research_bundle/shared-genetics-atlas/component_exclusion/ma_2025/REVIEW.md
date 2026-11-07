# Ma et al. frailty supplementary-data audit

Checked 6 November 2026. The [article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11979021/) analyzes the full five-component Fried phenotype, including exhaustion. Its data statement identifies restricted UK Biobank data. The supplementary workbook is publicly downloadable without an account from the [publisher](https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs11357-024-01395-7/MediaObjects/11357_2024_1395_MOESM2_ESM.xlsx), despite a browser challenge at the separate PMC mirror.

The complete downloaded workbook is retained as `source_supplement.xlsx`; it was read without modification. Two separate publisher retrievals yielded 3,668,249 bytes with identical SHA-256 `40222eac737c6806e1978f7e8835dce36be44120916ecd3070a9385186a75f3c`. The [receipt](receipt.json) records all 26 sheet names, worksheet row counts, content classifications, exhaustion-label locations and the source URL.

OOXML row and cell inspection found:

| Sheets | Contents | Role in this project |
|---|---|---|
| 1–2 | Frailty criteria and descriptive counts | Confirms exhaustion remains part of the phenotype; its label is at Table 1 cell A5 and Table 2 cell A12. |
| 3–4 | Exposure correlations and coding | Defines epidemiological exposures, rather than signed SNP effects. |
| 5–22 | Exposure-wide associations and subgroups | These are environmental/exposure analyses, not genome-wide per-variant association statistics. Stratification by polygenic score does not turn them into SNP-level data. |
| 23–26 | Attributable-fraction summaries | Population-level exposure summaries; cannot support LDSC, locus colocalization or modified-score reconstruction. |

No inspected sheet contains a full per-variant GWAS table. There are no cells consisting of an rsID and no external OOXML relationships. Absence of rsIDs alone would not establish absence of GWAS data; the content and headers across the worksheets provide the substantive distinction. No external table-download relationship was found in this workbook. This audit does not establish that no other public release exists.

**Decision:** Do not add this workbook to genetic-correlation estimation or ML training as a physical-frailty GWAS. It does not resolve the matched exhaustion-deleted phenotype requirement. Use the existing [explicit physical-function index](../physical_factor_workaround/INDEX_RESULTS.md) as the distinct-construct sensitivity, while retaining the missing-data limitation. The [linear-reconstruction audit](../identifiability/LINEAR_ROUTE.md) specifies what compatible summary information would make an additive reconstruction possible.
