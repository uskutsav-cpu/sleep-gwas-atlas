# MVP insomnia lookup of two earlier FI variants

The two variants highlighted by [Mekli et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6258126/) were each present once by rsID in the previously verified MVP clinical-record insomnia file. Neither showed nominal insomnia association in this lookup:

| Prior FI variant | Earlier FI result | MVP insomnia p | MVP source alleles and effect |
|---|---|---:|---|
| rs6765037 | HRS discovery p = 2.14 × 10⁻⁸; ELSA p = 0.204 | 0.871 | C versus T, OR 1.002 |
| rs7134291 | HRS p = 1.85 × 10⁻⁶; ELSA p = 0.0343 | 0.465 | G versus A, OR 0.9948 |

This is a post-selected historical-locus lookup, not a prespecified two-trait replication. Both variants were chosen from an earlier FI paper, and its second cohort's rs7134291 result was nominal; neither result establishes shared genetic causation. The reported p-values do not rule out an insomnia association below the available power. Sleep ascertainment and sample composition differ, and participant or family overlap has not been independently audited.

The original FI analysis reported the A allele for both variants but the second allele and strand conventions were not independently established here. The MVP file reports C/T at rs6765037 and G/A at rs7134291, so no cross-study effect-direction or exact molecular-identity claim is made. The source's build-specific rsID and position are recorded in the [rows](mvp_source_rows.tsv); an exact independent reference/build check and full allele reconciliation would be needed before directional interpretation. The MVP odds ratios are rounded, and its supplied SE fields for these records are unavailable. The FI in both earlier cohorts includes restless sleep, so this is not an item-excluded frailty comparison or a mechanism for GRIN2B.

[Protocol](PROTOCOL.md), [scan code](scan.py) and [receipt](receipt.json) record the lookup. The full 19,703,815-row MVP file passed gzip reading and its compressed SHA-256 matched the previously saved receipt. A matching file hash does not establish participant independence, source-model validity, causal sharing or a passed local QC gate.
