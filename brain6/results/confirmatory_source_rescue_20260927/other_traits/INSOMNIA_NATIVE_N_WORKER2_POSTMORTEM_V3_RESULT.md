# Worker-2 postmortem v3 result: structural empty-input locus

**Diagnostic-only result.** The separately frozen and committed worker-2 replay (`c5fb4d7c`) reproduced the v2 failure at exactly one locus: **950**. Its failure text is `Less than 3 SNPs shared across data sets; make sure you have matching SNP ID formats across sumstats / reference data sets`. The diagnostic TSV contains all 22 rows: 10 `TESTED`, 11 `NOT_RUN`, and 1 `FAILED`. The R worker wrote those rows and its summary, then exited 1 as predeclared. No v2 output was overwritten, no workers 1/3/4 were replayed, and this result is not inserted into v2's missing worker-2 receipt.

The original canonical insomnia shard for locus 950 had **zero data rows**, as did the v2 native-N materialized shard. Their gzip hashes are identical: `572dbf5336c52498fc275cb5567ea69fd69e92b6d18994b8e2d9fdbd0954867b`. The canonical status was `NOT_RUN` with `FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS`. This is a structural empty-input condition in the existing filtered row set. There is no N value in this locus to replace, so its failure does not test the Jansen per-SNP N representation or its case-fraction approximation. The v2 runner classified LAVA's empty-input error as `FAILED`, triggering its frozen zero-failure hold. That rule remains final for v2; the postmortem cannot turn the cell into a v2 `NOT_RUN` after seeing outcomes.

The source remains `ADMISSIBLE_SENSITIVITY_ONLY` for technical consideration. This postmortem supplies no new confirmatory model or overlap evidence, no full-family outcome, and no protected candidate result.

## Receipts

V3 SSD root: `/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-worker2-postmortem-v3`

- Frozen v3 config SHA-256 `c6a3a913f3e237b357bfcf22dea925487a355979f3467f7d7d27019e556fa022`.
- Launch SHA-256 `670f06abadc0e1561e7c977e902c7eac4ae4d03250b8f697e6d1431196d54e68`.
- Diagnostic TSV SHA-256 `10e54d4808ed89c36664e7160ca7b8e4b38e805455afd1ea63ace6a0b5d4dff1`.
- Diagnostic summary SHA-256 `28b538dfaecbae537db62d19e69447efe8268ecaebe8211593c411e5b1829c68`.
- Diagnostic receipt SHA-256 `a023f33abfd63bbdcda1f781b206a61e0f3629672152775970dd099afcbd0457`; R exit code 1.
