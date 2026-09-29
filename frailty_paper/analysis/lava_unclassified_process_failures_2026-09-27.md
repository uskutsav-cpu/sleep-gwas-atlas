# Two unclassified LAVA process failures — 2026-09-27

The 20:01 UTC and 20:06 UTC read-only full-family audits report two process failures outside the two previously recognized categories. Read-only inspection of the corresponding receipt/log pairs identifies these records:

| Trait / locus | Attempt in worker log (UTC) | Receipt error | Log evidence | Receipt SHA-256 | Log SHA-256 |
|---|---|---|---|---|---|
| `longsleep` / 1545 | 2026-09-26 18:54:53 | `No write permission for directory: /var/folders/f4/8spnm17n3g9cvbpydqnw8m380000gn/T/RtmpVqt24S` | R warned about `writeBin` connection writes and `file.access`/append access to that temp path; no additional explicit error line. | `8dbe4e8270308accd789bd743b5aca6926e0dd08df1f96d95416b15226f70fe7` | `312fecb31bfbcc24b471baa347d4bfd28e6c813e0c6ef6a6ec93048f4db18716` |
| `sleep_apnea` / 1997 | 2026-09-26 19:57:50 | `process.locus returned NULL` | Last error line: `Error: invalid 'type' (complex) of argument`; this followed extraction/alignment of 179,137 shared SNPs. | `3e31a536d258acbb221cff9569f72884b1d98df8d2edeeaaa15b7ed338054ae5` | `821c61b68e7bb8472b6904bf37815728361c3bdfcd552881c46f8f877323b999` |

These are process failures, not receipt-integrity failures. The full-family auditor validates both receipt/log pairs; its receipt/claim/duplicate checks remain clean. The current run continues at three workers. Both slots count against the frozen trait-level failure gates, which all 12 traits still fail. This inspection did not remove, rewrite, or retry either receipt, change the frozen family, or modify the live run. The two stored messages do not by themselves establish a root cause or justify a retry.
