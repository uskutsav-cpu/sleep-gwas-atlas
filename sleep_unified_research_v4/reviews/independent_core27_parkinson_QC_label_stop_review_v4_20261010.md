# Independent Parkinson QC-label stop adjudication

2026-10-10. **The actual scientific outputs and QC counts match; the preserved stop is caused solely by two QC-reason specificity labels.** Receipt status remains `FAILED_PRESERVED_NO_AUTOMATIC_RETRY`, with `COMPLETE_CORE_SERIALIZED_CONTENT_OR_QC_DIFFERS`. This review preserves that stop and supplies a prospective comparator recommendation; it grants no retry or resume authority.

Only the admitted frozen plan, the compact Parkinson comparison receipt and new `parkinson.qc.txt` were read. Historical QC is taken from the frozen plan, without reading historical/raw/output bodies, code, workers or other traits. Plan SHA-256: `8a839122ef7a9a0faefb6c3971f04abd8fa755ed630a6db1efe8f725ca54c6e3`.

The exact ordered-step differences are:

| Ordered step | Frozen reason | New reason | Dropped, both | Remaining, both |
|---:|---|---|---:|---:|
| 2 | `alleles conflict with pinned GRCh37 HapMap3 map` | `alleles conflict with pinned GRCh37 EUR HapMap3 map` | 0 | 1,182,445 |
| 3 | `source coordinate conflicts with pinned GRCh37 HapMap3 map` | `source coordinate conflicts with pinned GRCh37 EUR HapMap3 map` | 0 | 1,182,445 |

All **20** dropped/remaining pairs agree; the other **18** reason strings are identical. The already EUR-specific first mapping label is unchanged. Sequential arithmetic independently reconciles **17,510,617 input = 16,378,539 rejected + 1,132,078 harmonized retained rows**. The common nonzero losses are 16,328,172 absent from the pinned map, 1,749 frequency/MAF failures and 48,618 below the sample-size cutoff. These two zero-drop label changes have no retained-row or QC-count effect.

Every original scientific metadata field matches, including hg19 input/output build, `BY_COORD_ALLELES`, map SHA-256 `6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4`, raw SHA-256, per-variant case/control sample-size mode, row totals and 6.47% retention. Only `infile` and `variant_map` execution paths differ. Five additive metadata fields report the map scope/schema/size/provenance and absent prefilter; these are distinct from the historical metadata and were not newly certified against a map body. All 13 scalar checks in the executed receipt are true; only the literal `ordered_filter_steps` check is false.

The completed full-stream comparisons in that failed receipt establish **zero unequal harmonized rows**, identical ordered decompressed harmonized SHA-256 `b69764fa51cd330dc48c1ecad2f150a27bf6565a98bbea56813f61935240e84b`, and identical full HM3 SNP/allele/N/Z/missingness content with decompressed SHA-256 `3bdec94ca74364afcd0f69e75fbb1efed9aadcb09dcf62bacd0de9fcc06322f7`. Both comparisons completed gzip CRC/EOF validation. The **1,217,311-row HM3 template has 1,131,977 finite N/Z rows and 85,334 missing/nonfinite rows**, identical to the archive. No tolerance was introduced. Compressed hashes differ for both outputs; their cause is not inferred. This consumes executed stream evidence and does not repeat its body scan.

The result supports a **label-specificity discrepancy with no observed scientific-content/filter effect**, rather than a biological or numeric replay discrepancy. It does not independently prove historical code/runtime identity or remove source qualifications: INFO remains absent, per-SNP effective N remains the original NCASE/NCONTROL transformation, and the Nalls public release’s UKB proxy-case composition and overlap/construct limitations persist. No source-validity, fit, discovery or replication admission follows.

Recommendation: an existing-pipeline resume should first use a prospectively reviewed comparator correction that recognizes only these **two exact old/new reason aliases, at ordered positions 2 and 3**, under `parkinson`, `BY_COORD_ALLELES` and the pinned map identity above. Keep exact step number/order, zero dropped counts, 1,182,445 remaining counts, every other label/metadata/count comparison, full-stream SNP/allele/N/Z/missingness equality and CRC/EOF requirements. Record the accepted aliases explicitly while preserving both literal QC texts and the original failed receipt. Do not globally strip “EUR,” weaken numeric/QC equality, mutate outputs, or retroactively relabel the failed run as successful. A new plan/comparator identity and separately authorized continuation remain necessary; this adjudication launches nothing.

Evidence identities:

| Compact evidence | SHA-256 |
|---|---|
| Actual comparison receipt | `5777cca448d21b3eafcff99fd70d0fb56f2e76902963cdba1a4a401f10860dec` |
| New Parkinson QC | `dbc1d62c32b596e02be5691b7a3fcc875378998507c2edd81446467e27955815` |
| Frozen historical QC identity, inherited from plan | `1c797b936c87020036aa2728a9535407df0120786644fc6fe074cb466224b83d` |
