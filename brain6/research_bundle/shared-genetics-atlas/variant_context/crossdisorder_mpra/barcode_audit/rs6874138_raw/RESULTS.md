# rs6874138: raw-barcode and aggregate-count diagnostic

The released mapping associates 125 REF and 131 ALT barcode labels with rs6874138. Exact matched labels found in the 44 raw count files produce the following totals, compared with the released aggregate matrices. The raw-to-matrix ratios are sums over 22 files or columns; pairing of individual replicates below assumes tar order and is not independently established.

| Material and allele | Exact-matched raw reads | Aggregate-matrix reads | Ratio | Pearson r under assumed order |
| --- | ---: | ---: | ---: | ---: |
| DNA REF | 43,032 | 43,796 | 0.983 | 0.9995 |
| DNA ALT | 33,410 | 40,948 | 0.816 | 0.9952 |
| RNA REF | 49,231 | 50,243 | 0.980 | 0.9998 |
| RNA ALT | 42,272 | 76,640 | 0.552 | 0.6749 |

Under the assumed file order, every raw total is below its corresponding matrix count. Reference counts track the matrix closely, while the matched ALT subset accounts for less of the aggregate, especially for RNA. These patterns are consistent with an incomplete or differently processed barcode mapping, but **do not identify the cause**. A high correlation does not prove a replicate crosswalk, and the raw files, mapping and aggregate matrices may reflect different filtering or library processing. These partial totals are not a substitute for the authors' full MPRA model or a quantitative estimate of the variant effect. They also do not prove that the published effect or construct is wrong.

The results make the [conditional ALT-direction interpretation](../../published_table/allele_direction_audit/RESULTS.md) more qualified: its reconstruction matches the aggregate matrices, but this raw-barcode check has not independently reconstructed the ALT matrix or confirmed barcode/construct provenance. The other active reporter candidate, rs325485, had no exact raw-label matches in the prior 44-file scan. Independent validation needs the authors' complete barcode and sample crosswalk and preprocessing rules, or an independent assay.

The [pre-result protocol](PROTOCOL.md), [44-row table](counts.tsv), [scan code](audit.py), [source receipt](receipt.json), and [verification](verification.json) document the diagnostic. All four compressed source hashes match earlier project receipts, all 44 raw streams were read to completion, and summary counts were independently checked from the saved table. The matching uses a known selected variant and does not provide independent GWAS replication, a causal gene, a sleep-independent frailty signal, or an AlphaGenome model-training label.

A [GEO supplementary-directory inventory](geo_inventory.json) lists the raw archive, two count matrices, exclusion list and one `unique_bcdat` mapping; no separate sample-to-matrix or full barcode-processing crosswalk was identified in that listing. Its displayed modification timestamps differ between the matrices and mapping/archive, but timestamps do not establish source lineage or explain this discrepancy. The listing does not rule out a crosswalk in another repository or supplementary source.
