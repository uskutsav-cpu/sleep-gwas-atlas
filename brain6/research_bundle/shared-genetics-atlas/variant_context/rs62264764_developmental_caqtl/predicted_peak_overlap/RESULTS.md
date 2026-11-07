# AlphaGenome prediction over the developmental ATAC peak

The saved adult frontal-cortex DNase prediction decreases for the alternate allele in both the portion of the 501-base window that overlaps the developmental ATAC peak and the portion outside it. Under the stated coordinate assumption, the peak-overlapping portion accounts for **−0.778** of the full window's net **−2.685** predicted ALT-minus-REF sum, about **29%** of that net difference. This is a post-hoc decomposition of model output, not the fraction of a biological effect caused by the peak.

| Segment | Bases | Predicted REF sum | Predicted ALT sum | ALT − REF | Log₂ ratio with offset 1 |
|---|---:|---:|---:|---:|---:|
| Complete 501-base window | 501 | 14.551 | 11.867 | −2.685 | −0.2734 |
| Developmental peak overlap within that window | 153 | 5.301 | 4.523 | −0.778 | −0.1901 |
| Remaining window outside the peak | 348 | 9.250 | 7.343 | −1.906 | −0.2969 |
| Entire developmental peak | 449 | 14.110 | 12.611 | −1.498 | −0.1507 |

The two disjoint raw differences add to the full-window difference. Their log ratios do **not** add: each segment uses its own offset of 1 and the segments have different lengths. The 29% is a ratio of net model-signal sums and can be affected by cancellation within segments, the chosen window and the peak definition. It is not a probability, risk effect, statistical significance test, or attribution of causality. Shifting both peak endpoints by one base in either direction changes the overlap-segment log ratio only slightly, from −0.1901 to approximately −0.1900 or −0.1903; this does not resolve larger annotation or cell-context differences.

The source is the previously saved, user-authorized JBrowse AlphaGenome adult frontal-cortex DNase track for chr3:117920728 G>A (GRCh38); the model release was not exposed. The separate [measured ATAC peak](../RESULTS.md) comes from primary developmental progenitors and differentiated neurons. We assume its peak identifier is zero-based and half-open; the source convention is not independently confirmed. The two rs62264764/rs62264767 markers had no records in the selected caQTL files. Thus this is **not a measured allele-specific accessibility effect**, proof the adult prediction corresponds to the developmental peak, a gene assignment, or evidence of sleep–ADHD causation. The strongest next step is allele-specific measurement and perturbation in matched, relevant cells with replicated donors.

The [protocol](PROTOCOL.md), [segment table](segment_scores.tsv), [calculation](analyze.py), [receipt](receipt.json) and [independent numerical verification](verification.json) record the coordinate mapping and calculations. The 501-base output reproduces the earlier saved score, raw segment deltas sum exactly, and an independent NumPy implementation reproduces all segment sums. These checks establish arithmetic and source-file consistency, not biological validity.
