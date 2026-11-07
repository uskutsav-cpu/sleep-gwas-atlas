# Regional AlphaGenome peaks versus shared genetic evidence

These are the largest absolute predicted adult frontal-cortex DNase changes among the already selected variants **within each region**. The prediction is ALT versus REF over a 501-base window with offset 1; it is not a measured effect, gene-level expression, disease-risk direction, p-value or calibrated causal probability. The region-level PP4 is a single-signal trait–trait colocalization posterior, not a posterior for the named peak variant.

| Region (GRCh37) | Peak variant (GRCh38) | Predicted change | In insomnia / ADHD selected sets? | Regional PP4 at p12=10⁻⁵ / 10⁻⁶ | Positional gene context |
|---|---|---:|---|---:|---|
| chr11:28076261-29076261 | rs10835373; chr11:28622395:T>C | +0.092 | no / yes | 0.00288 / 0.000288 | Nearest LINC02758 (56788 bp) |
| chr14:98137479-99137479 | rs79470265; chr14:98182792:A>C | +0.085 | no / yes | 0.545 / 0.107 | Overlaps ENSG00000259097 |
| chr3:117142005-118142005 | rs62264764; chr3:117920728:G>A | -0.273 | yes / yes | 0.921 / 0.537 | Overlaps LINC03051, LSAMP |
| chr5:103447968-104481726 | rs2447828; chr5:104611504:C>T | -0.059 | yes / no | 0.991 / 0.917 | Overlaps ENSG00000251574 |
| chr5:87185500-88185500 | rs324886; chr5:88600784:C>T | -0.014 | no / yes | 0.138 / 0.0157 | Overlaps MIR9-2HG |
| chr7:113559156-114571035 | rs2189015; chr7:114508188:A>G | +0.050 | yes / no | 0.843 / 0.35 | Overlaps FOXP2 |

**Interpretation.** rs62264764 has the largest absolute prediction across the selected 425 variants and is in both admitted sets, but its regional PP4 falls under the lower sharing prior; neither establishes LSAMP or LINC03051 as the effector gene. The chr5:103 Mb region has the highest PP4 of these regions under both tested priors, but its peak-prediction variant rs2447828 is not in the admitted ADHD set; region-level sharing cannot be transferred to that particular marker. The chr11 and chr5:87 Mb regional models provide little H4 support at these priors despite having a predicted peak. A “no” means the variant is outside the project’s selected, admitted credible set, not that it has no association or cannot be causal.

Genes are overlap or nearest-gene annotations from the recorded Ensembl queries, not experimentally validated targets; numeric gene identifiers are kept when no gene symbol was returned. Gene boundaries, selected sets and 1-Mb windows are model- and release-dependent. No independent two-trait locus replication was verified, the original LAVA gates failed, and all entries are hypotheses for follow-up. Cross-variant prediction ranking is sensitive to window and offset choices; the JBrowse model release was not exposed.

Sources: [six-region locus atlas](../LOCUS_ATLAS.md), [prediction ranking](../../prediction_rank_comparison/REPORT.md), [window sensitivity](../../alphagenome_window_sensitivity/RESULTS.md), and [LSAMP evidence review](../../variant_context/rs62264764_lsamp_prior/RESULTS.md). [Code](build.py) and [input receipt](receipt.json) reproduce this selected-view table. The user-authorized AlphaGenome outputs are used for non-commercial research, not model training.
