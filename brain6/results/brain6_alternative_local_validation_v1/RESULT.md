# Brain6 secondary local covariance result

The frozen five-pair family contains 8,465 non-MHC LDetect block slots
(5 pairs × 1,693 EUR GRCh37 blocks). The Bonferroni threshold remained
`0.05 / 8,465 = 5.906674542232723e-06` throughout. Two source-eligible
insomnia pairs ran with four workers each, sequentially; all 1,693 block slots
per eligible pair are receipt-bound. Three long-sleep pairs were prospectively
marked `METHOD_INAPPLICABLE` under `source_eligibility_lock_v1.json` and were
never run, because the exact Dashti export's model and per-SNP analyzed-N
mapping remain unresolved. Their rows contain no estimated p value.

| Pair | Source eligibility | Estimated blocks | Low-SNP | Low-rank | FWER-significant blocks |
|---|---|---:|---:|---:|---:|
| insomnia–ADHD | ELIGIBLE_SECONDARY | 1,646 | 38 | 9 | 1 |
| insomnia–MDD | ELIGIBLE_SECONDARY | 1,658 | 23 | 12 | 2 |
| long sleep–SCZ | METHOD_INAPPLICABLE | 0 | 0 | 0 | 0 |
| long sleep–bipolar | METHOD_INAPPLICABLE | 0 | 0 | 0 | 0 |
| long sleep–Parkinson | METHOD_INAPPLICABLE | 0 | 0 | 0 | 0 |

The three FWER-significant genome-wide blocks do not overlap any of the
protected pair-specific candidate intervals. Of the 25 protected candidates,
12 were evaluated and `NOT_SUPPORTED`; 13 are `METHOD_INAPPLICABLE`. Of the 20
geographic regions, 10 are `NOT_SUPPORTED` and 10 are `METHOD_INAPPLICABLE`.
No protected candidate receives secondary local covariance support.

The separate receipt-bound `candidate_direction_audit.tsv` compares the best
estimated block's covariance sign with reference-aligned lead-SNP Z products
where both lead SNPs survived the method input QC. Four eligible candidates
have a concordant comparison; eight eligible lead SNP comparisons are
unavailable after reference/input matching. The 13 inapplicable candidates
have no local direction comparison. Sign concordance is descriptive and does
not change support status.

These are same-source methodological checks, not independent cohort
replication. The scalar study-total N approximation and meta-analysis
heterogeneity remain limitations even for the two eligible pairs. These
results do not alter the frozen primary LAVA failure (13,745 `TESTED`, 3,720
`NOT_RUN`, 0 `FAILED`; ceiling 873) or promote any protected candidate.

`run_provenance.json` binds the frozen protocol, source eligibility lock, two
pair receipts, and SHA-256 hashes of all three tabular outputs. The complete
block-level rows are in `genome_wide_blocks.tsv`; the pair-specific and region
joins are in `candidate_crosswalk.tsv` and `region_crosswalk.tsv`.
