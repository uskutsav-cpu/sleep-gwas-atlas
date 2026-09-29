# Brain6 PLACO+ family v1

This is the frozen, pre-result production plan for the four new PLACO+ follow-up pairs. It excludes `insomnia__adhd`, whose Track B output is protected and must not be rerun. Pair selection is inherited from `extensions/brain6/work/overnight-v03/brain6-pairs.lock.json`; the planned family denominator remains **five selected tracks**, including protected Track B, as required by the frozen correction policy.

## Frozen test and QC policy

- Method: official PLACO 0.2.0 `PLACO_PLUS`, source pinned in `brain6/manifests/placo_source_v0.2.0.json`.
- Pair inputs: checksum-verified prepared genome-wide pair joins and deterministic 20,000-row chunks under `/Volumes/Extreme SSD/brain6-work/preparation-v1`.
- Inherited variant QC: every selected dense source is accompanied by a source-processing QC log that records the CDG3 `invalid FRQ or MAF <= 0.01` filter. The source-level output row counts match the current locked dense inputs for all seven traits; the upstream-QC SHA256s, source SHA256s, and row-count checks are frozen in `brain6/qc/upstream_frequency_filter_audit.tsv`. This establishes inherited frequency filtering for the present pair inputs; it does not independently recover the historical source transformation code or frequency orientation.
- PLACO+ global null parameter estimation uses genome-wide variants with `Z1² <= 80` and `Z2² <= 80`, with the pinned implementation's `p_threshold=1e-4`, minimum 1,000,000 rows, seed 20260908, and absolute tolerance 1e-13. Any pair failing parameter estimation or minimum-count gates remains failed; do not lower gates.
- Chunk-level results are collated once per pair. Within-pair BH uses the full row denominator, with numerical failures conservatively retained as p=1 for denominator accounting and separately reported. Results exceeding the predeclared 0.1% failure-rate ceiling fail QC and are not consumed.
- Across-pair correction is Bonferroni across **five** selected tracks: headline threshold `5e-8 / 5 = 1e-8`; `Q_BONFERRONI_ACROSS_PAIRS = min(1, Q_WITHIN_PAIR * 5)`. Track B is retained in the correction denominator although it is not rerun here.
- Independent-locus definition and evidence tier are applied only after a passing pair-level result and are not inferred from raw variant hits.

The four `*/plan.json` files bind each plan to its pair lock, pair file, chunk manifest, and pinned PLACO source. Their own `plan_sha256` guards the settings against edits. The first two v1 chunks for insomnia–MDD are preserved under `native_results/`; all observed numerical failures were exact-zero-Z cases and exceeded the frozen 0.1% pair failure ceiling. The v1 execution was stopped before collation or interpretation. See `brain6/manifests/placo_family_v1_run.json`. The separately frozen v2 adapter applies the exact mathematical P=1 boundary and retains all other settings and QC thresholds.

## Pair denominators

| Pair | Prepared variant rows | Chunks |
|---|---:|---:|
| insomnia–MDD | 1,122,359 | 57 |
| long sleep–SCZ | 6,296,500 | 315 |
| long sleep–bipolar disorder | 6,305,673 | 316 |
| long sleep–Parkinson's disease | 1,131,725 | 57 |

Across the four unprotected pairs, the prepared inputs contain 14,856,257 pair-variant rows. Each pair has its own full-variant denominator for BH; this sum is not one pooled BH family.

## Interpretation limits

No row-wise `INFO` field is present in the Brain6 dense inputs. Upstream QC applies `INFO > 0.6` where the source provides it, and records INFO absent for MDD and Parkinson's disease; this plan does not claim uniform INFO quality. PLACO+ sharing is statistical cross-trait evidence, not proof of biological pleiotropy or causality. Frequency orientation is not needed for the PLACO+ Z-statistic calculation and is explicitly not asserted by the filter-lineage artifact.
