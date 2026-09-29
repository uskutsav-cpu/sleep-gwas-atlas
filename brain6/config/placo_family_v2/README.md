# Brain6 PLACO+ family v2

This immutable v2 execution plan corrects one identified numerical boundary case in the adapter. The v1 attempt is preserved under `config/placo_family_v1/`: its first two completed chunks show that every numerical failure occurred when either Z score was exactly zero. Since the PLACO statistic is `Z1 × Z2`, its two-sided tail probability at an observed statistic of zero is exactly 1. The official PLACO+ integral has an integrable singularity at zero and returned an invalid numerical value for that case. The v2 adapter assigns the exact mathematical value 1 only for this boundary and calls the unchanged, hash-pinned upstream PLACO+ implementation for every nonzero pair of Z scores. This fix does not change the frozen 0.1% pair failure-rate ceiling or any significance threshold.

The adapter fix passed `extensions/brain6/tests/test_placo_zero_boundary.R`: both zero-coordinate boundary cases return exactly 1, and a nonzero case matches the official function exactly at the frozen settings. The four plans below were generated after that test and before inspecting any v2 result. To preserve workspace space, each plan’s `native_results` directory resolves to the corresponding external scratch run root on `/Volumes/Extreme SSD/brain6-work/placo-family-v2-results/`; these paths are listed in the family lock.

## Frozen analysis policy

- Pinned PLACO 0.2.0 source and pair lock are identical to v1.
- The five-track correction, including protected Track B, remains frozen; the headline threshold is `5e-8/5 = 1e-8`.
- Pair-specific full-row denominators and chunk counts: insomnia–MDD 1,122,359 (57); long sleep–SCZ 6,296,500 (315); long sleep–bipolar disorder 6,305,673 (316); long sleep–Parkinson's disease 1,131,725 (57).
- Global null estimation retains `Z1² <= 80` and `Z2² <= 80`, `p_threshold=1e-4`, at least 1,000,000 rows, seed 20260908, and absolute tolerance 1e-13.
- Within-pair BH uses the complete pair denominator; numerical failures conservatively count as p=1 for denominator accounting, remain marked failed, and fail pair QC above 0.1%.
- Pair plans bind the same source, pair lock, pair joins, and chunk manifests as v1. The adapter hash in each native job receipt identifies the v2 zero-boundary fix.

## Limits

The upstream source transformation code is unavailable, but upstream filter lineage is checksum- and row-count-audited in `qc/upstream_frequency_filter_audit.tsv`. No row-wise INFO field exists in the current dense inputs; INFO availability remains nonuniform. PLACO+ results are cross-trait statistical evidence and do not establish causality.

Execution starts only after freezing this README and the four hashed plan files. No v2 results have been reviewed at freeze time.
