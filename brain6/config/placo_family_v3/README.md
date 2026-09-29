# Brain6 PLACO+ family v3

V3 is the frozen production rerun after preserving two earlier execution attempts. V1 exposed a genuine PLACO+ numerical boundary case: rows with `Z1 == 0` or `Z2 == 0` have a product statistic of zero and exact tail probability 1, but the official integral can return an invalid value at its integrable singularity. The adapter now handles only that exact boundary analytically and calls the unchanged, hash-pinned upstream implementation for every nonzero pair. A focused R test confirmed the boundary outputs and exact parity with the official function away from zero. The frozen numerical-failure ceiling remains 0.1%.

V2 parameter estimates passed, but the transactional artifact system rejected active chunks because the bundled extension code tree changed during execution when the adapter copy was synchronized. Those attempts remain preserved in the v2 run roots; no v2 chunk output was promoted or interpreted. V3 starts from separate scratch roots and records a complete implementation file-hash inventory in `family_lock.json` so code stays fixed throughout the run.

## Frozen analysis policy

- Selection: four new pairs from the reviewed five-track lock; protected Track B (`insomnia__adhd`) is not rerun and remains included in across-track correction.
- Method: official PLACO 0.2.0 `PLACO_PLUS`, source hash in the family lock.
- Global null estimation: all genome-wide rows with `Z1² <= 80` and `Z2² <= 80`; `p_threshold=1e-4`; minimum 1,000,000 variants; seed 20260908; absolute tolerance 1e-13.
- Pair denominators and deterministic 20,000-row chunks: insomnia–MDD 1,122,359 (57); long sleep–SCZ 6,296,500 (315); long sleep–bipolar disorder 6,305,673 (316); long sleep–Parkinson's disease 1,131,725 (57).
- Pair-level QC: full-row denominator; within-pair BH; numerical failures are counted conservatively in the denominator and remain failed; pair fails if numerical failure rate exceeds 0.1%.
- Across-pair correction: Bonferroni across five selected tracks, including protected Track B; headline threshold `5e-8/5 = 1e-8`.
- Inputs inherit the checksum- and row-count-audited upstream MAF filtering in `qc/upstream_frequency_filter_audit.tsv`. INFO availability remains nonuniform: no row-wise INFO column is present in the dense inputs, and upstream logs record INFO absent for MDD and Parkinson's disease.
- Outputs write to per-pair paths on `/Volumes/Extreme SSD/brain6-work/placo-family-v3-results/`; symlink targets are recorded in the family lock and run manifest.

All plan hashes, code inventory, selected-pair lock, source code, runtime, and scratch roots were frozen before v3 results. No v3 results have been reviewed at freeze time. PLACO+ is statistical cross-trait evidence and does not establish causality.
