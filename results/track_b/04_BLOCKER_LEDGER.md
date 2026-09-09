# Track B compute and data blocker ledger

This is an execution-readiness artifact that binds terminal global-replication classifications. It is not a local or mechanistic analysis result. It records terminal data limitations and the exact gates that prevent downstream stages from being run or interpreted.

Current machine snapshot:

- free storage: 16445267968 bytes (15.316 GiB)
- physical memory visible to the process: 8589934592 bytes (8.000 GiB)
- LAVA RAM assessment: `RAM_BENCHMARK_REQUIRED` from 0 measured discovery loci; maximum observed peak NA GiB plus a 1073741824 byte non-worker reserve
- MATLAB: `NONE`

Status counts:

- `BLOCKED_BY_DATA`: 2
- `BLOCKED_BY_DATA_AND_COMPUTE`: 2
- `BLOCKED_UPSTREAM`: 9
- `COMPLETE`: 1
- `NO_VALID_REPLICATION`: 1
- `READY_FOR_RAM_BENCHMARK`: 1

The Pair A replication outcome is terminal `NO_VALID_REPLICATION`; it is not a negative replication result. Pair B completed external-cohort LDSC with a significant concordant global genetic correlation (`rg=0.3817`, `SE=0.0521`, `P=2.32742962810806e-13`), earning `DIRECTIONAL_REPLICATION` with the frozen Finnish-founder, phenotype-equivalence, and unverified individual-overlap caveats. This clears only the global replication gate. No local correlation, pleiotropic locus, colocalization, tissue/cell mechanism, molecular-QTL, pathway, functional, causal, or final mechanistic result is present or promoted by this ledger.

After storage or runtime conditions change, rerun `python3 scripts/117_build_track_b_blocker_ledger.py` before executing the next stage. Empty or synthetic downstream tables are forbidden as substitutes for real analyses.
