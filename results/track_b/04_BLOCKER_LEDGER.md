# Track B compute and data blocker ledger

This is an execution-readiness artifact, not an analysis result. It records terminal data limitations and the exact gates that prevent downstream stages from being run or interpreted.

Current machine snapshot:

- free storage: 244744192 bytes (0.228 GiB)
- physical memory visible to the process: 8589934592 bytes (8.000 GiB)
- MATLAB: `NONE`

Status counts:

- `BLOCKED_BY_COMPUTE`: 1
- `BLOCKED_BY_DATA`: 2
- `BLOCKED_BY_DATA_AND_COMPUTE`: 3
- `BLOCKED_UPSTREAM`: 9
- `NO_VALID_REPLICATION`: 1

The Pair A replication outcome is terminal `NO_VALID_REPLICATION`; it is not a negative replication result. Pair B has a frozen independent-cohort source but no result. No local correlation, pleiotropic locus, colocalization, tissue/cell mechanism, molecular-QTL, pathway, functional, causal, or final mechanistic claim may be made from this ledger.

After storage or runtime conditions change, rerun `python3 scripts/117_build_track_b_blocker_ledger.py` before executing the next stage. Empty or synthetic downstream tables are forbidden as substitutes for real analyses.
