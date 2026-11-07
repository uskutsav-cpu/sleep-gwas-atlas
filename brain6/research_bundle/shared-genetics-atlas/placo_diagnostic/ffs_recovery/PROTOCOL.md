# FFS PLACO+ rerun after gzip integrity failure

The original FFS worker exited zero but its `FFS_screen.tsv.gz` fails a CRC check. The combined audit also found an invalid status. Preserve that file as evidence; do not interpret its results or quietly repair its gzip checksum. FI's completed gzip passes the integrity check and is audited separately with the original two-pair family denominator.

Rerun all FFS rows, with the same pinned PLACO+ implementation, Z input, fitted parameters, exclusion rule, integration tolerance and output fields as the original protocol. Partition the ordered input into two non-overlapping contiguous ranges, executed independently; preserve completed chunks and logs on the external volume. Validate that the ranges cover all 1,138,262 rows exactly. Each chunk closes and atomically renames only after its range finishes. Verify each gzip CRC before streaming the two ordered chunks into a single new gzip file under this recovery folder, writing the common header once; do not replace or delete the damaged original. The combined output must itself pass gzip integrity and the full completed-scan audit, including row identities, counts and all control values, before it is used. Capture source, chunk and combined hashes. Interruption, an error, a missing chunk or insufficient storage is not success. This rerun does not resolve calibration, p-value precision, external sample overlap or LD limitations.

Once both chunk workers have exited successfully and both completed chunk files exist, run:

```sh
python3 outputs/shared-genetics-atlas/placo_diagnostic/ffs_recovery/merge_chunks.py
python3 outputs/shared-genetics-atlas/placo_diagnostic/audit_completed.py \
  --input '/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support' \
  --results outputs/shared-genetics-atlas/placo_diagnostic \
  --ffs-results outputs/shared-genetics-atlas/placo_diagnostic/ffs_recovery \
  --out outputs/shared-genetics-atlas/placo_diagnostic/completed_audit \
  --controls outputs/shared-genetics-atlas/placo_diagnostic/control_results.tsv
```

The merge checks each input gzip through end of file, validates its row count and identical headers, writes one header to a new gzip, rechecks that gzip's CRC, and retains all source files. It refuses to overwrite a partial or final output. `test_merge.py` checks header, row-count and termination behavior in memory. The subsequent audit remains necessary to compare every row and control against the original aligned inputs; a CRC only validates the compressed stream's integrity, not the correctness of its contents.
