# Storage and SSD execution

Status: **COMPLETE_WITH_QUALIFICATIONS** for storage recovery and operational launch. Scientific campaign execution remains active. Latest resource observation: 2026-10-09T19:56:06.358040+00:00.

The physical Extreme 55AE is a4 TB USB solid-state drive, mounted read/write at `/Volumes/Extreme SSD`, device `/dev/disk4s1`, whole disk4, ExFAT, volume UUID `77FD98CC-B09E-3DAA-ACC0-82FF2B776B28`. Diskutil reports SMART Verified. A65,536-byte write/fsync/readback probe passed its exact hash check. No unmounted full filesystem check is claimed. ExFAT metadata semantics are qualified separately from byte identity.

Initial internal free storage was604,614,656 bytes, below the unchanged3 GiB guard. The user explicitly authorized **“Approve this exact relocation”** after seven files had verified byte-identical SSD copies. Exactly4,355,618,590 logical bytes in the protected shared checkout's `data/raw.local-preserved/.archives` were relocated; that directory now resolves through an absolute symlink to the versioned SSD copy. All original hashes remain preserved. Git status SHA was identical before/after, and no original open file handle remained at the switch. Internal free storage increased from563,666,944 to4,950,347,776 bytes. No other change to the frailty/Brain6 checkout is authorized by this operation.

At the latest observation, free internal storage is4,528,574,464 bytes (4.218 GiB), and SSD free is1,035,289,493,504 bytes (964.188 GiB). These are dated observations, not promises of future capacity. RAM is8 GiB. Current swap report:

```
hw.memsize: 8589934592
vm.swapusage: total = 3328.00M  used = 2945.19M  free = 382.81M  (encrypted)
```

The unchanged original190-job campaign writes support files, native captures, logs and block vectors under `/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/native/original_runner_support_v4_2`. TMPDIR and Python/research caches are on the SSD. Essential macOS files and virtual-memory swap are unchanged. One heavy scientific worker and one BLAS thread are allowed; owned RSS and storage are sampled every2 seconds. The original3 GiB internal and5 GiB SSD floors,2 GiB observed owned RSS cap,36hour stage deadline and4 GiB native-output cap remain. Resource breaches stop and preserve the affected stage for review. Observed RSS sampling does not guarantee an instantaneous peak limit.

Raw public source transfers are separately scoped and sequential, with hard per-source curl byte limits and final original checksum checks. The original100 extension bodies total227,610,388,647 bytes (211.979 GiB), with a frozen300 GiB SSD reservation. The corrected acquisition uses an exclusive family lock and immutable attempt/receipt paths. Its original interrupted partial and every superseded plan/failure remain preserved. No automatic retry is permitted. The single new FinnGen R13 insomnia body is separately admitted only for qualified source feasibility.

Evidence: `manifests/ssd_operational_preflight_v4.json`, `logs/internal_archive_relocation_candidate_receipt_v4.json`, `logs/authorized_archive_relocation_receipt_v4.json`, `logs/storage_state_after_native_core_v4.json`, `manifests/ssd_native_execution_plan_v4_3.json`, `logs/core_native_monitor_receipt_v4.json`, and the versioned acquisition reviews/admission. The exact relocation rollback copies into a new internal sibling, verifies all seven hashes, then replaces only the symlink; it retains the SSD originals and requires adequate internal capacity.
