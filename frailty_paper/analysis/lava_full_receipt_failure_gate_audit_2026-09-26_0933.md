# FI×sleep LAVA full-family audit — 2026-09-26 09:33 UTC

The read-only full-family auditor validated 11,875/29,940 receipt files from 09:33:16 to 09:33:21 UTC. The inventory grew by one during the scan; the coordinator reached 11,878 at 09:33:31. Remaining: 18,062.

Receipt integrity, duplicate claim identities, duplicate trait/locus identities, and claim issues were zero. Workers 78548 and 78549 were live on distinct claims. Launch failures and stale recoveries were zero. Six workers remain requested; two are effective under the existing adaptive safeguard, whose persisted trigger is 2,772/3,072 MiB swap (90%). Lock hash is `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; input-manifest hash is `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

From the 09:25 audit (11,783 receipts at 09:25:22) to this audit (11,875 at 09:33:16), observed throughput is about 700 receipts/hour; approximately 26 hours remain at that rate. This is a provisional measured estimate. The current all-traits gate remains failed: 2,025 all-phenotype negative-variance outcomes and 242 no-reference-SNP outcomes; no local-sharing inference is admissible.

## Host resource sample — 09:28 UTC

The authorized read-only sample found 8 logical CPU cores; load averages 6.51/6.19/6.06; CPU use 60.84% user, 23.11% system, 16.3% idle. `memory_pressure` reported 50% free; `top` reported 7,481 MiB physical memory used, 152 MiB unused, 1,435 MiB wired, and 2,164 MiB compressed. Swap was 1,519.31/3,072 MiB used (1,552.69 MiB free). At the process sample, worker 2's R job used 108% CPU and 393,360 KiB RSS; worker 1's R job used 5% CPU and 130,432 KiB RSS. The first `iostat` row showed 100.36 MB/s on internal disk and 32.23 MB/s on the external disk; subsequent one-second intervals fell below 0.5 MB/s on each. Internal free disk was 6.5 GiB and external free disk was 1.6 TiB.

The worker count was not raised: swap was only about 17 MiB below the scheduler's 1,536 MiB downshift rung, CPU idle was 16.3%, and prior campaign observations had already triggered the 6→4→3→2 memory safeguards. The current worker claims were left undisturbed. This is a resource-conservative continuation of the existing run, not a restart.

The locked plan validator passes (12 traits, multiplicity 396); the cohort ledger validates structurally (34 rows), while exact participant intersections remain unknown. The reviewer queue remains at 56,117 with zero decisions; both 113-batch packet sets pass dry-run import without decisions applied. Human screening/adjudication, licensed database exports, exact HFRS summary-statistic identity/access, exact overlap, and independent replication remain unresolved. Overall readiness remains **NO-GO**.
