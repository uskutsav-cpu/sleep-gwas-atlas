# FI×sleep LAVA full-family audit — 2026-09-26 09:42 UTC

The read-only full-family auditor validated 11,988/29,940 receipt files from 09:41:58 to 09:42:03 UTC. The inventory grew from 11,988 to 11,990 during scanning; the coordinator reported 11,991 verified receipts at 09:42:07. Remaining: 17,949.

Both worker PIDs passed liveness checks on distinct active jobs. Receipt integrity, claim issues, duplicate identities, and duplicate claims were zero; launch failures and stale recoveries were zero. Six workers remain requested and two effective. The persisted safeguard records a prior trigger at 2,772/3,072 MiB swap (90%), reducing three workers to two. Frozen lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; input-manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

From the 09:36 scan (11,924 validated at 09:36:58) to this scan (11,988 validated at 09:41:58), measured throughput is approximately 768 receipts/hour. About 17,949 remain, corresponding to roughly 23.4 hours at this short-window rate; this estimate is provisional. Failure counts are 2,042 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 locked 1% gates fail; no local-sharing inference is admissible.

## Resource check

The authorized 09:41 host sample found 8 CPU cores, load averages 6.21/6.33/6.22, and CPU utilization 37.93% user / 13.79% system / 48.27% idle. `memory_pressure` reported 53% free. `top` reported 7,412 MiB memory used, 218 MiB unused, 1,397 MiB wired, and 2,146 MiB compressed. Swap remained 1,519.31/3,072 MiB used. The active R child for worker 2 (`sleepdur`, locus 1096 at sample time) used 549.2% CPU and 914,576 KiB RSS; worker 1 (`sleep_apnea`, locus 667) used 130.7% CPU and 436,928 KiB RSS. The next one-second disk interval was 0.06 MB/s internal and 0.09 MB/s external after the initial cumulative rates of 99.54/31.95 MB/s. The external SSD retains 1.6 TiB free; internal storage is at 6.5 GiB free.

The active R jobs together occupied about 1.3 GiB RSS, while swap was only 16.7 MiB under the 1,536 MiB downshift rung. I left the current two-worker run undisturbed; adding workers would risk immediate fallback and more swapping. The previous measured two-worker rate is materially higher than earlier four-worker intervals.

The 12-trait/396-test analysis plan passes validation. The 34-row cohort ledger is structurally valid but exact participant intersections remain unknown. Review validation passes structurally, but both reviewers have returned zero decisions; licensed database exports, adjudication, exact HFRS file identity/access, and independent replication remain unresolved. Overall readiness remains **NO-GO**.
