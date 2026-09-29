# FI×sleep LAVA full-family audit — 2026-09-26 09:17 UTC

The read-only auditor validated 11,685 of 29,940 immutable receipt files (scan 09:16:47–09:16:53 UTC); inventory grew by one while scanning and the coordinator advanced to 11,688 by 09:17:02. The coordinator snapshot has 18,252 receipts remaining.

Receipt integrity, claim identity, and duplicate-job checks report zero issues. Two workers are active on distinct loci (PIDs 78548 and 78549); launch failures and stale recoveries are zero. Six are requested and two effective under the existing memory safeguard, triggered at 2,772/3,072 MiB swap (90%). The run remains live and resumable, with no scientific-input or configuration changes.

From the 09:06 full audit (11,565 validated) to this scan (11,685), 120 receipts were validated in about 10.1 minutes, approximately 715/hour; 18,252 remaining would take about 25.5 hours at that observed rate. This is a provisional estimate. The recent shorter 09:11–09:17 interval was slower, so throughput will continue to be monitored. Process-level CPU/RAM and disk-I/O readings remain unavailable because sandbox access to `ps`, `sysctl`, and `iostat` is denied.

The frozen analysis-plan validator passes (`version=1`, 12 sleep traits, multiplicity 396). The 34-row cohort ledger validates and frozen source IDs match, but exact cohort intersections are still unknown. The 1% LAVA failure gate remains failed for all 12 traits: 2,003 all-phenotype negative-variance process outcomes and 242 no-reference-SNP outcomes. No local-sharing inference is promoted; overall readiness remains **NO-GO** pending analysis gates, human review, licensed-source imports, exact HFRS summary-statistic access, cohort-overlap resolution, and independent replication.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Frozen input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0917.json`.
