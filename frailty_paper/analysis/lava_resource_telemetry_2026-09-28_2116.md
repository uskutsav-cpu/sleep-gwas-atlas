# FI × sleep LAVA resource telemetry — 2026-09-28 21:16 UTC

Read-only host and runner sample while the original campaign remained active. No coordinator settings, claims, receipts, or scientific inputs were changed.

- Runner: PID 7390 alive; state updated at 21:15:38 UTC with 27,828/29,940 verified receipts, 6 requested / 2 effective workers, zero launch failures, and zero stale recoveries.
- Active jobs: longsleep locus 196 (R PID 5041; process group 37878) and napping locus 196 (R PID 5076; process group 89389). Both R processes were running at 98.8% CPU in the process sample. Their RSS was approximately 1.3 MiB and 1.2 MiB, respectively; mapped VSZ was about 412 GiB each.
- Receipts: latest append-only receipt remains sleepdur locus 205 at 19:36:32 UTC; 2,112 remain. The 21:00:19–21:01:00 full audit independently validated all 27,828 inventoried receipts with stable inventory, no receipt/claim/duplicate issues, and both active claims alive. The short-window completed-receipt rate is zero; the workers remain CPU-active, so this is not evidence that their claims are stale and no completion ETA is defensible.
- Host: 8 logical CPUs; `memory_pressure` reported 57% free. Swap was 6,663.25/8,192 MiB (81.3% used; 1,528.75 MiB free). Free pages were 4,922 at 16 KiB/page; 1,004,687 pages were stored in the compressor.
- Storage/I/O: `/Volumes/Extreme SSD` had 1.3 TiB free. Two `iostat -d 2 2` samples showed disk4 at 4.54 then 0.03 MB/s and disk0 at 23.81 then 0.52 MB/s. This is system-wide sampled I/O, not a process-specific measurement.
- Frozen configuration: lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both match the audited frozen inputs. All 12 locked 1% trait gates fail (4,800 all-phenotype negative-variance, 252 no-reference-SNP, four other failures).

The existing two-worker safeguard remains in force. No worker was restarted, no claim was recovered, and no QC threshold was changed.
