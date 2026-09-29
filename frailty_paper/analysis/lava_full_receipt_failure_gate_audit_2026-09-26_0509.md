# Full LAVA receipt and locked failure-gate audit — 05:09 UTC

Audit snapshot: 2026-09-26T05:09:09.908871+00:00. The coordinator advanced from 9,844 to 9,845 while the read-only scan ran.

## Receipt integrity

All 9,845 receipt JSON files present at scan time passed `frailty_paper/scripts/45_run_lava_sensitivity.py::verify_receipt`, matched their frozen trait/pair/locus identity, and had a corresponding nonempty log. Receipt-file count reconciled with the coordinator's end state. Structural issues: **0**; duplicate trait/locus identities: **0**. Process statuses: 7,938 `PROCESSED` and 1,907 `PROCESS_FAILED`. The ordered receipt hash stream is `66ace45740fd4e8788a7f5890366c6a205fe9122ee920ad6bd706740f8e4f321` (frozen trait order, ascending locus index; each relative path, NUL, file contents, NUL).

The coordinator remained active with two workers and distinct claims: sleep_efficiency locus 562 (PID 78548) and sleepiness locus 987 (PID 78549). Launch failures and stale recoveries were both zero. Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. This audit did not interrupt workers, rerun loci, change settings, or collate/overwrite analysis results.

## Locked 1% locus-failure gate

The frozen collator's 1% locus-failure threshold is exceeded for every trait based on `PROCESS_FAILED` receipts alone. This remains a lower bound because the collator also counts non-tested FI/sleep univariate statuses.

| Sleep trait | Process failures | Lower-bound fraction | Negative-variance outcomes | No SNP IDs in reference | Passes 1% gate? |
|---|---:|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 23.53% | 566 | 21 | No |
| chronotype | 304 | 12.18% | 283 | 21 | No |
| insomnia | 80 | 3.21% | 60 | 20 | No |
| longsleep | 118 | 4.73% | 98 | 20 | No |
| napping | 86 | 3.45% | 66 | 20 | No |
| shortsleep | 105 | 4.21% | 85 | 20 | No |
| sleep_apnea | 46 | 1.84% | 26 | 20 | No |
| sleep_efficiency | 145 | 5.81% | 125 | 20 | No |
| sleep_timing | 146 | 5.85% | 126 | 20 | No |
| sleepdur | 85 | 3.41% | 65 | 20 | No |
| sleepiness | 110 | 4.41% | 90 | 20 | No |
| snoring | 95 | 3.81% | 75 | 20 | No |

Across the audited family, 1,665 logs reported the recognized all-phenotype negative-variance outcome; 242 reported that none of the specified SNP IDs were present in the reference data. Those 242 outcomes map to loci 950–969 across the 12 traits (chr6:25,684,630–33,864,262, spanning the MHC region) and locus 1484 for two traits (chr9:140,097,760–141,146,682). These remain recorded process failures; this audit did not establish why the reference intersection is empty and did not reclassify or repair them.

The campaign should continue under the existing atomic-claim runner to preserve planned receipts. Local genetic-sharing and downstream PLACO, fine-mapping, colocalization, or molecular claims remain unjustified under the frozen QC gate. Detailed counts and provenance are in `lava_receipt_full_audit_2026-09-26_0509.json`.
