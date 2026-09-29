# Full LAVA receipt and locked failure-gate audit

Audit snapshot: 2026-09-26T05:00:05.147599+00:00

## Receipt integrity

The coordinator state reported 9,819/29,940 receipts at the start and end of the audit. All 9,819 existing receipt JSON files passed the locked runner's receipt validator, matched the expected trait, pair-order, locus-index, and frozen locus coordinates, and had a corresponding nonempty log. The receipt-file count after the audit was 9,819; receipt JSON hash-stream SHA-256: `9888d7f0e5295fc738dc63e66e1c579ad5bd90b2ef2e5d869ec63ae98153aa6f`. Structural issues: 0; duplicate `(trait,locus)` identities: 0. The family remains incomplete at 9,819/29,940; 20,121 scheduled receipts remain.

The current coordinator state listed two distinct active claims (shortsleep locus 561 (PID 78548), sleep_apnea locus 986 (PID 78549)), zero launch failures, and zero stale recoveries. The frozen analysis-lock and input-manifest hashes remain unchanged. This audit did not stop workers, rerun loci, change settings, or collate/overwrite analysis results.

## Locked 1% locus-failure gate

The frozen collator's 1% locus-failure limit is exceeded for every one of the 12 sleep traits based on `PROCESS_FAILED` receipts alone. Its full locus-failure predicate also counts non-tested FI/sleep univariate statuses, so actual gate-failure fractions can be higher than these lower bounds. No threshold was changed.

| Sleep trait | Process-failed receipts | Lower-bound failure fraction | Negative-variance outcomes | No SNP IDs in reference | Passes 1% gate from process failures alone? |
|---|---:|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 23.53% | 566 | 21 | False |
| chronotype | 304 | 12.18% | 283 | 21 | False |
| insomnia | 79 | 3.17% | 59 | 20 | False |
| longsleep | 117 | 4.69% | 97 | 20 | False |
| napping | 86 | 3.45% | 66 | 20 | False |
| shortsleep | 103 | 4.13% | 83 | 20 | False |
| sleep_apnea | 46 | 1.84% | 26 | 20 | False |
| sleep_efficiency | 144 | 5.77% | 124 | 20 | False |
| sleep_timing | 145 | 5.81% | 125 | 20 | False |
| sleepdur | 85 | 3.41% | 65 | 20 | False |
| sleepiness | 109 | 4.37% | 89 | 20 | False |
| snoring | 94 | 3.77% | 74 | 20 | False |

Across the 9,819 receipts, 1,899 were marked `PROCESS_FAILED`: 1,657 logs contain the recognized “Negative variance estimate for all phenotypes” outcome, while 242 logs report “none of specified SNP IDs are present in reference data.” Those are valid recorded locus outcomes, not launch failures. The 242 no-reference-SNP cases map to loci 950–969 (all 12 traits; chr6:25,684,630–33,864,262, spanning the MHC region) and locus 1484 for two traits (chr9:140,097,760–141,146,682). This pattern is reported, not repaired or reclassified.

## Interpretation

The full family should continue to completion under the existing atomic-claim runner so the planned receipts and their QC diagnostics are preserved. Local genetic sharing and downstream PLACO, fine-mapping, colocalization, or molecular claims remain unjustified under the frozen failure gate. The first-pair and chronotype gates remain failed.

Detailed machine-readable counts and hashes: `frailty_paper/analysis/lava_receipt_full_audit_2026-09-26_0458.json`. Validator: `frailty_paper/scripts/45_run_lava_sensitivity.py::verify_receipt`; frozen failure rule: `frailty_paper/scripts/46_collate_lava_sensitivity.py` (`MAX_LOCUS_FAILURE=0.01`).
