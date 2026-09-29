# Frailty project checkpoint — 2026-09-27 19:40 UTC

The existing FI × sleep LAVA run remains active on the external SSD. Runner state at 19:40:12 UTC reports 20,636/29,940 receipts, three workers, zero launch failures and zero stale-claim recoveries. The rolling 30.1-minute throughput estimate remains 615.83 receipts/hour; the remaining-runtime estimate is provisional. The full read-only audit at 19:39:47–19:40:04 validated 20,633 receipts and saw three arrivals during scanning (20,636 inventory afterward). Receipt, claim, and duplicate-identity checks were clean; all three workers were live. All 12 locked 1% trait gates still fail (3,625 negative-variance, 252 no-reference-SNP, one other process failure). Evidence: `lava_full_receipt_audit_2026-09-27_1940.json` and `lava_runner_state_snapshot_2026-09-27_1940.json`.

Both distinct reviewers are ready. Their canonical archives still match the recorded SHA-256 values, each contains 113 batches, and both pass ZIP integrity checks. The title/abstract queue remains 56,117 records with zero decisions; no reviewer messages were sent and no decisions were applied. The plan, overlap, and review validators passed at 19:30 UTC; see `core_gate_validation_2026-09-27_1930.log`.

External licensed exports/full texts, independent decisions/adjudication, exact HFRS object access and exact cohort intersections remain outstanding. Current LAVA QC does not authorize local-sharing or downstream molecular analyses.

Figure 4 was regenerated from this read-only audit; the rendered panel was visually inspected and its evidence counts, failure wording and no-inference caveat match the audit. Renderer provenance records the SVG/PDF/PNG and source-data hashes.
