# Frailty project progress — 2026-09-27 19:30 UTC

## Current operational evidence

The existing FI × sleep LAVA sensitivity run remains active on the external analysis SSD. At 19:30:50 UTC, its runner state reported 20,528 of 29,940 validated receipts, three effective workers, zero launch failures, and zero stale claims. The rolling throughput estimate was 615.83 receipts/hour over 1,806 seconds. This is an interim operational rate, not a completion guarantee. The full captured runner state and source hash are in `lava_runner_state_snapshot_2026-09-27_1930.json`.

## Reviewer handoff

The user confirmed that two distinct human reviewers are ready. The canonical title/abstract queue remains at 56,117 records with zero decisions; the seven-record supplemental queue remains separate. No reviewer messages were sent and no decisions were imported during this update.

## Reproducibility gates

The locked analysis plan, cohort overlap ledger, and review structure validators passed at 19:30 UTC. Results are captured in `core_gate_validation_2026-09-27_1930.log`: the plan remains 12 traits with multiplicity 396, the overlap ledger has 34 rows with exact intersections unknown, and review has 56,117 records with zero decisions and no full texts.

## Remaining gates

The LAVA run has not completed and its 12 trait-level QC gates remain failed per the latest full audit. Screening/adjudication and licensed database exports/full texts remain outstanding. The custom HFRS summary-statistics access route and exact cohort intersections/independent replication remain unresolved. No downstream local-sharing or molecular claims are authorized by this update.
