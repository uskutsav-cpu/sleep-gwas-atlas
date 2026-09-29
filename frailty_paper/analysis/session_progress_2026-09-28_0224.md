# Frailty LAVA recovery checkpoint — 2026-09-28 02:24 UTC

Coordinator state at 2026-09-28T02:24:21.900210+00:00: 21,865/29,940, two active claims (PIDs 15810, 15811), zero launch failures and stale recoveries. Full audit 2026-09-28T02:19:17+00:00–2026-09-28T02:20:10+00:00 validated 21,844; inventory reached 21,846; receipt/claim/duplicate issues zero and both PIDs live. All 12 locked gates fail (3,804 negative variance, 252 no-reference-SNP, 2 other failures).

There are 8,075 slots remaining. Receipt-event logs show 31 completions in the latest 30 wall-clock minutes (62/hour) and 30 in the latest 15 minutes (120/hour), with no duplicate events. The coordinator retains 223.64/hour from its last completed 30.0-minute window. Receipt-event gaps in the last six hours reach 44.0 minutes, so no current ETA is defensible. Details and source hash are in `lava_rolling_throughput_2026-09-28_0224.json`.

Figure 4 was refreshed from the full audit. Host memory free was 47%; external SSD free is 1.5 TiB. Current CPU/swap-total metrics remain unavailable in the restricted shell. Frozen lock/input hashes are unchanged. Reviewer archives remain ready; delivery routing is still pending.
