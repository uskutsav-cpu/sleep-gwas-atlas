# Insomnia native-N pilot v4: terminal result

The v4 protocol was frozen in commit `4302d55a` before the corrected worker ran. It reused the exact v2 outcomes for workers 1, 3, and 4 and reran only worker 2, after checking the original source, README, 88 selected loci, all materialized shards, reference, software, rules, and successful-worker receipts. V2 and v3 remain immutable. The new output root is `/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v4`.

The zero-SNP locus 950 is now **`NOT_RUN_NO_SNPS`**, with `n_snps=0`; it is never called TESTED. Worker 2 wrote 22 rows, including 10 TESTED, 12 NOT_RUN, and **0 FAILED**, and exited 0. Tests cover zero-row, malformed-header, missing-file, inconsistent row-count, valid-shard, and genuine LAVA failure paths. The integration failure test confirmed that a nonempty one-SNP shard yields a structured FAILED row and a nonzero process exit.

The first v4 controller invocation wrote both combined tables but failed at its final receipt lookup because it searched for the v2 materialization receipt under the v4 root. No worker was rerun or output overwritten. `aggregation_repair_v1` links to that original receipt and the four verified worker results, calls the **unchanged v2 aggregation and decision routine**, and verifies that its combined-table hashes are byte-identical to the first v4 tables. The separate repair receipt records this path-only correction.

| Metric | Frozen canonical on 88 loci | V4 native-N pilot | Change |
| --- | ---: | ---: | ---: |
| TESTED | 44 | 47 | +3 |
| NOT_RUN | 44 | 41 | −3 |
| FAILED | 0 | 0 | 0 |
| `LOW_LOCAL_H2_UNDERPOWERED` NOT_RUN | 43 | 40 | −3 (6.98% relative) |

All three conversions are canonical `NOT_RUN` to v4 `TESTED` at LAVA loci 2099, 2374, and 2462. There are no losses of canonical TESTED loci. The TESTED fraction gain is **3/88 = 3.41 percentage points**. It fails both original advancement alternatives (≥10 percentage points TESTED gain or ≥25% relative low-h2 NOT_RUN reduction), with zero failures. The terminal rule result is **`NO_MATERIAL_PILOT_IMPROVEMENT`**. No full 2,495-locus insomnia screen is authorized by this pilot; the Jansen input remains `ADMISSIBLE_SENSITIVITY_ONLY`, and no candidate or primary LAVA cell is promoted.

The pilot sampled two hash-ranked loci per chromosome and prior canonical status, excluding all 38 candidate-overlapping LAVA blocks. A design-weighted **diagnostic** extrapolation of the three recovered cells over the 670 eligible canonical insomnia NOT_RUN blocks is 17.5 cells (successes occurred on chromosomes 15, 20, and 22). This is not a validated or unbiased estimate for candidate blocks, and it is not credited against the primary family gate. **Validated recovery remains 0; the canonical family projection remains 3,720 NOT_RUN against the frozen ceiling 873.** The exact calculation and hashes are in `insomnia_native_n_pilot_v4_audit/summary.json`.

Key SHA-256 receipts: v4 config `9cd1f0a18a17e370e68e1e6170f8c9edda1a7b7ad86e5887505bd46e0d14d47c`; v4 worker-2 TSV `c51fce87d9f3f04d5f710be62b70ae8629f3fc8a21ad3fd56c6b8f71d610df89`; complete aggregate `129744175ffd97c2bff399befa39d3b815da4e4105d2c83d15492e3381fc5ff9`; comparison `d937e1ff6280383ee567628e43363aa2cc11982cc027a7777f2dc19f0ed538be`. The frozen canonical aggregate remains `ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352`.
