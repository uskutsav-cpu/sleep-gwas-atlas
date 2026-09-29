# Long-sleep model-alignment pilot: terminal decision (2026-09-27)

**Decision: `NO_FULL_SCREEN_JUSTIFICATION`.** This is a blinded, trait-only diagnostic under the [pre-run addendum](../lava_confirmatory_protocol_v1/longsleep_linear_model_pilot_addendum_20260927.md). It does not admit a new source, authorize a full long-sleep screen, change the canonical LAVA family, or promote any candidate.

## Frozen inputs and execution

- The original Dashti long-sleep archive was source-verified before this run; its SHA-256 is `0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885`.
- The 88 loci were chosen prospectively: 44 canonical `TESTED` and 44 canonical `NOT_RUN`, balanced by chromosome, excluding the 25 PLACO candidate intervals. Selection used no association outcomes. The fixed input retained 209,273 source-verified SNP rows; 87/88 selected loci retained at least 80% of the original shard rows. One originally empty shard stayed in both denominators.
- Frozen configuration SHA-256: `fd3989c7f89b1f1e1a8ee73df3c01bed475e335cc3ced79799a0cb30c10b5079`; materialization receipt SHA-256: `03d7dcb9ae2fed7d86bfdfa41b6cfe19adf78259ab7c8831b2cd23027c5c1462`.
- Two arms ran sequentially with **four independent workers per arm**. All eight worker processes exited 0; each arm's terminal receipt, four partitions, 88 rows, output summaries, logs and SHA-256 digests passed the frozen adjudicator. The source-verified linear-model arm was primary; the same SNPs with LAVA's binary reconstruction were a predeclared sensitivity.
- Original arm outputs, launch configurations and exit receipts are under `/Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1/longsleep_linear_88_v1/results/`. The [paired comparison](</Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1/longsleep_linear_88_v1/results/paired_comparison.tsv>) has SHA-256 `f81ad6d5a80982d1401bbe70cc7f19fa3deb392cf96dff50cb36cb11c7b64edb`. The [terminal decision](</Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1/longsleep_linear_88_v1/results/pilot_decision.json>) has SHA-256 `bc0b88468e693ba081a4cf3e6b19ab291694b68eeecc8dadb76a8a90345418f0`.

| 88 selected loci | Canonical | Primary linear arm | Secondary binary arm |
|---|---:|---:|---:|
| `TESTED` | 44 | 44 | 47 |
| `NOT_RUN` | 44 | 44 | 41 |
| `FAILED` | 0 | 0 | 0 |
| Strict univariate gate pass | — | 0 | 0 |

The primary arm gained **0/88 loci, or 0 percentage points**, against the prospectively required gain of at least 10 percentage points (at least 9 net loci). It had five `NOT_RUN→TESTED` and five `TESTED→NOT_RUN` transitions; 43 of its 44 `NOT_RUN` calls were `LOW_LOCAL_H2_UNDERPOWERED`, and one was `NO_SOURCE_VERIFIED_SNPS`. The binary sensitivity had eight and five respective transitions, for a net gain of three loci; its 41 `NOT_RUN` calls comprised 40 low local h² and one empty source-verified shard. Exactly three loci differed in status between arms; none differed in strict-gate pass. The secondary arm cannot replace the frozen primary by post-outcome selection, and the original source's per-SNP N and model calibration remain unresolved.

## Confirmatory implication

The canonical v3 family remains `FAILED_QC_NOT_PROMOTED`: 3,720/17,465 `NOT_RUN` against a frozen ceiling of 873. Long sleep alone contributes 1,291 `NOT_RUN`; even if every other trait's cells were repaired, 1,291 would exceed the ceiling by 418. The 88-locus pilot is not an estimate of the full-trait rescue rate, but it gives no justification for a larger original-source screen under its prespecified criterion. No full-family or bivariate LAVA run was launched.

The [source-search ledger](../lava_confirmatory_source_triage_v1/source_search_and_admission_20260927.md) documents why the inspected 2023 UKB/MVP ≥10-hour and 2026 SleepChart >8-hour sources are not stronger phenotype-equivalent replacements. A continued primary-source search on 2026-09-27 found no denser, adequately powered, matching ≥9-hour versus 7–8-hour public release. A new source would need a compatible phenotype, enough long-sleep cases, genome-wide statistics with defensible per-variant or model-matched N, source-verified variant and allele semantics, and pair-specific overlap covariance before any new pilot or family run. None is currently admitted.

Frozen canonical outputs, v2, roundoff and exploratory outputs, hashes, and candidate slots were left unchanged.
