# Brain6 confirmatory terminal status — source-rescue audit

**Canonical LAVA v3:** `FAILED_QC_NOT_PROMOTED`. The frozen seven-trait family has **13,745 TESTED**, **3,720 NOT_RUN**, **0 FAILED** across 17,465 planned trait–locus cells; the locked maximum is **873 NOT_RUN**. No protected region was promoted. The v3, v2, roundoff, pilot receipts, hashes, and output directories were preserved.

| Trait | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|
| Long sleep | 1,204 | 1,291 | 0 |
| Insomnia | 1,824 | 671 | 0 |
| Parkinson disease | 1,859 | 636 | 0 |
| MDD | 1,906 | 589 | 0 |
| ADHD | 2,272 | 223 | 0 |
| Bipolar disorder | 2,327 | 168 | 0 |
| Schizophrenia | 2,353 | 142 | 0 |

Source: [frozen canonical trait receipt](../lava/canonical_v3_not_run_by_trait_v1.tsv). The exact Dashti public release cannot yet be represented with verified per-variant N and model-matched LAVA semantics; its official EBI mirror is byte-identical. The prespecified 88-locus total-N linear pilot remained **44/88 TESTED**, equal to canonical, and did not advance. The binary sensitivity's **47/88** is not a source-admitted alternative. A separate MDD2025 pilot reached **108/111** numerically but remains on scientific sample-size interpretation hold. No full-family rescue was launched.

Passing requires at least **2,847** fewer NOT_RUN cells overall. Long sleep alone exceeds the frozen ceiling by **418**, even if all other traits were repaired; perfectly repairing long sleep alone would still leave **2,429** NOT_RUN. The local [rescue arithmetic](../lava_confirmatory_rescue_plan_v1/planning_summary.json) identifies a theoretical four-trait perfect-repair combination leaving 533, but none is an admitted input set. Other binary source headers were [triaged](OTHER_BINARY_SOURCE_TRIAGE.md); further model and N validation would precede any pilot.

**External action to reopen exact long sleep:** custodian confirmation or a variant-keyed N/NMISS file plus exact released BOLT-model interpretation, as specified in the [unsent request](LONG_SLEEP_DATA_REQUEST.md). Current [blocker](LONG_SLEEP_EXTERNAL_BLOCKER.md) gives the next validation and frozen-run sequence. Exploratory Yale results remain separate and cannot change this status.
