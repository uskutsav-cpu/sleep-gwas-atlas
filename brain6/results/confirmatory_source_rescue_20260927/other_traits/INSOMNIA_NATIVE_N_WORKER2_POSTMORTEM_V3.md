# Insomnia native-N worker-2 postmortem v3

**Purpose:** identify the exact locus and error hidden by the v2 runner's pre-write zero-failure assertion. This postmortem is a separate diagnostic execution after the v2 terminal failure. It cannot satisfy v2's 88-locus advancement rule, replace a v2 worker result, authorize a full trait screen, amend the family gate, or promote a source/candidate.

Replay **only v2 worker 2's 22 loci** against the unchanged, SHA-256-verified v2 materialized input shards, selected.loci, per-chromosome UKB reference, binary source-cohort cases/controls (109,402/277,131), LAVA 0.1.5, random seed 20260928, canonical locus/univariate settings, exact roundoff patch, and strict p threshold. Workers 1, 3, and 4 are not rerun. The R diagnostic records all 22 status rows and a summary in a new exclusive v3 SSD output root before returning nonzero if any locus failed. It makes no conditional code, N, model, QC, threshold, or source selection changes based on the v2 outcome.

The postmortem is allowed to describe the error and source/method implication only. If the error is transient or cannot be reproduced, v2 still fails its frozen rule. Diagnostic outputs are marked `POSTMORTEM_ONLY_NO_ADVANCEMENT` and are never merged into the v2 aggregate.
