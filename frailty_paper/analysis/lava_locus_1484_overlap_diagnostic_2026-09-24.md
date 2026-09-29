# LAVA locus 1484 input-overlap diagnosis

**Generated:** 2026-09-24 11:38 UTC
**Scope:** read-only inspection of the active `accel_sleep_duration` × FI first pair; no receipts, inputs, references, or thresholds were changed.

The 11:32 LAVA log for locus 1484 (GRCh37 chr9:140,097,760–141,146,682) reports that `process.input` found 250,205 variants shared across the pair inputs genome-wide, then `process.locus` returned NULL with `none of specified SNP IDs are present in reference data`. The global shared count does not imply any shared variants in each locus.

A direct ID-set intersection against the locked UK Biobank LAVA v1.1 chr9 `.info` SNP IDs found:

| Set within locus 1484 | SNP IDs |
|---|---:|
| Reference | 1,501 |
| `accel_sleep_duration` input shard | 1,483 |
| FI input shard | 0 |
| Joint FI × sleep × reference | 0 |

As a neighboring-locus check, locus 1483 (chr9:138,995,792–140,097,759) had 1,870 FI IDs, 2,433 sleep IDs and 1,870 joint IDs against the reference. This supports a locus-specific empty FI/reference intersection at 1484 rather than a general chromosome-9 reference read failure. LAVA's pinned `process.locus` implementation intersects locus reference SNP IDs with `input$analysis.snps` before calling `read.ld`; with no pairwise IDs in this interval, the returned empty-locus error is consistent with the observed inputs.

The frozen interim collator continues to count locus 1484 as a process failure. No exception or threshold change is made. The global 1% family-completeness gate had already failed from other flagged loci; local genetic-sharing inference remains unavailable.

## Input identity

| Input | SHA-256 |
|---|---|
| Frozen locus file `ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile` | `462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882` |
| LAVA reference `lava-ukb-v1.1_chr9.info` | `8b2c4b0762e573f46b13f9cc198856490a2a356fc4be39319ab312db40a8a6e0` |
| FI `chr09.sumstats.tsv.gz` | `f45aef246e04c86cbab9ba49b5ccee459e0a6d34e12a1b03d5e89de6d9e5909d` |
| Sleep `chr09.sumstats.tsv.gz` | `572f7078554763fbda5f5b6d35bb1426d2106dfb61e4c1cdcb0d8e74ed74f226` |
| Pair `input_info_chr09.tsv` | `646db62fdda880482f416a08596a9cf36d08c78eb4ac77840a6d9f3ddaaa26bd` |
| Active run lock | `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` |

## Chronotype pair replication of the locus-specific empty-input failure

**Checked:** 2026-09-24 18:18 UTC. The active `chronotype` pair reached the
same chr9 interval and returned the same `process.locus returned NULL` error
at locus 1484. Its 18:15:36 UTC log reports 250,218 variants shared across
the pair inputs, followed by the empty reference-ID error. The larger
genome-wide shared count does not imply an intersection in this particular
locus.

A direct interval scan of the active pair's pinned chr9 reference and summary
input files reproduced the same locus-level condition:

| Set within locus 1484 | SNP IDs |
|---|---:|
| Reference | 1,501 |
| Frailty input | 0 |
| Chronotype input | 1,483 |
| Joint frailty × chronotype × reference | 0 |

The reference `.info` SHA-256 remains
`8b2c4b0762e573f46b13f9cc198856490a2a356fc4be39319ab312db40a8a6e0`. The
frailty input checksum is
`f45aef246e04c86cbab9ba49b5ccee459e0a6d34e12a1b03d5e89de6d9e5909d`; the
chronotype input checksum is
`af954957636d9e9e7b400f9d138096ddf2667a5ad8e8f565874b7608b4a64d77`. This
supports the same missing FI/reference overlap mechanism previously observed
for the completed first pair. The status remains `PROCESS_FAILED` and counted
under the frozen v1 rule. In the latest chronotype prefix, this is a 21st
unexpected process failure in addition to the 20 chr6 MHC-exclusion failures;
no threshold or inference status changes.
