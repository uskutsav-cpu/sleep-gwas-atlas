# Long-sleep 88-locus pilot: source-verified input preparation

The [pre-run protocol](../lava_confirmatory_protocol_v1/longsleep_linear_model_pilot_addendum_20260927.md) and [88-locus selection](longsleep_linear_88_loci.tsv) were committed before these inputs were created. The selection includes 44 frozen canonical `TESTED` and 44 `NOT_RUN` loci, two of each per autosome. All 25 existing pair-specific candidate intervals were excluded by coordinate overlap, removing 38 LAVA blocks from the selection pool. The immutable canonical family was not edited.

The [materialization receipt](</Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1/longsleep_linear_88_v1/materialization.receipt.json>) has SHA-256 `03d7dcb9ae2fed7d86bfdfa41b6cfe19adf78259ab7c8831b2cd23027c5c1462`. Preparation rehashed every selected original canonical shard, scanned all **14,661,601** rows of the original Dashti archive, and matched **219,357** unique IDs from **219,649** selected shard rows. It quarantined **292** source-duplicate SNP IDs and excluded **10,084** more rows under the frozen INFO `>=0.95` and MAF `>=0.01` filter. No selected row failed source coordinate, allele, finite-value or signed `Z=BETA/SE` checks. Both arms received byte-identical, newly versioned SNP/Z/`N=339926` shards with separate model declarations and file hashes.

| Canonical stratum | Loci | Original shard rows | Retained rows |
|---|---:|---:|---:|
| TESTED | 44 | 109,724 | 104,752 |
| NOT_RUN | 44 | 109,925 | 104,521 |
| Total | 88 | 219,649 | 209,273 |

One preselected NOT_RUN block, locus 967, had zero SNPs in its original long-sleep shard. It remains in the fixed denominator as `NO_SOURCE_VERIFIED_SNPS` for both arms. Across all 88 loci, **87/88** retained at least 80% of the original shard rows; the sole exception is the already empty block. This source-only coverage result does not establish LAVA local-h² testability or justify full-family admission.
