# Five-track candidate-region reconciliation

The checksum-bound five-track PLACO output contains 25 pair-specific candidate
intervals. Grouping intervals that overlap on the same chromosome produces **20
geographic candidate regions**. Four groups contain more than one pair-specific
candidate; all 25 original candidates remain individually recorded in
`pair_candidate_members.tsv`.

| Geographic group | Pair-specific members | Existing frozen UKB lead LD |
| --- | ---: | --- |
| chr5:103447968-104481726 | insomnia–ADHD; insomnia–MDD | rs2431108–rs40465, r² = 0.98473097 |
| chr7:113559156-114571035 | insomnia–ADHD; insomnia–MDD | rs1476535–rs2894699, r² = 0.96252403 |
| chr15:84700520-85700520 | long sleep–bipolar; long sleep–SCZ | same lead rs1051168 |
| chr17:43049526-45343136 | two long sleep–Parkinson candidates; long sleep–SCZ | overlapping edges have r² = 0.76631925 and 0.88608216 for listed lead pairs |

The LD values are read from the already frozen five-track candidate-variant
table. No reference data were recomputed or substituted. This is a **geographic
reconciliation**, not an independent-signal count: the original clumping rule
defines leads within each pair, and local LAVA v3 remains failed QC. The 20
groups cannot be promoted for downstream inference on this evidence alone.
`geographic_region_groups.tsv`, `cross_pair_overlap_edges.tsv`, and
`provenance.json` give the complete membership, overlap, and source hashes.
