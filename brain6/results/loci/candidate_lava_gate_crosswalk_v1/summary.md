# Candidate-region LAVA gate crosswalk (diagnostic)

The 19 post-Track-B blocked region rows contain 21 lead variants. Every lead
maps uniquely by its GRCh37 position to a frozen LAVA block. Only three region
rows have both canonical trait univariate P values strictly below
`0.05 / 17,465`; they represent two unique pair–block combinations, both at
block 2207 on chromosome 17. Two long-sleep–Parkinson candidate rows share
one pair–block combination. The third row is long sleep–schizophrenia.

| Pair–block | Candidate region rows | Existing diagnostic local-rg P | Existing diagnostic full-family q |
| --- | ---: | ---: | ---: |
| longsleep__parkinson, LAVA 2207 | 2 | 2.21035e-05 | 0.2757411625 |
| longsleep__scz, LAVA 2207 | 1 | 0.120592 | 1 |

These pairwise values are from the **pre-existing diagnostic probe**; this
crosswalk did not run bivariate LAVA. Neither probe passes the recorded
full-family correction. The canonical LAVA v3 family remains
`FAILED_QC_NOT_PROMOTED`, and full-family PLACO locus deduplication and exact
independent replication remain unresolved. No region is promoted or rejected
by this crosswalk. See `crosswalk.tsv` and `provenance.json` for every region,
lead coordinate, canonical univariate cell, strict gate, and input hash.
