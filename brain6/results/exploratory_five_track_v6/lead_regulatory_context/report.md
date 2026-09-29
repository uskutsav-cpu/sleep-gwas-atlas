# Five-track lead regulatory overlap — exploratory

All 27 pair-specific lead rows (26 unique variants) were queried at their exact
GRCh37 coordinates using the official Ensembl GRCh37 regulatory overlap endpoint.
Six protected Track B responses were reused byte-for-byte; the other 20 unique
queries are individually cached and hashed. A returned feature means coordinate
overlap only, not regulatory function, variant causality, or shared-trait support.
An empty response cannot exclude regulation by other variants in the locus.

Exact-lead regulatory feature rows: 6.
Lead rows with at least one overlap: 5 / 27.
All LAVA local-rg and final tiers remain BLOCKED_LAVA.
