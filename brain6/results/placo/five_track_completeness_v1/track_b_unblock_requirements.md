# Brain6 Track B: exact unblock requirements

Track B is the legacy **Pair B: insomnia–ADHD** analysis. Its locked sources are the publicly released Jansen 2019 UK Biobank-only insomnia GWAS and Demontis 2023 iPSYCH/deCODE/PGC European ADHD GWAS. Both original aggregates and harmonized GRCh37 files are present and SHA256-verified; their prepared 5,514,402-row pair input is present and verified. The access restriction is on admission of a protected historical result under the frozen Brain6 v3 design, not on obtaining these public GWAS summaries.

The frozen Brain6 v3 lock names this pair as `protected_legacy_pair` and the checked-in `plan-placo` and `execute-placo` functions explicitly forbid rerunning it into v3. The four other selected pairs have complete, receipt-bound QC-passing outputs. Their five-track Bonferroni multiplier is already applied, but the family remains incomplete because Track B has no admissible terminal output.

## Required legacy output package

- `results/track_b/pleiotropy/results/placo/B.full.tsv.gz` and `B.provenance.json`, with full genome-wide rows, per-row P and within-pair BH q, exact allele/coordinate schema, full family counts and QC status. The canonical expected directory is fixed by `config/track_b_pleiotropy_policy.json`.
- Matching `contract.lock.json` and `input_gate.lock.json`; SHA256 references in the sidecar must match the supplied bytes. The current contract lock is `28869d2206c8930c68f253cfa3cc6efb093abc89d666ab6d797abe8e88eb3c95`; the current gate lock is `66916e3416a37ce3c31ffb12d1ccf51a59b1f065ed002efbe7b196f5ffa7c402`.
- The exact aligned input and its checksum, pinned PLACO 0.2.0 source, nuisance parameter checkpoint, runner and materializer source files, complete shard receipts, and successful full-scan process/RAM monitor and terminal exit receipts. The sidecar must bind them and the result output hash. See `result_provenance_required` in the frozen Track B policy.

## Why the mounted failed-partial archive is inadmissible

The mounted archive's `B.full.tsv.gz` is intact (SHA256 `f0eb8b9e8a5b3b4a0ea6151ed75473b0be3d311e463824889a050d83beaaac2e`, 5514399 rows), but its sidecar states `input_gate_lock_sha256=9e9a40aa47711fdfb8191ca8cdfcf84828ae5d3f2715578b9bf91e939260ac97`; the archived gate file hashes to `79d2ae0faa0beb6368496617ab4c88a580a34e7f459f4825d5586e31eecb2101`, and the current gate to `66916e3416a37ce3c31ffb12d1ccf51a59b1f065ed002efbe7b196f5ffa7c402`. The aligned input named in that sidecar (SHA256 `59c21df0bf237e8ede77cfcdbf7b37b54444bc852d39991c88ba74465b2f02e7`) is not supplied with an independently verifiable terminal package. Its runner/materializer hashes differ from the available sources. The archive's V4 amendment explicitly rejects old shards without successful full-scan monitor evidence and requires recomputation under a new implementation fingerprint. The old result cannot be imported into the protected slot.

## Acquisition route and next computation

The source GWAS summaries themselves are public: [Jansen insomnia release](https://cncr.nl/research/summary_statistics/) and [Demontis ADHD release](https://doi.org/10.6084/m9.figshare.22564390.v1). No login, license, or manual GWAS download is needed for the verified local copies. To complete the frozen protected slot, provide a different terminal-valid legacy Track B result package from the authorized run with the exact files and matching receipts above. If none exists, the frozen v3 family cannot be completed as designed. A newly executed public-data result can be analyzed as a separately named sensitivity, but cannot retroactively become the protected legacy output.

Once an admissible protected package exists, verify all hashes and terminal receipts, import without changing frozen outputs, rebuild complete five-track family accounting under the pinned BH and Bonferroni policy, then recompute candidate significance and independent loci before downstream work.
