# Track B recovery summary

**Decision: BLOCKED_V3_ADMISSION. No protected-slot promotion was performed.**

The archived B ledger is intact and its row counts, schema, P/q ranges, family counts, provenance hash, exact source-code hashes, V2 terminal gate, B-specific successful full-scan monitor, and cleanup attestation reconcile. The 276 original shard files and aligned input were intentionally cleaned after publication; the original V4 clean-package validator explicitly allows this state when its cleanup receipts and retained canonical output pass. Exact regeneration of the aligned gzip from frozen dense inputs remains unverified.

The prior archive audit overstates two blockers. In the V2 bridge, `scripts/140_materialize_track_b_placo_pair_v2.py:375` intentionally sets the provenance key `input_gate_lock_sha256` to the **V2 terminal gate**, whose archived SHA is `9e9a40aa...`. The B-specific V4 monitor succeeded; the V4 amendment's failed-monitor discussion concerns superseded **Pair A** shards. Exact B materializer and runner sources were recovered from the backup tar.

The live frozen Brain6 v3 audit still explicitly marks B inadmissible and expects the current V1 input-gate identity (`66916e34...`). The archived V1 gate differs (`79d2ae0f...`) because its readiness-gate hash differs, although the other V1 gate fields match. The archived sidecar uses V2 terminal-gate semantics. The legacy 17-field ledger also differs from the 10-field v3 pair schema and needs an explicit checked adapter. No frozen Brain6 v3 importer was found that authorizes this alias/cleanup attestation as a protected legacy package. The existing five-track auditor writes into an already populated output directory and hard-codes the archived package as rejected; executing it would neither validate nor safely promote this package. The original V4 `--verify` path also performs lock/repair/rebuild operations and requires original filesystem identities, so it cannot be run unchanged against this copied archive in a read-only audit.

## Exact remaining requirements

1. Freeze a narrow, prospective admission amendment for protected *legacy* B that explicitly recognizes the V2 terminal-gate alias, separately verifies the V1 gate/source identity, accepts the original V4 cleanup attestation after independent hash, monitor, code, and output checks, and defines a deterministic legacy-to-v3 row-schema adapter. Do not change statistical thresholds or the five-track family.
2. Provide a read-only admission validator or restore the original V4 checkout with stable identities in a separate authorized environment; validate the complete clean package, including all originally retained A/B/CONTROL family receipts.
3. If that validator rejects the historical package, obtain a different terminal-valid authorized legacy package or approve a versioned replacement-run protocol; the completed prospective public-input sensitivity must remain sensitivity evidence.

The original aligned gzip may be reproducible from the frozen dense files and recovered code, but byte identity has not been demonstrated. Original shard bytes cannot be reconstructed without recomputing statistics. Under the original V4 clean-package rule, neither was required to remain after successful attested cleanup.
