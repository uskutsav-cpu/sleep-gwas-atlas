# Independent screening archive handoff validation — 2026-09-26 21:02 UTC

The two reviewer-specific ZIP archives were extracted into a temporary staging copy of the original packet tree. The packet importer ran in dry-run mode for each slot. This did not apply decisions or write to the canonical review queue.

| Check | Result |
|---|---|
| Queue SHA-256 | `3ae3fb3b4765c85536c3ccf64d20297a632af892e2c66e622ee665db1ea43b2e` (matches packet manifest) |
| Protocol SHA-256 | `1013526fa40de6708ed9e547e09c1cb0ae3de0c27d8cbe2d5d04a99ffa796abb` (matches packet manifest) |
| Reviewer 1 archive | 36,268,755 bytes; SHA-256 `88bfc194c6c900e87f335d67a2dbda2c3e04daf65f6110e7a009ec9c1ddf9d41`; 113 TSV batches; ZIP integrity passed |
| Reviewer 2 archive | 36,268,749 bytes; SHA-256 `0b87aaebfdea35022ccf59aefb648d06f40e5ece665b2e9921fa35a19b42f4f6`; 113 TSV batches; ZIP integrity passed |
| Reviewer 1 importer dry run | Exit 0; 56,117 queue rows; 113 packet files; zero decisions imported; `applied=false` |
| Reviewer 2 importer dry run | Exit 0; 56,117 queue rows; 113 packet files; zero decisions imported; `applied=false` |
| Canonical review queue | Remains at the manifest SHA-256 above; no reviewer decisions or full-text decisions |

Archives are local and have not been sent to anyone. Return handling and import commands are documented in `review/screening/REVIEWER_PACKET_HANDOFF.md`.
