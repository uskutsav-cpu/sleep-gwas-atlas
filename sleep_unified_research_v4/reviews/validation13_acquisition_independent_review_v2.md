# Independent original-13 acquisition successor review v2

**Disposition: the two v1 family-primary corrections pass; script 82 v2 is not admitted for execution.** Two measured issues remain: seven prior-review artifacts are omitted from the frozen dependency graph, and private terminal-seal content can change before its first hash capture without invalidating success. Preserve this failed qualification and use a separately frozen successor.

Frozen identities:

- `scripts/82_acquire_validation_raw_sources_v2.py`: SHA256 `5b8fc6ac796ebc9dabd783d52fa0069ea498fb9904164b53dffa494cf263e783`.
- `manifests/validation_raw_acquisition_plan_v4_2.json`: SHA256 `fe63b4c911ea85c748141477f6f6769c2588f24590f205b9b738bd45693f26f7`.
- Prior v1 review seal: SHA256 `d40c0524acfe924c3ae86cbdfec1fe8d7df30afee3c5aa32cb377e6361223de9`.

Paths above are relative to `sleep_unified_research_v4/`. Production 81/82 and both plans remain unchanged. This is a narrow metadata/code and tiny-fixture review; no actual curl process, network transfer, family flock, GWAS/reference body or scientific estimate was used.

All 67 plan-bound metadata/code/runtime identities match. The original 53 metadata dependencies and six design-seal artifact identities remain consumed. All seven artifacts in the v1 review seal also independently hash-match in this review. Original source order, identifiers, URLs, SHA256/MD5/size and transport expectations are unchanged; only body/mirror/terminal output routes move to a new v2 epoch. The exact source family still contains eleven generation-pinned FinnGen R13 objects and the two original MVP outcome objects, totaling 10,058,648,185 bytes. MVP current-versus-historical ETag drift remains explicitly qualified; no source body is substituted. Fourteen unchanged transfer/ownership helpers retain the v1-reviewed AST identities. Shared-family locking, cleanup/quarantine, deferred signals and resource/deadline limits are unchanged.

The v1 seal's seven review bindings are under `review_artifact_sha256`. Script 82 instead reads `prior_seal.get('file_sha256', prior_seal.get('artifacts', {}))`. For the actual frozen seal that expression returns an empty mapping. The plan therefore binds the v1 seal itself, but none of its seven review artifacts. This is a schema-consumption failure, not a hash mismatch or permission issue. The checker records all seven missing paths and verifies their present hashes separately. A successor must explicitly consume the actual mapping and bind its full artifact set, together with this v2 seal and its artifacts.

The primary-family correction precomputes the intended master JSON digest before either receipt write, compares both copies against that fixed intent, and carries the fixed hashes into every terminal identity callback. Both family primary leaves and all ancestors are checked for symlinks; private PENDING and present seal paths receive regular-file checks. The two original false-PASS witnesses now reject correctly.

Twenty-one tiny controls use the same v1 fixture method: production function bodies, six aggregate-budget constants scaled to a synthetic total, local fixture routes/resource snapshots and inert process/lock stubs. The production main entry point and real process/lock APIs are not invoked. Results are **one healthy success, nineteen rejections/cleanup controls, and one new vulnerable success**.

| Boundary tested | Measured result |
|---|---|
| Healthy 13-source synthetic acquisition | Success; two intended family-primary digests match; PENDING removed. |
| Both primary copies changed before first hash capture | Rejected against precomputed master intent; PENDING retained. This closes v1 witness one. |
| Same-bytes repository logs-parent symlink after seal persistence | Rejected by terminal callback ancestor checks; PENDING retained. This closes v1 witness two. |
| One primary changed before capture; primary missing before capture | Both rejected; PENDING retained. |
| Primary leaf symlink, primary corruption or missing primary after seal persistence | All rejected; PENDING retained. |
| Changed/missing current body or changed source-receipt mirror after seal persistence | All rejected; PENDING retained. |
| Source/family failure sidecar after seal persistence | Both rejected; PENDING retained. |
| Deferred signal, failed internal floor or wrong frozen transport header | All rejected; PENDING retained. |
| PENDING leaf symlink/corruption or seal leaf symlink after persistence | All rejected; PENDING retained. |
| Cleanup failure plus receipt persistence failure | Failure remains preserved; quarantine spy is reached twice; PENDING retained. |
| Seal content changed immediately after write, before helper's first seal hash | **GAP:** changing `pending_sha256` to 64 zeros is accepted. Terminal2 removes PENDING, and its read-only `require_committed` consumer also accepts when supplied the fixed intended family-primary digests. |

The last witness distinguishes file-path regularity from the seal's intended content. Shared Terminal2 captures `seal_sha` after `save_new`. The callback validates the seal path but does not compare its bytes or schema to a seal payload fixed before persistence. A changed pending digest therefore becomes part of the newly observed seal hash. The helper's later same-hash comparison does not restore the original intention. The separate `validation13_private_seal_consumer_probe_v2.json` records acceptance by both the producer and read-only consumer. This is solely synthetic corruption evidence; it is not a real acquired source or scientific failure.

A successor should precompute the exact intended Terminal2 seal payload, including the constructor-frozen binding, original `terminal.pending_sha`, fixed intended family-primary digest mapping and the existing terminal fields. At each callback, any present regular seal must match those expected bytes; never derive expected seal content from the persisted file. Preserve the current regular leaf/ancestor checks, full current13 body hashes, both individual receipt copies, signal/resource/deadline/failure vetoes, quarantine and no-retry policy. Fix the prior-review mapping mechanically and freeze the corrected code, plan and complete transitive review graph before another narrow delta review. No unchanged transfer campaign needs repeating for this review.

The checker and fixture receipts document zero real transfers or locks. All production code/plan and current metadata/runtime hashes are checked again before sealing this report. Neither a v1 nor v2 production acquisition is launched or cleared here; completion of raw100 remains a separate prerequisite after an eligible successor is admitted. No new replication, phenotype filter, testing family, source right or gzip/native-pipeline claim is introduced.
