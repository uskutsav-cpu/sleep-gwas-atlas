# Independent original-13 acquisition correction review, v3

Verdict: **QUALIFIED_PASS_90_V3_PRELAUNCH_TERMINAL_DELTA_ONLY**. The measured 81 family-primary intent/parent gaps and the measured 82 review-map/private-seal intent gaps are closed for the tested boundaries. This review grants no source-body or scientific admission. Execution still requires the parent operational admission and the existing transfer-family lock after the running raw100 acquisition ends.

The frozen executor is `90_acquire_validation_raw_sources_v3.py`, SHA256 `8d0765930a346d1097d15a7b5fc3bea294f782b4eb8f52a98f894bbf693152cc`. The frozen plan is `validation_raw_acquisition_plan_v4_3.json`, SHA256 `8e93d548163bda284b6292b9780962060bbc3ff4afde44e5a916c80f3b69b708`. Both were rehashed after the controls. Production 81, 82, both earlier plans, both rejection seals and their artifacts remain unchanged. No candidate, scientific source, shared Terminal2 helper, historical output, threshold or filter was edited.

## Prior evidence actually consumed

90's actual `prepare()` reader uses the nonempty `review_artifact_sha256` maps from v1 (seven artifacts) and v2 (six artifacts), verifies their digests and adds them to the frozen dependency map. Both rejection-seal digests and the v2 production executor/plan are also bound. The current plan has 83 dependencies. All 83 metadata/code/runtime dependencies were checked before and after the terminal controls, including every one of those 13 prior review artifacts.

The preparation control executes the original candidate reader with real metadata hashes while suppressing every mkdir, write, physical-mount query and runtime subprocess. Its healthy result reconstructs exactly the current plan's dependency map and members and captures two identical intended plan copies. Four negative controls inject an empty v1 map, a wrong v1 artifact digest, an empty v2 map, and a wrong v2 artifact digest into the metadata reader. All fail before any intended plan write. The original metadata bytes are never changed.

## Fixed terminal intent

Before family receipt persistence, the caller computes the intended complete family JSON digest and freezes both primary paths to it. Before `terminal.commit()`, the caller computes the complete Terminal2 seal serialization using the fixed original binding, fixed primary digests, exact pending path, constructor-frozen `terminal.pending_sha`, status and success qualification. Its serialization matches Terminal2's `indent=2, allow_nan=False` JSON plus newline. Every identity callback, when the seal exists or is a symlink, demands regular leaves and all regular parents and that precomputed digest. Thus an observed persisted digest cannot become a new acceptance baseline.

The healthy terminal control succeeds, removes PENDING and is accepted by the read-only consumer with the fixed primary map. Its persisted seal digest equals the independently captured intended pre-write digest. The prior witness changing `pending_sha256` to 64 zeroes after seal persistence but before first seal hash capture now fails with PENDING retained. Changes to the seal binding, primary map, success qualification and whitespace at the same boundary also fail. The two prior 81 witnesses—both family receipt copies changed before first receipt hash capture and same-byte log-parent symlink after seal persistence—also fail.

There are 25 terminal controls: one healthy success and 24 rejections/quarantine controls. Negative coverage also includes missing/corrupt family receipts, source body/receipt drift, source/family failure addenda, pending/seal leaf symlinks, pending byte drift, deferred signal, resource failure, wrong transport headers and cleanup plus receipt-persistence failure. Every negative retains PENDING and lacks a successful read-only consumption result. Cleanup-persistence failure exercises quarantine spies rather than live processes. Results and exact error messages are in `validation13_acquisition_independent_receipt_v3.json`.

## Unchanged transfer and source scope

All 13 historical source identifiers, order, URLs, original SHA256/MD5/size and frozen transport fields are unchanged from 82. Only the private body route changes to `validation_raw_replay_v3`; receipt/admission/plan/terminal output routes are versioned. The total remains 10,058,648,185 bytes: 11 exact FinnGen R13 generation URLs and the original two MVP accessions. Both MVP current ETags remain explicitly qualified as different from historical transport metadata; the original full SHA256/MD5/size and official YAML MD5 remain required. Current HEAD evidence does not replace the original body identity gate.

All 14 inherited transfer/ownership helper ASTs equal 52 v8. The acquisition function equals 82 after normalizing its new receipt-output directory, and transport-header checks are identical. The single transfer-family mutex, no retry/reuse/resume behavior, 3 GiB internal and 5 GiB SSD runtime floors, 2 GiB owned transfer RSS limit, two-hour per-body limit, 96-hour family limit, 300 GiB namespace meter and v4_4 resource ledger are unchanged. Current full SHA256/MD5/size of all 13 bodies and both individual receipt copies remain required by terminal identity callbacks before and after persistence.

## Bounds and qualifications

This is a narrow code/metadata and synthetic-fixture review. It reads only the explicit frozen metadata/code/runtime dependency list, with a 20 MiB per-file metadata cap; the checker hashes with 64 KiB buffers. The terminal fixtures retain 1,812 file/symlink artifacts and 3,244,200 regular bytes. Historical `.gz` filenames in these fixtures contain tiny plain ASCII synthetic strings. AST adaptation changes only six aggregate body-byte constants to the tiny synthetic total; actual executor function bodies otherwise run unchanged. Curl/runtime and ownership are inert mocks, and the physical family flock is never opened or acquired. No source or reference body, network transfer, genuine worker, fit or scientific estimate is used.

This PASS covers the witnessed correction boundaries under the existing owned/private immutable-metadata contract. It does not certify crash recovery, hostile arbitrary changes after the last identity callback or after PENDING removal, real curl termination under load, actual SSD capacity, source-content correctness, EOF/replay, filtering, statistical validity or independent replication. The shared read-only consumer still relies on the owned private seal's immutability; this caller-side correction does not change that helper or retroactively strengthen previous completed stages. Public source rights and individual cohort-overlap questions are not inferred from this transport review. Later work needs its own concrete source hashes and stage admission.

## Review artifacts

- `validation13_acquisition_independent_checker_v3.py`: exact metadata, preparation and tiny terminal controls.
- `validation13_acquisition_independent_receipt_v3.json`: 83 dependency identities, five preparation-reader controls and 25 terminal results.
- `validation13_original_source_transport_delta_receipt_v3.json`: unchanged source/transport/guard comparison.
- `validation13_acquisition_fixture_manifest_v3.json`: exact retained synthetic artifact hashes and symlink targets.
- `validation13_acquisition_independent_review_seal_v3.json`: immutable review artifact map and production binding.
