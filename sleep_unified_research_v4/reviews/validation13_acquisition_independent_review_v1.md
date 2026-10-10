# Independent original-13 source acquisition review v1

Scope: metadata/code/runtime identities and tiny synthetic controls only. No GWAS or reference body, actual network transfer, actual curl process, actual family lock, estimator or new outcome was used. Production script 81 and its frozen plan were preserved. The review does not admit scientific replication, change the 217-candidate family, select sources by outcomes or adjudicate source licenses.

**Verdict: `MEASURED_TERMINAL_GAPS_DO_NOT_LAUNCH_81_V1`.** Source bindings and the inherited owned-transfer safeguards pass this bounded review, but two family terminal proofs can certify evidence that no longer matches the intended regular-file contract. A separately versioned successor and frozen plan need review before any launch. Waiting for the running raw100 family does not resolve either gap.

Exact production identities:

- Script: `scripts/81_acquire_validation_raw_sources_v1.py`, SHA256 `61893f3a6a1f35daf06bb50b77643d3411c94a684b16d0e5032d02ba1b0e9450`.
- Plan: `manifests/validation_raw_acquisition_plan_v4_1.json`, SHA256 `fd753352967f01b2cfda3747044158540ce29103ee7725071789efb4855965e5`.
- The checker/receipt/seal and all fixture artifacts are bound in the accompanying independent review seal. Prefixes above are relative to `sleep_unified_research_v4/`.

The metadata checker consumes all 53 sealed design dependencies and all 64 unique plan-bound metadata/code/runtime identities. Their current hashes match. The six named design-seal artifact hashes match, and the design's historical source-body records equal the original per-source streaming receipts. The 13 IDs/order also match the original 13-row source-h2 inventory and first occurrence of those sources in the original candidate queue. Exact expected body bytes sum to **10,058,648,185**. This is provenance metadata verification, not a fresh body checksum.

For eleven FinnGen R13 inputs, the original URL includes the exact recorded GCS generation. Frozen HTTP200, size, ETag, Last-Modified, generation, GCS MD5 header and stored length agree with the recorded HEAD probes. Decoded GCS MD5 headers equal original expected MD5s. For the two original MVP outcome sources, GCST90479148 and GCST90479330, current official YAML MD5s match original full-body MD5s. Their current HEAD ETags differ from historical transport metadata; that drift is explicitly qualified. Script 81 requires the current frozen headers and the original full SHA256, MD5 and size before renaming a partial source. HEAD or YAML agreement alone cannot establish current full-body identity.

Fourteen helper function ASTs match reviewed immutable script 52v8 exactly, including owned cleanup, fallback cleanup, quarantine, deferred signals, binding checks, namespace metering and the same persistent network-family flock. The plan retains one network family, no resume/reuse or automatic retry, the 3 GiB internal and 5 GiB SSD runtime floors, 2 GiB owned-transfer RSS stop, 2-hour body deadline, 96-hour family deadline and 300 GiB global meter. Ledger v4_4 reserves the exact original-13 body total. Terminal callbacks recompute all 13 current synthetic body identities and consume both copies of each source receipt in the fixture tests; production callbacks specify the corresponding full-body requirements.

The independent checker executes script 81's function bodies using tiny synthetic files. It changes only six occurrences of the hardcoded aggregate byte budget to the synthetic aggregate, substitutes local fixture paths/resource snapshots/process identities and replaces process launch and the logical lock context with inert stubs. `__main__` is never executed. There are zero real curl or flock calls. All actual body hashing in these controls concerns synthetic text, not GWAS data. The production code and plan are hash-checked again after the controls.

| Control | Measured behavior |
|---|---|
| Normal synthetic family | All 13 sources/receipts accepted; two primary family copies agree; Terminal2 completes and removes PENDING. |
| Mutate both family copies before first hash capture | **GAP:** both copies are changed after persistence to a different status and count12. Their hashes differ from the intended master-payload hashes. Terminal2 nevertheless completes and removes PENDING because those newly observed hashes become the expected hashes. |
| Same-bytes repository logs-parent symlink after seal persistence | **GAP:** moving the logs directory and replacing it with a symlink preserves bytes but changes the parent contract. Terminal2 nevertheless completes and removes PENDING. |
| Mutate only one family copy | Rejected; differing primary hashes; PENDING retained. |
| Same-bytes family-primary leaf symlink | Rejected; PENDING retained. |
| Change current source body after seal persistence | Rejected; PENDING retained. |
| Change source-receipt mirror after seal persistence | Rejected; PENDING retained. |
| Add source or family failure sidecar after seal persistence | Both rejected; PENDING retained. |
| Deferred signal or failed resource floor after seal persistence | Both rejected; PENDING retained. |
| Wrong frozen transport header | First source fails; family stops and PENDING remains. |
| Cleanup failure plus source-receipt persistence failure | Family stops; the quarantine spy is reached twice despite persistence failure; PENDING remains. No actual process group or lock was used. |

There are 13 controls: one normal success, two witnessed vulnerable successes and ten rejection/cleanup controls. Their individual receipts record current/intended primary hashes, PENDING/seal/failure state, synthetic launches and the applied mutation. The synthetic success artifacts are test evidence only and must not be consumed as real source acquisition receipts.

The first gap is at script 81 lines 508–510: primary hashes are captured after `write_new`, without comparison to a digest of the intended master payload frozen before either write. The second is in the terminal identity boundary: `regular()` validates primary paths before Terminal2, but the callback at lines 479–496 rechecks source receipt/body paths rather than both family primary parent chains; the shared helper checks primary leaves. Same-byte ancestor substitution therefore escapes that check.

A successor must precompute the exact intended family-primary payload SHA before persistence, require both copies to equal that fixed digest, and recheck both regular/no-symlink parent chains inside every Terminal2 identity callback. It should also bind regular private PENDING/seal paths and their ancestors. Keep the current all-13 body SHA/MD5/size, both source-receipt copies, deferred termination, resource/deadline, failure-sidecar, cleanup and quarantine requirements. Freeze the successor's code/plan/review identities; do not edit or retrofit an apparent PASS onto script 81 or this failed review.

The initial checker metadata allowlist omitted a `.md` design-review dependency and rejected it before any fixture was created. A separate repair receipt preserves that draft checker hash and the metadata-only failure. The corrected checker then completed all dependency checks and controls. This reviewer performed no production preparation, transfer, lock acquisition, body reuse, investigator contact or scientific analysis.
