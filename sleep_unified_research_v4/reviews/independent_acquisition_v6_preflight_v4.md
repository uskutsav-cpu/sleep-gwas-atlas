# Independent acquisition v6 preparation review

Verdict: **PASS for the exact prepared v6 executor and plan, subject to the complete terminal consumer contract below. No material blocker found.** This is a preparation and operational identity review; it creates no execution admission and makes no acquired-family, gzip/EOF, pipeline, replication or biological claim. Root admission was absent at the audit. No transfers, estimators, source decompression or additional real GWAS body reads occurred.

Reviewed executor `scripts/52_acquire_extension_raw_sources_v6.py`, SHA256 `d408b245c362c274e95f8c33727ba016367ed14b1e0cd0a77adfd29ae577bf8d`, and plan `manifests/extension_raw_acquisition_plan_v4_6.json`, SHA256 `393af7ffc3231ceddf829ae669558001a7411bb4a7917b950a713c41546237fe`. All 137 bound dependencies and the SSD plan copy match; 141 metadata files, totaling 4,413,500 bytes per verification pass, were independently hashed in 64 KiB buffers before and after the controls and were unchanged. Exact earlier failures, code versions, receipts, prefixes and review seals remain bound and preserved.

## Identity and resource design

All 100 ordered original members retain their URL, immutable S3 version, filename, original byte count, MD5 and full SHA256. Members 1–2 retain their original v4 body and receipt paths. The unchanged `extension_replay_common_v2.py` receipt/header/stat gate independently accepts their bound metadata; this review reuses the prior two-pass whole-body checkpoint proof rather than re-reading those bodies. Members 3–100 use 98 distinct new v6 body paths. No source substitutions, options, phenotype membership, scientific thresholds or retry policy changed.

Source 3 retains the original 1,276,116,992-byte prefix, SHA256 `27e643fe9fe7d5a98fbfa21ceecb5961605975a2743a4d4df6e1f82c2fa6d681`. The executor copies rather than edits that prefix, requires its exact byte/hash readback, then requires HTTP 206 with the original immutable version/ETag, exact suffix Content-Range and Content-Length, and the final original full byte count, MD5 and SHA256. Source admission still does not certify gzip EOF or preprocessing identity. The two historical bodies total 4,539,067,770 bytes; the remaining network requirement after source 3's resume is 221,795,203,885 bytes. All original final bodies total 227,610,388,647 bytes. Both preserved extra prefixes total 2,405,433,344 bytes.

The global ledger SHA256 `e2432da81b96aaf22dd12c5fb552f385f119b0ca6dfc1ef4f028bfaf790a401a` reconciles its components to 318,518,274,622 bytes under the 322,122,547,200-byte (300 GiB) actual namespace ceiling, leaving 3,604,272,578 reserved bytes unallocated. It includes both prefixes, planned derivatives, failed/temporary products, other declared namespaces and a 1 GiB metadata allowance. The meter includes regular non-symlink files and AppleDouble sidecars. The physical 5 GiB SSD free floor is a separate guard. The original 3 GiB internal floor, one worker, 2 GiB owned-worker RSS limit, 7,200-second body deadline, 345,600-second family deadline, fixed curl/binary hash and no automatic retry remain unchanged. At audit completion, internal free space was 34,501,054,464 bytes and SSD free space was 1,013,674,147,840 bytes; these are observations, not future guarantees.

## Finalization and ownership

The v5 rejection is preserved, including both independently demonstrated provisional-success persistence failures. V6 creates and file/directory-fsyncs PENDING before family ownership and source work (executor lines 539–544). Its new protected finalizer covers both primary copies, their directory fsyncs, exact current readback hashes, source receipt hashes, the terminal seal, complete bindings/mount/namespace checks and final resource/deadline/deferred-signal checks (lines 457–501). PENDING remains the authoritative veto if any primary, terminal or supplemental failure write fails. Its unlink is the last commit action after all checks; no fallible diagnostic, hash, filesystem query or resource check follows that unlink within finalization.

The final unlink is deliberately not followed by a directory fsync: a crash restoring PENDING conservatively rejects an otherwise completed attempt. This review does not simulate power failure or establish hardware durability. The bound actual-SSD tiny directory-fsync control SHA256 is `b0d93669d5154e691a40fe7dd8eb5fd156133f62af9c5d37b29c13f99699431b`; the new controls also use real tiny file/directory fsync operations.

Source 3 now repeats full bindings, correct mounted-volume identity, global namespace cap, internal/SSD free floors, body/family deadlines and deferred termination checks after prefix copy/readback and immediately before Popen (lines 326–335). Both independently injected post-copy binding and free-space failures preserved the invented prefix and a failed source receipt with zero worker creation.

Eight earlier real helper controls, including owned cleanup fallback, deferred signals, receipt-write failure and early-status failure under retained flock, remain applicable through exact helper AST identity and pinned control-script/executor/receipt hashes. They were not rerun here. `family_lock`, deferred signal registration, cleanup proof, independent fallback, quarantine and safe diagnostics are unchanged. The source-level unconditional ownership cleanup still encloses receipt inspection/persistence; a group whose disappearance is not verified cannot release the family flock through those paths. No current real group was started by this review.

## Independent bounded controls

`reviews/independent_acquisition_v6_metadata_controls_v4.py` uses the exact imported v6 family execution/finalization functions with mocked acquisition and identity prerequisites, real private SSD tiny JSON files, fsyncs and a temporary flock. Its source 3 tests use only invented `abc` bytes and the exact pre-worker acquire control flow. A separately written consumer oracle decides admission from complete artifact state, independently of the executor's internal committed flag.

| Control | Observed terminal behavior |
|---|---|
| Second primary copy fails; every supplemental failure write fails | PENDING retained; provisional first ALL100 rejected |
| Post-terminal free floor fails; every supplemental failure write fails | Both ALL100 copies and seal remain provisional; PENDING vetoes them |
| Terminal seal write fails; every supplemental failure write fails | PENDING retained; both primaries rejected |
| PENDING unlink fails | PENDING retained; failure markers written; rejected |
| Initial PENDING directory fsync fails | No ownership/source work; PENDING retained; rejected |
| Deferred SIGTERM request after terminal persistence | Final signal gate fails; PENDING retained; rejected |
| A primary changes after terminal persistence | Current hash readback fails; PENDING retained; rejected |
| Family deadline expires after terminal persistence | Final deadline gate fails; PENDING retained; rejected |
| Healthy synthetic 100-member completion | Two identical primary copies and exact terminal seal; PENDING removed; admitted |
| Duplicate invocation after healthy completion | Immutable ownership rejects before acquisition; successful artifacts unchanged |
| Source 3 post-copy binding changes | Failed receipt and three-byte invented prefix preserved; Popen count zero |
| Source 3 post-copy internal free floor fails | Failed receipt and three-byte invented prefix preserved; Popen count zero |

All 12 controls passed in 2.43 seconds with 38,289,408-byte peak reviewer RSS. There were zero real curl/native workers and zero real GWAS body reads. Mocked identity and resource prerequisites in these controls are explicitly distinct from the actual frozen metadata checks and future runtime gates. The tests establish these bounded failure behaviors, not scientific reproduction or untested hardware/runtime outcomes.

## Mandatory consumer contract and scope limit

Any subsequent acquisition or pipeline consumer must require the exact terminal seal plus **both** primary copies with their current bound identical hashes, exact reviewed plan/executor identities, 100 ordered original receipt identities and the exact reused-2/new-98 map. It must reject any PENDING marker or any primary/terminal `.failure.json` marker. A primary `ALL100_EXACT_SOURCE_BODIES_ACQUIRED` status alone is insufficient. This review certifies the v6 producer preparation and tests that complete contract using an independent oracle; future consumer implementation is not certified here.

The hashed JSON receipt provides every metadata binding, each control observation, precise qualification and the resource snapshot. Its prior checkpoint body proof remains the authority for the two old bodies and retained source 3 prefix. Root must separately bind this exact review and receipt in an admission, retain all earlier failures, enforce current runtime gates and verify eventual real receipts; this preflight does not launch or declare completion of the 98 pending transfers.
