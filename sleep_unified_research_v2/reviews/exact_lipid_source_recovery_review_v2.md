# Exact lipid source recovery: bounded independent technical review, v2

Verdict: **the revised gates and frozen resource plan are adequate for the current three pinned source identities; acquisition pending or incomplete; no three-source recovery claim yet.**

Reviewed acquisition script SHA-256: `8e627294ee22fdf5931c8b9ea8a78c0fb49999ebc534ef1dc4bae01e28c1bf8a`. Frozen resource plan SHA-256: `554a6a740cdbf553c1a101ab3dd9f75c63c548250be54f2ac672c340521d3a7a`. No acquisition/native function was executed by this reviewer. The Python review uses independent 64 KiB-buffer hashing and writes only its own files in `sleep_unified_research_v2/reviews/`.

All 16 metadata/receipt/registry checks pass. The registry is byte-identical to the frozen checkpoint (`4b194d25641c2308a84428a7d127969e2eeec0b1884626ab0b8f5a49e1b7cc93`); it contains exactly HDL, LDL and log-triglycerides, all PUBLIC DIRECT_GZIP sources. Current HEAD receipts contain three distinct matching trait/URL/size/SHA identities, successful curl return codes and exact Content-Length values. Their raw header hashes were checked independently. HEAD is a preflight gate; it is not a content-identity test.

The frozen plan permits one worker, no body retries, 6,844,892,917 total source bytes, and a 5,368,709,120-byte SSD reserve. SSD free space before acquisition was 1,131,610,636,288 bytes. The source bodies go directly to a new SSD directory; the script hashes with 4 MiB reads. Actual curl is 8.7.1, and the revised command enforces a per-source --max-filesize cap. curl documents in-transfer size-limit enforcement from 8.4.0; both the installed local manual and [official curl manual](https://curl.se/docs/manpage.html#--max-filesize) support the version gate. This is a source-storage plan, not clearance for later preprocessing/LDSC memory or internal-disk use.

Two material faults identified in the first inspected script were corrected before acquisition: an empty or unrelated HEAD list could pass the original all() test, and expected retained bytes were not enforced during GET. The frozen current plan binds the revised script, registry, original HEAD receipt and original filename-search receipt by SHA-256. The earlier HEAD/search receipts retain their older script SHA, which truthfully identifies their generating version.

The original search receipt is preserved with its inaccessible `/Volumes/Extreme SSD/.TemporaryItems` entry. The amended gate accepts only that path and the recorded permission error; other errors and any lipid filename candidate still stop acquisition. The new plan explicitly qualifies the inaccessible namespace. There are no lipid filename candidates in the scanned readable namespaces. A separate named-archive receipt reports no target candidates. Neither proves that no differently named source bytes exist elsewhere.

Original source locations and frozen v1 outputs are only read. JSON evidence uses exclusive creation, and existing destination source/partial files cause an abort. Successful files are admitted only when curl returns zero, actual bytes exactly equal the expected count, and streaming SHA-256 equals the frozen source hash; a mismatch is preserved and blocks downstream use. The destination retains canonical archive filenames in a new directory rather than recreating old raw aliases.

Residual hardening items do not invalidate the current acquisition: (1) the generalized length/set gate still allows three duplicate registry and corresponding HEAD rows if both are changed together; an isolated predicate counterexample confirms this, but the real frozen registry/HEAD are unique; (2) download-header text lacks a preexisting-file guard before --dump-header, unlike source/partial/JSON files. No header collision was observed. Explicit uniqueness assertions and an exclusive header-path check would make those preservation promises hold in these edge cases.

| Trait | Acquisition receipt | Independent acquired-file hash |
|---|---|---|
| hdl | failed receipt preserved | not performed |
| ldl | pending | not performed |
| triglycerides | pending | not performed |

Acquisition receipt success can support exact historical upstream byte recovery. It does not demonstrate a raw-to-harmonized-to-munged rerun, allele/build/coding correctness, native LDSC estimates, all 396/1,200 results, biological validity or independent replication. No scientific reproduction is inferred from source hashes.

Run this read-only audit with `python3 sleep_unified_research_v2/reviews/exact_lipid_source_recovery_review_v2.py`; add `--hash-acquired` only after final source receipts exist and resource headroom permits sequential source reads. The JSON contains file paths/hashes, deterministic gate counterexample and receipt details. The code/report/receipt are sealed by the matching .sha256 file.
