# Exact lipid ranged-transfer review: separate v2 addendum

Verdict: **the frozen ranged-transfer design passes source-recovery preflight; all three source receipts are complete.** No material source-admission fault was found.

Reviewed script SHA-256: `9e6eadd08d7a1ece5a6ad7cf6cf27dfd2ebacf1da4001c739043b81a2f23263a`. Frozen resource plan SHA-256: `cec636526e997e4dea84a89d054095431ad26fdd1ec29265ef650781920435cb`. This reviewer did not import or run the acquisition script and launched no downloads/native analyses. All 19 independent small-receipt/preflight checks pass.

The three selected traits are unique and match the exact frozen source registry. HEAD identities and raw-header SHA-256 are bound into preflight. The preserved 1 MiB HDL probe has HTTP206, exact Content-Range/Content-Length and the original HEAD ETag. The failed sequential partial remains 225,443,840 bytes with no admitted SHA; old source files, partials and sealed v1 outputs are preserved.

The separate frozen plan permits eight curl workers, 64 MiB chunks and 600 seconds per chunk. Expected new source network bytes are 6,844,892,917; retained source/chunk/final/previous-partial/probe footprint is 13,916,278,250 bytes, with a 5 GiB reserve against 1,129,735,782,400 SSD bytes free. Internal source bytes are zero. Hash/assembly buffers are 4 MiB; this independent reviewer optionally rehashes final files with 64 KiB reads. Chunk header/receipt overhead fits within the reserve; this does not authorize native-analysis resource use.

Each chunk uses If-Match plus a byte range and a per-range max-filesize cap on curl>=8.4. Admission requires successful curl, exact received bytes, HTTP206, exact Content-Range including full source size, matching Content-Length and original ETag. The fresh destination, exclusive JSON creation and prior chunk/header/body checks prevent ordinary rerun overwrite. Every passing chunk has a body/header hash receipt.

Independent interval construction produces 34 HDL, 35 LDL and 35 triglycerides chunks with no gaps or overlaps and exact total coverage. Assembly sorts by start, compares the full interval sequence, checks summed bytes, rehashes each chunk against its receipt, then writes via exclusive creation. Final admission requires exact full-file size and the original historical SHA-256. ETag/protocol checks are preconditions; final SHA-256 establishes content identity. Retained chunks support auditability and explain the two-copy storage budget.

A failed future does not cancel already queued transfers: ThreadPoolExecutor's context waits for submitted work. Thus failure can continue bounded chunk acquisition until pending workers finish. Each transfer remains capped by chunk bytes/time and eight workers; failure exits before assembly/admission. This is an operational limitation, not a false source-admission pathway. A later plan could cancel unstarted futures for faster failure termination.

The permitted SSD temporary-namespace search gap remains explicitly recorded; filename/one-archive searches do not prove all possible source bytes absent. Success would restore upstream archive bytes in a new directory, not missing aliases or completed raw-to-munged/native inference. No biological finding, full 396/1,200-family reproduction, or independent replication follows from these hashes.

| Trait | Final receipt review | Independent full-file hash |
|---|---|---|
| hdl | exact source/chunk metadata passes | original SHA matches |
| ldl | exact source/chunk metadata passes | original SHA matches |
| triglycerides | exact source/chunk metadata passes | original SHA matches |

The corresponding JSON records paths/hashes, independent coverage counts and acquisition state. Re-run `python3 sleep_unified_research_v2/reviews/exact_lipid_ranged_review_v2.py` after receipts arrive; `--hash-acquired` additionally verifies completed full-file bytes. Code/report/JSON are sealed by the matching .sha256 file.
