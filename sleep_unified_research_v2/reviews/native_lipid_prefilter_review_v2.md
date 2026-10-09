# Lightweight native lipid prefilter review, v2

Verdict: **the corrected wrapper is acceptable for the scoped native rsID prefilter; 3/3 native prefilter receipts are independently verified.** No native worker or scientific analysis was launched by this reviewer.

Reviewed wrapper SHA-256: `3962ba4fc9b932cf83048ee9b06fe5dbca0886636659b4ecffeb37ad9fd7569d`. The actual native script, current HEAD blob and original-main 659d01cf blob are all `aae8f23d5c417641d529287c8e7509ea788455e3`. Native script SHA-256 is `13e05a44d28d4433f4109d7ed67c812de6bfd4c7e94a4f2fd91cf0e7477db36d`. The actual HapMap3 allowlist SHA-256 is `ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed`, matching all three historical QC records.

The wrapper gates exact acquired-source SHA against historical prefilter source hash, checks the original allowlist and native-script identity, uses recovered Python 3.9, freezes a per-trait command/resource plan, and binds source/allowlist/native-code/Python-binary/historical-QC/acquisition-receipt hashes before and after execution. Original data and sealed v1 outputs are read; the new per-trait SSD output directory must not preexist.

Two material guard faults were corrected before launch: the native script writes output.tmp.<PID> until its final atomic rename, so monitoring only final output missed growing bytes; the revised wrapper includes known temporary/final output and provenance paths and a post-exit size/floor/time check. The revised termination handler escalates a 30-second terminate timeout to kill/wait for only the owned worker, then preserves a failure receipt. The current monitor also tolerates FileNotFoundError during the worker's atomic temporary-to-final rename.

The lightweight plan requires 512 MiB internal free space before launch, stops below 128 MiB internal free space, observes worker RSS against 768 MiB, observes combined worker output/provenance against 128 MiB, permits one worker and 3,600 seconds, and places temporary source/output work on SSD. The full-native LDSC 3 GiB internal guard remains unchanged. Sampled 2-second stop thresholds can briefly overshoot and are not strict OS-enforced quotas; exact expected output is far below the source size. The orchestrator must keep trait launches sequential.

| Trait | Historical source rows | Historical retained rows | Historical prefilter gzip SHA-256 | Receipt verification |
|---|---:|---:|---|---|
| hdl | 46150908 | 1223518 | `a68d17d076d2f07e972b9c994df471b8078155f0d021516f6ff3f8fc03864d93` | exact native bytes/counts verified |
| ldl | 47006483 | 1223471 | `86ed65e10e0b6aa4a31da13b8cea35ab39c1964305741c87ba400d63c5c2b56e` | exact native bytes/counts verified |
| triglycerides | 47196261 | 1223515 | `e4c8fb0e31cd443a906d2d293005086738a19b6ebbaecf480fc3f787e4211ae1` | exact native bytes/counts verified |

Success requires unchanged bound hashes, zero exit/no resource stop, native provenance, exact historical source/retained counts, allowlist SHA and gzip SHA. This reviewer independently hashes completed gzip outputs with 64 KiB reads and checks provenance/plan receipts; pending files are not counted as success.

This reproduces one preparatory rsID row-selection step only. It does not complete raw harmonization, MAF/INFO screening, allele/effect/build QC, HapMap3 munging, heritability, pairwise LDSC, separate 396/1,200 testing families, biological inference or independent replication. Hash-identical gzip is a reproducible data-processing result, not a new sleep-genetics finding.

The JSON contains actual source/QC/reference/code receipts and per-result checks. Re-run `python3 sleep_unified_research_v2/reviews/native_lipid_prefilter_review_v2.py` after execution receipts appear. Code/report/JSON are sealed by the matching .sha256 file.
