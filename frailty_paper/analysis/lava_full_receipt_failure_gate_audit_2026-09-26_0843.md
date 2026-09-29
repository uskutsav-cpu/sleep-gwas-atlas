# FI×sleep LAVA full receipt and failure classification — 2026-09-26 08:43 UTC

Read-only full-family scan completed 2026-09-26T08:43:47–08:43:52 UTC.

- 11,281 receipts were validated; the inventory grew 11,281→11,283 while the coordinator advanced from 11,281 to 11,282. No malformed receipts, duplicate identities, claim issues, or audit mutations were found.
- The coordinator requested six workers but kept two effective under the existing swap safeguard (2,772/3,072 MiB, at least 90%). PIDs 78548 and 78549 were alive at audit time, with distinct atomic claims. The analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failures: 2,148 total: 1,906 all-phenotype negative-variance failures and 242 no-specified-SNP-in-reference failures. All 12 traits fail the locked 1% process/univariate gate.

## Failure interpretation

The worker is pinned to R 4.3.3 and LAVA 0.1.5 at commit `e729a245f7b6923967a96804fbf5246eadf2d6c6`, verified by the runtime validator in the runner. At that exact LAVA source, `process.locus` returns `NULL` when both phenotypes have negative local variance estimates (`drop.failed=TRUE`); when only one phenotype fails, it drops that phenotype and retains the other. Thus the 1,906 failures are expected LAVA method outcomes, but they still count as untested/failing slots under the frozen completeness gate. The frozen runner settings remain unchanged. Primary source: [pinned LAVA input processing implementation](https://github.com/josefin-werme/LAVA/blob/e729a245f7b6923967a96804fbf5246eadf2d6c6/R/input_processing.R).

Of the 242 no-reference-SNP failures, 240 occur at chr6 loci 950–969 for every sleep trait. These are the frozen MHC-excluded intervals (chr6:25,684,630–33,864,262); a direct read-only comparison confirmed no prepared FI SNP IDs in that interval despite 43,819 reference SNPs there. The two remaining failures are chr9 locus 1484 (140,097,760–141,146,682) for accel sleep duration and chronotype. Earlier provenance records attribute the MHC absence to the locked QC exclusion; no inputs or thresholds were changed to resolve these failures.

## Live progress

The 08:43 coordinator snapshot was 11,282/29,940. At 08:44:40 it had advanced to 11,293/29,940 with two active workers and zero launch failures or stale recoveries. Since the 07:17:51 checkpoint at 10,227 receipts, the observed rate is about 736 receipts/hour; 18,647 remain, giving a provisional estimate of about 25.3 hours at that measured rate. Six-worker scaling remains unsafe given swap saturation and the previously recorded six-worker resource trial; the current two-worker safeguards remain in place.

No analysis settings, reference data, hashes, worker behavior, or scientific QC criteria were changed. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted by the failed family gates.
