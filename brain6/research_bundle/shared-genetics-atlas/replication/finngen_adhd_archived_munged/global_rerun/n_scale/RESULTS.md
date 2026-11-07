# Sample-size scale changes the reported ADHD heritability

Every one of the **1,156,825** rows in the archived ADHD-derived file has N=**21,986.309212**. The historical source record lists 5,559 cases and 489,493 controls; the common formula 4×cases×controls/(cases+controls) gives **21,986.30921196**, matching the file's value to its displayed precision. This strongly suggests a constant effective-sample-size convention, but the raw source, the actual per-variant N and preprocessing remain unverified; the numerical match is not source authentication.

A controlled software sensitivity changed **only** that N column to the historical total of 495,052 in a synthetic copy. SNP IDs, alleles and Z values were checked row by row. No source statistics were changed.

| LDSC output | Archived N=21,986.31 | Synthetic N=495,052 |
|---|---:|---:|
| Standalone h² label printed by LDSC | 0.1892 (SE 0.0238) | 0.0084 (SE 0.0011) |
| Cross-trait ADHD h² label printed by LDSC | 0.1904 (SE 0.0250) | 0.0085 (SE 0.0011) |
| Cross-trait genetic correlation | 0.3817 (SE 0.0521) | 0.3817 (SE 0.0521) |

The sample-size ratio is about 22.5. In this specific constant-rescaling experiment the printed h² and its SE scale approximately inversely with N, while the correlation and SE agree at displayed precision. This demonstrates why the original **0.1892 must be described as the program's N-dependent output**, rather than interpreted as a validated observed-scale disease heritability. The synthetic result is equally **not** a corrected estimate: using total N does not by itself establish the right binary-trait model, scale, prevalence or source convention. No liability-scale conversion is justified. The fact that r_g is unchanged here does not prove its robustness to different per-SNP Ns, ascertainment, sample overlap, allele errors or ancestry mismatch.

The [protocol](PROTOCOL.md), [copy and full row verification](build.py), [receipt](build_receipt.json), [run commands](run.sh), [standalone log](adhd_h2_totalN.log), and [cross-trait log](rg_totalN.log) preserve the sensitivity. The synthetic copy is included solely for reproducibility; it is not a new source or a replication dataset. The original data's provenance and local LAVA gates remain unresolved.
