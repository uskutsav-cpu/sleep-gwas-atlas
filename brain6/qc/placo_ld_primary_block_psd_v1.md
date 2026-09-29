# PLACO primary-block PSD diagnostic

**Diagnostic only. No LD values or candidate/locus results were modified or promoted.**

This analysis extends the exception-edge factor reconstruction to all 19 primary eigen-decomposed LD blocks (250 SNPs each) containing the 12,317 unique out-of-range PLACO edges. For each block it reconstructs the stored factor Gram matrix, a LAVA-style matrix with the same off-diagonal entries and the diagonal fixed to one, and a diagonal-normalized Gram matrix. All source `.info`/`.bcor` files are verified against the checksums used by the archived pinned-LAVA replay and the prior edge audit.

## Findings

- All 19 blocks have positive stored eigenvalues (minimum across blocks `0.00964597`; maximum `223.9961`); factor ranks range from 16 to 107.
- The raw factor Gram diagonals range from `0.985396` to `1.022450`. The factor Gram matrices and diagonal-normalized matrices have minimum eigenvalues within floating-point error of zero, as expected for their low-rank positive-semidefinite constructions.
- Replacing each block's factor diagonal with one while leaving its off-diagonal products unchanged makes all 19 block matrices non-PSD. Their minimum eigenvalues range from `-0.022055` to `-0.00001258` (median `-0.002045`). Across these blocks, 29,304 raw off-diagonal factor products exceed `|r|=1`; none do after factor-diagonal normalization.

This demonstrates that the unit-diagonal substitution is sufficient to create the observed PSD/range pathology in the stored-factor representation. It strongly supports the proposed scale-mismatch mechanism for LAVA's returned matrices. It still does not establish the genotype-level LD that should be used for PLACO clumping: normalized low-rank reference factors are not a substitute for raw-genotype reconstruction or a vendor-validated corrected reference. Original raw edges remain excluded.

## Reproduction

The authoritative code-bound output pair is `qc/placo_ld_primary_block_psd_v2.tsv` and `qc/placo_ld_primary_block_psd_v2.provenance.json`. Reproduce with `python brain6/scripts/audit_placo_ld_primary_block_psd_v1.py --output brain6/qc/placo_ld_primary_block_psd_reproduction.tsv --provenance brain6/qc/placo_ld_primary_block_psd_reproduction.provenance.json`. The script verifies prior factor-audit provenance and reference hashes, records its own script hash, and refuses to overwrite existing outputs. This audit does not change the four-pair partial PLACO status or downstream eligibility gates.
