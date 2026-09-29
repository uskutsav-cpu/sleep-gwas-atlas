# PLACO LD factor reconstruction diagnostic

**Status: diagnostic only; no LD values or candidate/locus results promoted.**

The preserved exception table contains 12,321 pair-specific rows representing 12,317 unique chromosome/SNP edges on chromosomes 2, 5, 7, 11, 13, 14, and 17. Every raw value remains excluded from LD clumping because it lies outside the valid correlation range. The existing pinned LAVA 0.1.5 replay reproduced these raw values; this separate audit asks whether the official UK Biobank v1.1 format-14 eigen-decomposed factor products reproduce them and what a unit-diagonal rescaling would look like for diagnosis.

## Findings

- All 12,317 unique exception edges lie within primary eigen-decomposed blocks. The stored float32 eigenvectors/eigenvalues reproduce each raw edge with maximum absolute difference `5.107e-15` from the archived LAVA replay inputs.
- The factor-product diagonal at exception endpoints ranges from `0.9996840901` to `1.0017012812` (median `1.0000042858`), while the LAVA reader sets the returned matrix diagonal to one.
- Dividing each factor product by the square root of its endpoint factor diagonals yields diagnostic absolute values in `[0.9992110145, 1.0]`; none exceed one. All 12,317 raw edges have absolute value at least `0.999`.

This is consistent with a diagonal scale mismatch in low-rank reconstruction, but it does not prove that this is the source-reference defect or establish that rescaled edges are suitable for clumping. The rescaled values are not used for analysis. The native replay remains authoritative evidence that the out-of-range values are emitted by LAVA; no raw-genotype reconstruction was available to independently establish the correct correlations. The underlying source of the scale behavior is therefore unresolved.

## Reproduction and provenance

The authoritative code-bound output pair for this diagnostic is `qc/placo_ld_factor_normalization_v3.tsv` and `qc/placo_ld_factor_normalization_v3.provenance.json`; v1 and v2 exploratory writes are retained but are not used by current-output validation. Reproduce with `python brain6/scripts/audit_placo_ld_factor_normalization_v1.py --output brain6/qc/placo_ld_factor_normalization_reproduction.tsv --provenance brain6/qc/placo_ld_factor_normalization_reproduction.provenance.json`. The script verifies the archived exception and replay checksums, verifies the sizes and SHA-256 values of each affected reference `.info` and `.bcor` file against the earlier replay provenance, decodes the format-14 block indexes, checks edge-level factor-product agreement, records its own code hash, and refuses to overwrite prior diagnostic outputs. The upstream format-level implementation is [LAVA `load_ld.cpp`](https://github.com/josefin-werme/LAVA/blob/main/src/load_ld.cpp); current upstream source does not replace the pinned 0.1.5 runtime replay as evidence about the historical run.

This diagnostic does not change the four-pair partial PLACO status, the absence of protected Track B, the 12,317-edge exclusion, or the downstream gates for independent loci, fine-mapping, and molecular follow-up.
