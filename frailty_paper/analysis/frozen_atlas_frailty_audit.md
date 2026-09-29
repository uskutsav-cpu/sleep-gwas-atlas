# Reuse audit: frozen sleep–frailty global LDSC results

Audit date: 2026-09-23. This is a read-only audit/extract of the already-frozen core atlas, not a new global-rg LDSC run. Recreate and verify it from the repository root with `python frailty_paper/scripts/11_audit_frozen_atlas_frailty.py`; publication is idempotent and non-identical existing outputs are refused.

## Frozen evidence

`results/atlas/core.provenance.json` declares atlas-v1.0 complete for its core downstream outputs. SHA-256 checks passed for `results/atlas/traits.tsv`, `results/atlas/trait_pairs.tsv`, and `results/tables/rg_matrix.tsv` against that provenance. The locked family has 396 sleep-by-non-sleep pairs, of which 372 are primary. The frozen Frailty Index row is European ancestry, GRCh37/hg19, N=175,226, with source accession GCST90020053. Its reported h2 is 0.1093 (SE 0.0050; Z=21.86), LDSC intercept 1.02, ratio 0.0509, and Phase 1 verdict PASS. The raw `data/raw/frailty.txt.gz` passed a full gzip read and its recomputed SHA-256 is `a6702bc1c70bd1f2f3ec650d192dcdad72f7596bb716f8e6ef6c72359923e693`, matching the previously registered checksum.

The frozen all-396-family q-values support 9 of the 12 sleep–frailty global rg rows at q<=0.05. The 12 unchanged rows are extracted to `frozen_atlas_frailty_global_rg.tsv`; their original estimates and q-values are retained. This is the result family for potential reuse in the manuscript, subject to the provenance limitation below.

## Heritability reproduction and remaining limitation

The current versioned `frailty_v1` workspace contains a separate harmonized and munged FI pipeline output. We reran LDSC h² from that munged input using the repository LDSC 3.0.1 environment and the EUR LD-score reference. The rerun produced observed-scale h²=0.1093 (SE=0.0050; Z=21.86), intercept=1.0203 (SE=0.0086), ratio=0.0509, and 1,167,050 SNPs after reference merges. These match the frozen-atlas FI h² values at the precision reported. The command log is `../results/frailty_v1/logs_fi_verification/h2_frailty_reproduction_2026-09-23.log`; the manifest records hashes for the source, versioned harmonized/munged inputs, LDSC code, all 44 chromosome LD-score/M reference files, and the run log in `frailty_index_h2_reproduction_manifest_2026-09-23.json`.

This is an independent rerun of the h² calculation on the current versioned frailty pipeline input, not an exact replay of the original frozen-atlas build. The legacy provenance expects `results/tables/h2_summary.tsv` with SHA-256 `1c5c2a2e0b1716fba8c907f7e57492db1cf24b1f0d0d052c2d4518fd3482be67` and a harmonized FI file with SHA-256 `2b18346e024040ee8cbebade0888292cc971e5f1656591351c1698af88626e3c`; these exact artifacts are not present at their legacy paths. The current versioned harmonized FI hash differs. Thus, the frozen summary outputs and h² value are corroborated, while exact legacy preprocessing and the original atlas summary-table byte-level provenance remain unreproduced. Do not label this a complete end-to-end reproduction of atlas-v1.0.

## Interpretation boundary

These are global genetic correlations from the frozen atlas. They do not establish local sharing, pleiotropic loci, a shared causal variant, colocalization, molecular mechanism, or causal direction. UK Biobank sample overlap is expected for most sleep traits and is quantified only at the cohort level in `sample_overlap_assessment.tsv`; exact participant intersections remain unknown.

## Current source-chain recheck (2026-09-24 01:00 UTC)

Revalidated every path in `frailty_index_h2_reproduction_manifest_2026-09-23.json`: all six declared inputs (raw FI, versioned harmonized and munged FI, LDSC runner/core, and rerun log) plus all 44 chromosome LD-score/M reference files are present and match their recorded byte counts and SHA-256 values (50/50; zero mismatches). The versioned FI h² rerun remains 0.1093 (SE 0.0050), with 1,167,050 SNPs after reference merges. This is read-only manifest verification, not a new LDSC run, and it does not reproduce the absent legacy atlas harmonized file or h² summary.
