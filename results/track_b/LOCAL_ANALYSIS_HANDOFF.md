# Track B local-analysis production handoff

Current status: `BLOCKED_BY_DATA` and `BLOCKED_BY_COMPUTE`.

The three frozen pairs, eight required analysis traits, 2,495 loci, LDSC overlap submatrix, conditional covariates, correction families, and input hashes are prepared and locked. The official LAVA UK Biobank European v1.1 LD payload is absent.

Current free space at preparation: 217157632 bytes (0.202 GiB).
Policy minimum for checksum-ledgered acquisition and extraction: 37580963840 bytes (35 GiB).

## Minimum production environment

- x86_64 Linux recommended
- R 4.3.x with LAVA 0.1.5
- at least 16 GiB RAM; 32 GiB recommended for parallel or conditional work
- at least 35 GiB free for reference acquisition; 60 GiB recommended for outputs/checkpoints
- the exact repository commit and ignored dense/munged inputs whose hashes are in `local_analysis_input.lock.json`

## Exact execution order

```bash
python3 scripts/110_track_b_checkpoint.py --verify
python3 scripts/111_freeze_track_b_pairs.py --verify
python3 scripts/113_freeze_track_b_replication_sources.py --verify
python3 scripts/115_build_track_b_dense_qc.py --verify
python3 scripts/116_prepare_track_b_local_inputs.py --verify
bash scripts/32_download_lava_reference.sh --download
python3 scripts/lava_contract.py --verify-reference --rehash
```

Do not run the existing 45-trait/396-pair LAVA result family as a substitute for Track B. A Track B runner must consume `results/track_b/lava_pair_manifest.tsv`, use all 2,495 loci, run univariate h2 for all eight predeclared traits, and apply BH FDR across every actually tested locus row for the three frozen pairs. Conditional tests must use only the covariates in `local_conditional_manifest.tsv`.

No local result exists yet. No empty table or synthetic output is presented as science.
