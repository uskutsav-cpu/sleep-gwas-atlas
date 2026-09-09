# Track B local-analysis production handoff

Current data status: `READY`. Current compute status: `READY_FOR_RAM_BENCHMARK`.

The three frozen pairs, eight required analysis traits, 2,495 loci, LDSC overlap submatrix, conditional covariates, correction families, and input hashes are prepared and locked. The official LAVA UK Biobank European v1.1 LD payload is present and sealed.

Current free space at preparation: 17686564864 bytes (16.472 GiB).
Policy minimum for checksum-ledgered acquisition and extraction: 37580963840 bytes (35 GiB).
Physical memory at preparation: 8589934592 bytes (8.000 GiB).
Fingerprint-bound RAM evidence: `RAM_BENCHMARK_REQUIRED` (0 selected intact phase representatives; maximum observed peak NA GiB).
No assumed monolithic RAM minimum is used for the production decision; admission preserves a 1073741824 byte (1 GiB) non-worker reserve.

## Minimum production environment

- x86_64 Linux recommended
- R 4.3.x with LAVA 0.1.5
- enough RAM for the largest measured whole-locus worker; execute one fresh R process per locus
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
python3 scripts/lava_contract.py --verify-reference
python3 scripts/129_prepare_track_b_lava_chromosome_inputs.py
python3 scripts/129_prepare_track_b_lava_chromosome_inputs.py --verify
python3 scripts/119_track_b_lava_contract.py --preflight
python3 scripts/131_run_track_b_lava_sequential.py --benchmark
python3 scripts/131_run_track_b_lava_sequential.py --run
python3 scripts/121_validate_track_b_lava.py
```

Do not run the existing 45-trait/396-pair LAVA result family as a substitute for Track B. The dedicated Track B supervisor consumes `results/track_b/lava_pair_manifest.tsv`, uses all 2,495 complete loci, runs each locus in a fresh R process, and applies BH FDR only after collating the full frozen family. Conditional tests run in a second per-locus pass after the discovery-family BH gate and use only the separately frozen BMI-only, sleep-apnea-only, and MDD-only models in `local_conditional_manifest.tsv`.

No local result exists yet. No empty table or synthetic output is presented as science.
