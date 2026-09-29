# MDD2025 balanced-equivalent sensitivity arm: live launch receipt

The original 111-locus MDD2025 pilot and canonical LAVA outputs remain unchanged. A separate source-model sensitivity arm was frozen in `brain6/config/lava_multitrait_feasibility_v1/mdd2025_balanced_equivalent_pilot_v1.json` (SHA-256 `226d9cc24c1713967ddfed561e49cdf4a3e93d93036625d75dccd1ba7ba5e772`) and prepared after all original source-shard hashes matched the prior materialization receipt. New sumstats files are read-only symlinks to those verified shards; new input-info files encode a 0.5 modeled case fraction.

The exact launch command from the repository root was:

```sh
python3 brain6/scripts/run_mdd2025_balanced_pilot_v1.py run
```

The SSD output root is `/Volumes/Extreme SSD/brain6-work/lava-multitrait-feasibility-v1/mdd2025_balanced_equivalent_pilot_v1`. Its `results/launch.json` has SHA-256 `45a232c8e3cdd7cb40e418396de550d5774da520393bde38dc09b3eb7da409ec` and records four R worker PIDs: 78484, 78485, 78486, and 78491. Worker logs are `results/worker_1.log` through `worker_4.log`; expected terminal files are one TSV and summary JSON per worker, followed by `results/exit.json` and the separate model-sensitivity adjudication. The frozen partitions are 28, 28, 28, and 27 loci.

At the initial checkpoint the four logs were receiving SNP-harmonization messages, with no terminal shards yet. Do not infer TESTED counts from these messages or start a duplicate launch. If interrupted, inspect all files and exact hashes first; preserve completed worker outputs and the launch receipt. After all four exit successfully, run:

```sh
python3 brain6/scripts/run_mdd2025_balanced_pilot_v1.py adjudicate
```

This trait-only arm cannot promote candidates or authorize a full family rescue. The frozen family remains 3,720 `NOT_RUN` against a maximum of 873, and the unresolved long-sleep source contributes 1,291.
