# Reproduce the saved-summary workaround locally

For the separate, SSD-dependent LDSC rerun, use `python3 replay_ldsc.py --output-dir /absolute/path/to/a/new-directory`. The script requires the retained source paths and original LDSC environment documented inside it. It refuses an existing destination, checks archived common-support input hashes, records reference/code hashes, reruns nine fits, and compares all point estimates, SEs and 27 block-delete vectors with the saved originals. It preserves failures and partial output. The [6 November replay receipt](ldsc_replay_20261106/receipt.json) records a completed matching run; [verification notes](VERIFICATION.md) define its scope. This route starts with retained munged inputs, rather than original GWAS production or raw-summary ingestion.

The [reproduction driver](reproduce.py) reconstructs the nine-moment sampling covariance, refits the joint and residual-covariance models, checks identification ranges, recalculates the equal-weight index and weight sensitivity, and draws the scientific figure. It needs Python with NumPy and base R (`Rscript` on PATH). It requires neither SciPy, additional R packages, network access nor the external SSD.

From this analysis directory, supply a **new** output directory whose parent already exists:

```sh
python3 reproduce.py --output-dir /absolute/path/to/new-output-directory
```

The driver refuses an existing destination, copies only the required local logs, block-delete values, receipts and scripts, and runs in that destination. It preserves the original results and records input hashes and software details in `replay_receipt.json`. Subprocess outputs and errors are retained as logs; failures stop the run without deleting its files. The index, identification and weight-sensitivity JSON must exactly match the saved numerical outputs before the driver reports success. Inspect the regenerated joint and residual parameter tables and diagnostics as well.

## Verified local run

On 6 November 2026, the run at `/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/work/physical-index-replay-20261106/` completed successfully. All three numerical JSON outputs matched exactly. A separate comparison confirmed that `parameters.txt`, `diagnostics.txt`, `residual_parameters.txt` and `residual_diagnostics.txt` also matched byte for byte. The [reproduction receipt](../../PUBLICATION.md#external-dependencies) (external target: `../../../../work/physical-index-replay-20261106/replay_receipt.json`) records hashes and environment information.

This establishes reproducibility of calculations from the **saved LDSC summaries**. It does not rerun source GWAS ingestion or LDSC, authenticate individual-level phenotype construction, resolve sample overlap or LD-reference limitations, provide independent statistical validation, reproduce a matched exhaustion-deleted frailty GWAS, or train/validate the ML prototype. The original external source files remain necessary for end-to-end source reproduction. Author requests remain unsent pending explicit authorization.
