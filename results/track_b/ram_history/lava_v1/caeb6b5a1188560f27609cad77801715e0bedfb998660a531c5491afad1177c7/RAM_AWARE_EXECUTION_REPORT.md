# Track B RAM-aware execution report

This is an operational execution report, not a scientific result. Physical memory visible to the supervisor is 8589934592 bytes (8.000 GiB). The predeclared non-worker reserve is 1073741824 bytes (1.000 GiB). Measurements are OS-reported peak RSS around one fresh analysis process; no 16-GiB monolithic assumption is used.

## Scientifically equivalent decomposition decisions

- **LAVA — safely decomposable by complete predeclared locus.** Each discovery worker loads one exact chromosome summary-statistic shard and the complete reference chromosome, constructs one intact official locus, runs all eight local univariate tests, and runs every locally eligible frozen pair. All 2,495 locus results are collated before the unchanged family-wide BH corrections. The conditional pass starts only after that full discovery family is frozen and runs BMI-only, sleep-apnea-only, and MDD-only models separately. Failed and unstable loci remain explicit; no skipped locus is called null.
- **PLACO+ — the global nuisance/correlation fit is not decomposable.** It must use all valid genome-wide variants. Only single-variant testing after those global parameters are frozen may be chunked, and only under the exact same equations and final global multiple-testing family.
- **Fine-mapping and coloc — safely decomposable by validated whole locus, never within an LD block.** Every unit retains the complete frozen SNP set and full ancestry-matched LD matrix.
- **Indivisible units.** One complete LAVA locus, the PLACO+ genome-wide nuisance fit, and one whole fine-map/coloc LD locus are the smallest permitted memory units.

## Measured resource results

| Analysis | Scientifically completed units | Max peak RAM (GiB) | Median runtime (s) | Total observed runtime (h) |
|---|---:|---:|---:|---:|
| LAVA_DISCOVERY | 2432 | 0.470 | 5.31 | 4.080 |
| LAVA_DISCOVERY_AGGREGATION | 0 | NA | NA | NA |

Current LAVA assessment: **PROVISIONAL DISCOVERY PASS — conditional representative measurement remains pending.** Fingerprint-bound benchmark state: `DISCOVERY_REPRESENTATIVES_COMPLETE_CONDITIONAL_PENDING`. Maximum observed admission-relevant peak (all constructed units plus every attempted frozen benchmark candidate): 0.470 GiB. Expected sequential runtime from currently completed units: discovery ~3.68 h.

The benchmark selection is deterministic and frozen before each phase: discovery uses only sealed reference chromosome/locus SNP burden; conditional benchmarking is restricted only by the already-predeclared full-family FDR and conditioner-h2 eligibility gates, then uses the same input-burden ranking. Scientific failures are preserved but cannot serve as a RAM-pass representative. Every worker has an immutable receipt/log/result bundle committed by one atomic ready link.

## Remaining larger-server decisions

- LAVA needs a larger server only if a complete intact locus plus the reserve actually exceeds this machine, or if a non-locus aggregation/finalization unit cannot be made equivalently streaming. It is not blocked by an assumed monolithic requirement.
- PLACO+ remains dependent on the measured global all-variant nuisance/correlation fit; later exact variant shards do not replace that global fit.
- Fine-mapping/coloc requires a larger server only for a complete indivisible locus whose measured peak plus reserve exceeds this machine.

Machine-readable measurements are in `results/track_b/RAM_BENCHMARK.tsv`; their current fingerprint binding is in `results/track_b/RAM_BENCHMARK.provenance.json`.
