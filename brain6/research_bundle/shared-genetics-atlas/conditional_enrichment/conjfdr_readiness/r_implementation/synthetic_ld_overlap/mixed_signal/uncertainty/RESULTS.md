# Simulation uncertainty: post-hoc paired check

The 20 completed mixed-signal datasets yield these approximate 95% Monte Carlo t intervals for mean **marker** false-discovery proportion (FDP), based on ten independent seeds per setting:

| Setting | Residual correlation 0 | Residual correlation 0.3 |
|---|---:|---:|
| Smoothed | 0.0702 (0.0439–0.0965) | 0.0945 (0.0695–0.1196) |
| Unsmoothed | 0.0703 (0.0441–0.0964) | 0.0940 (0.0690–0.1191) |

These intervals quantify simulation sampling uncertainty under this artificial, post-selected model; they are approximate with only ten bounded observations and do not measure real-data FDR. The interval at rho=0 includes the nominal 0.05 level. The rho=0.3 interval lies above it under the t approximation, but this is not a universal validation or a formal claim across models or configurations. Marker FDP is distinct from the separately reported block FDP.

The same seeds generated paired datasets for both rho values. The observed mean paired rho=0.3 minus rho=0 change was **0.0243** (approximate 95% interval 0.0047–0.0440) for smoothed and **0.0238** (0.0049–0.0427) for unsmoothed. This controlled simulation contrast does not isolate sample overlap in actual studies or identify a correction; changing residual correlation also changes the generated statistics under this one model. The model and comparison were post-selected.

Neither max-p control had a false marker in ten simulations per setting. For each setting and control, a two-sided exact binomial interval for the *probability a run has at least one false marker* is **0 to 0.3085**. This is not an interval for marker FDR, and zero observed false selections does not prove zero risk. A t interval for the FDP is omitted for these degenerate all-zero samples; reporting [0,0] would falsely imply certainty.

[Protocol](PROTOCOL.md), [code](analyze.R), [per-setting output](per_setting.tsv), [paired output](paired.tsv), and [run log](run.log) document the calculations and verify the 80 unique method–dataset rows against the existing summary. The underlying data and simulation assumptions are described in the [parent report](../RESULTS.md); no new simulation was run and the original LAVA quality gates remain failed.
