# LD-component sensitivity with the same markers

For the post-selected FI–insomnia locus 464, we held all 2,741 replayed markers and the summary statistics fixed and changed only the LAVA LD-variance retention threshold. The 99% setting reproduces the prior negative FI variance and 216 components. The other settings show a strong numerical dependence on the retained LD representation:

| LD variance threshold | Components | FI genetic-variance diagonal | FI failed |
|---:|---:|---:|---|
| 95% | 73 | +3.437×10⁻⁶ | No |
| 97% | 107 | +2.080×10⁻⁶ | No |
| 98% | 142 | +1.446×10⁻⁶ | No |
| 99% | 216 | −4.191×10⁻⁹ | Yes |
| 99.5% | 306 | −3.140×10⁻⁷ | Yes |

The inspected code subtracts an estimated sampling term from average projected association signal. The [arithmetic verifier](verify.py) independently reproduces that calculation from the saved quantities. As more components were retained, the average projected signal fell relative to the sampling term in this example. This demonstrates that **LD truncation alone is sufficient to change the sign of this locus's estimate**, even without changing marker support. It does not identify a specific responsible marker, show why the projection behaves this way, or determine which setting would be statistically accurate. Sampling variability, mismatch between GWAS and reference LD, and other issues remain possible explanations.

A [component-band replay](component_bands/RESULTS.md) identifies the added 98–99% components as contributing less average squared FI projection than the estimated sampling term. It explains the arithmetic sign change but does not determine its statistical cause.

This is a single post-selected locus and the thresholds were chosen after seeing instability. Positive variance at a lower threshold is not a validated repair, biological evidence, or permission to select a favorable setting. No local correlation was calculated, the original 99% result is still failed, and none of the 12 failed family gates is rescued. A valid replacement would require a separately frozen method, study-matched LD, overlap handling and family-wide calibration.

The [pre-run protocol](PROTOCOL.md), [R script](run.R), [all numerical outputs](results.json), and [verification](verify.py) document the replay. The 99% run matches the previous saved FI diagonal within 10⁻¹⁵ and uses the same pinned R/LAVA environment and reference described in the [parent replay](../RESULTS.md). This is consistency and arithmetic verification, not independent validation of the historical input provenance or the statistical model.
