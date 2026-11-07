# Brain-cell LINC03051 colocalization: insomnia-source sensitivity

On the same 2,441 matched variants, changing the insomnia source from the original Brain6 input to FinnGen registry insomnia materially changes the **single-signal** posterior. The FinnGen fit is dominated by the hypothesis that the expression association is present but the insomnia association is not detectable in this window. It therefore does not resolve whether the molecular and disease signals share a causal variant.

| Insomnia source | Shared prior p12 | H1: eQTL only | H3: both, distinct | H4: shared |
|---|---:|---:|---:|---:|
| Original Brain6 | 10⁻⁵ | 0.00291 | 0.99708 | 0.00000892 |
| Original Brain6 | 10⁻⁶ | 0.00291 | 0.99709 | 0.000000892 |
| FinnGen R13 | 10⁻⁵ | 0.76550 | 0.11619 | 0.11832 |
| FinnGen R13 | 10⁻⁶ | 0.85673 | 0.13003 | 0.01324 |

H0 and H2 are negligible at the displayed precision. These are model-dependent posterior probabilities, not p-values or probabilities of gene mediation. The apparently larger FinnGen H4 is **not evidence that the sources demonstrate different causal biology**: its evidence for an insomnia signal is much weaker and H4 changes substantially with the prior. The two models share the same expression data, and the region and candidate support were chosen after earlier results. A one-signal model can miss secondary shared signals.

At rs62264764, FinnGen's A-allele beta is −0.02248 (SE 0.00813; source p = 0.00569); the original Brain6 input has beta −0.03999 (SE 0.00726). These are marginal effects, not evidence that this marker is causal. The variant's A frequency is 0.1343 in the selected expression data and 0.2145 in FinnGen. Across matched variants the median absolute frequency difference is 0.0210 and the 95th percentile is 0.0984; the maximum is 0.1662. This flags population/frequency differences for attention, not an allele-orientation diagnosis or a formal ancestry comparison. The excitatory-neuron expression sample comes from aging-brain donors, not a sleep-disease experiment.

A [reported-statistic and prior sensitivity](statistic_sensitivity/RESULTS.md) found a noticeable QTL p-versus-beta/SE discrepancy at strong markers. A subsequent t-distribution check explains the discrepancy numerically; the tested substitutions and disease-prior settings retain an expression-only-dominated FinnGen fit, but do not validate its likelihood or a causal conclusion.

## Methods and limitations

The earlier [excitatory-neuron dataset](../RESULTS.md) provided 2,458 already selected Brain6–QTL matches. We matched 2,441 to the FinnGen R13 F5_INSOMNIA source by GRCh38 position and nonpalindromic allele pair, retaining exactly the same markers in both fits. Seventeen had no unique matching FinnGen record. The region may truncate relevant cis signals, and this restricted set is not a complete FinnGen or QTL fine-mapping analysis. Raw beta/SE rather than source p-values determine the ABFs; source association tests can have different p-value constructions. Expression prior SD 0.15 and case-control prior SD 0.20 are assumptions; p1 = p2 = 10⁻⁴. Sample overlap and the case-control effect-prior scale are not fully verified.

The [protocol](PROTOCOL.md), [matching code](prepare.py), [matched table](matched.tsv), [source receipt](matching_receipt.json), [publisher index](../../../../PUBLICATION.md#external-dependencies) (external target: `finngen.tbi`), [R analysis](run.R), [posteriors](posteriors.tsv), and [independent arithmetic verification](verification.json) preserve the check. The complete local FinnGen file's SHA-256 matched its previous receipt; its index MD5 matched publisher metadata. A separate Python implementation reproduced all four posterior vectors to within 5 × 10⁻¹⁶. The underlying excitatory QTL was a partial remote extraction, not a fully hashed chromosome file. This sensitivity does not validate AlphaGenome output, identify LSAMP or LINC03051 as the effector, establish two-trait replication or repair the failed LAVA QC gates.

Data for this study were prepared, archived and distributed by NIAGADS at the University of Pennsylvania (U24-AG041689), funded by the National Institute on Aging. We acknowledge ADSP FunGen xQTL Atlas contributors and participants; accession **NG00184**.
