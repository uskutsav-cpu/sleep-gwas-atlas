# Retained-failure diagnostic replay for locus 464

We replayed the saved FI–insomnia and FI–sleep-apnea locus calculations using the SSD's R 4.3.3/LAVA 0.1.5 environment. The repository runtime validator passed its version and selected implementation checks. We used the archived settings and seed, except `drop.failed=FALSE` so the failed FI quantities could be inspected; we did not run a bivariate test. Both replays reproduced their archived SNP and component counts and the reported valid univariate heritability values to within 10⁻⁹. This is a diagnostic consistency check, not a validation of all historical inputs, sample overlap, calibration or biological inference.

| Run | Internal SNPs / components | FI estimated genetic-variance diagonal | Implied FI univariate estimate before capping | FI status |
|---|---:|---:|---:|---|
| Insomnia | 2,741 / 216 | −4.1914×10⁻⁹ | −9.0534×10⁻⁷ | Failed; displayed heritability capped at zero |
| Sleep apnea | 2,452 / 185 | +4.8878×10⁻⁷ | +9.0424×10⁻⁵ | Retained; bivariate analysis was not eligible |

The variance diagonal is the quantity checked for failure in the inspected implementation; multiplying it by the retained component count gives the uncapped local estimate here because the reference scaling factor was 1. The replay verifies why FI was dropped in the insomnia calculation: its estimated genetic-variance diagonal was negative. It does **not** show whether sampling variation, LD/input mismatch or another factor produced that value. The capped zero is a computational boundary, not proof of zero heritability.

The replay's internal SNP lists confirm the reconstructed sets exactly: all 2,452 sleep-apnea markers are in the insomnia run, which adds 289. Differences in local marker support and LD decomposition are therefore real in this replay, but comparing the two runs does not isolate their causal effect. The sleep traits and other inputs differ, and the original study-family gates remain failed. Neither result supplies a local genetic correlation or independently replicates Song et al.

The saved descriptor paths contain spaces, which the pinned LAVA reader misparses if unquoted. For this replay only, [replay.R](replay.R) writes quoted copies of the same descriptors and uses their saved absolute file paths directly; the original files were unchanged. A temporary path alias restored the SSD environment's original expected location. [Insomnia](insomnia.json) and [sleep-apnea](sleep_apnea.json) outputs retain the quantities and internal marker lists. [Verification code](verify.py) checks exact set nesting and variance arithmetic and hashes the diagnostic outputs, reference LD file and package code. The other input file hashes and limitations are documented in the [marker audit](../RESULTS.md). The installed package check is not a complete independent software or source-provenance audit.

A [controlled marker-support replay](restricted_support/RESULTS.md) retains the FI–insomnia inputs while restricting to the sleep-apnea marker set. It reproduces the positive FI variance of the apnea run, demonstrating support sensitivity for this locus; it does not validate the smaller support or repair the failed gate.

A [fixed-marker LD truncation diagnostic](pc_threshold/RESULTS.md) also changes the sign of the FI variance when varying only the LD-variance retention threshold. It does not establish the correct threshold or repair any failed gate.

A [60-marker independent-reference comparison](reference_comparison/RESULTS.md) reveals one reproducible out-of-range entry and a small negative eigenvalue in the selected UKB LD submatrix. It does not establish why these appear or whether they caused the FI failure.
