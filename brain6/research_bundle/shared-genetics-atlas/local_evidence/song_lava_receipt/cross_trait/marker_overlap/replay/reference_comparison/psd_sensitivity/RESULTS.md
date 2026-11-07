# Sensitivity to numerical LD-matrix adjustments

We repeated the post-selected FI–insomnia locus 464 calculation on the same 2,741 markers, using the original 99% LAVA retention rule and fixed inputs. Two different in-memory transformations made the extracted approximate reference matrix positive semidefinite within numerical tolerance. Neither transformation reconstructs the underlying genotypes or establishes a correct study-matched LD matrix.

| Matrix used | Relative Frobenius change | Retained components | FI variance diagonal | FI status |
|---|---:|---:|---:|---|
| Original extracted matrix | 0 | 216 | −4.191×10⁻⁹ | Failed |
| Clip negative eigenvalues, rescale diagonal | 0.00103 | 216 | +5.406×10⁻⁹ | Not dropped in this diagnostic |
| Add diagonal ridge, rescale diagonal | 0.10155 | 2,456 | −4.849×10⁻⁶ | Failed |

The unchanged calculation reproduced the earlier replay. For eigenvalue clipping, the largest absolute matrix-entry change was 0.0372; for ridge it was 0.1030. The ridge was set to 0.1133559 before rescaling, based on the magnitude of the original minimum eigenvalue. The adjusted matrices' minimum eigenvalues were approximately −3.9×10⁻¹⁴ and +9.0×10⁻⁹, respectively, compatible with positive semidefiniteness to numerical precision. Clipping changes the near-zero FI estimate's sign; the ridge alters the eigenvalue distribution and causes the 99% rule to retain far more components.

These outcomes show that the selected locus is **sensitive to how the LD matrix is modified**. They do not establish which matrix or FI estimate is valid, prove that the original matrix defect caused the failure, or justify choosing the positive result. The two transforms change the data differently and were assessed only at this post-selected locus. The original 99% calculation still fails, all frozen local family gates remain failed, and no bivariate local correlation was estimated. A defensible replacement needs study-appropriate LD and an independently specified, calibrated family-wide workflow.

The [pre-run protocol](PROTOCOL.md), [R code](run.R), [results](results.json), and [arithmetic verification](verify.py) document the test. The wrapper substitutes each matrix only within this R process and restores the original function afterward; no source reference files were changed. The verifier checks saved arithmetic and PSD tolerances, not independent reference validity or inference calibration. Full matrices were not saved, and the [parent replay's provenance limitations](../../RESULTS.md) still apply.
