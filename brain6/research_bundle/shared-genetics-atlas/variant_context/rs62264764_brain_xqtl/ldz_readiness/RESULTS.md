# rs62264764 brain-cell QTL: LD readiness diagnostic

The available European 1000 Genomes reference is not yet sufficient to validate a multi-signal brain-cell colocalization. Both LINC03051 cell datasets match all **1,548** markers in the earlier reference LD matrix, and reference versus QTL ALT-frequency differences have a median of **0.0130** and maximum of **0.0624**. However, the LD–association diagnostics show notable residuals in each cell type. Similar frequencies do not establish that the reference LD matches the QTL or GWAS cohorts.

| Trait and reference | Estimated inconsistency parameter | Largest absolute standardized residual | Residuals above 3 | Residuals above 5 |
|---|---:|---:|---:|---:|
| Oligodendrocyte LINC03051 eQTL | 0.0235 | 6.42 | 28 | 4 |
| Excitatory-neuron LINC03051 eQTL | 0.0177 | 5.04 | 25 | 1 |
| Prior insomnia GWAS check | 0.0168 | 6.60 | 35 | 2 |
| Prior ADHD GWAS check | 0.0236 | 11.27 | 22 | 6 |

The last two rows come from the [earlier adrenal-region pilot](../../rs62264764_dense_qtl/multisignal_pilot/RESULTS.md), using the same reference marker order and GWAS statistics. The 3 and 5 cutoffs are descriptive, post-hoc flags, not prespecified or calibrated pass/fail thresholds. Outliers can reflect reference mismatch, sampling, allele problems, actual signals or model misspecification; they do not identify their cause. We did not run a new multi-signal coloc or use these diagnostics to declare the previous single-signal models correct or incorrect.

All 1,548 reference keys and GWAS beta/SE values matched each brain-cell dataset after ALT alignment. The portal source records `num_samples=737` for all matched variants in both cell types; we used 737 as the one sample-size input for these diagnostics. That source field is not independently established as the appropriate effective N for the pseudo-bulk association model. The reference contains 503 European-labeled individuals and is singular on this many markers; it is not study-matched. The two brain cell datasets share ROSMAP donors, and their identical frequency comparisons are not independent evidence. The excitatory regional bytes were retrieved by indexed ranges rather than authenticated as a complete chromosome file, as described in its [retrieval report](../excitatory_coloc/RESULTS.md).

A [public ADSP reference assessment](adsp_public_reference/REFERENCE_ASSESSMENT.md) identifies larger thresholded pairwise LD files and controlled-access full-matrix resources. The public pairwise files do not provide a complete matrix; zero-filling missing correlations is not a validated repair. A [selected-pair lookup](adsp_public_reference/RESULTS.md) reports near-perfect tagging in the larger panel, but is not a full matrix or disease replication.

A [post-selected comparison of public ADSP LD pairs](adsp_public_reference/selected_pair_comparison/RESULTS.md) found noticeable differences from the 503-person reference for some reported pairs. The public file suppresses weaker pairs, so it is not a full-matrix replacement or evidence of the cause of the LD–Z flags.

A [public GCTB UK Biobank block check](gctb_ukb_block/RESULTS.md) matched 1,547 variants and verified both streamed ZIP members, but the normalized extracted matrix had negative eigenvalues and failed its pre-specified numerical gate. It was not used for a new LD–Z or fine-mapping fit.

In a subsequent [post-failure sensitivity analysis](gctb_ukb_block/REPAIR_SENSITIVITY_RESULTS.md), two separately specified numerical repairs reduced some residual summaries, but ADHD still showed a substantial mismatch with both repaired matrices. This does not validate the references or change the failed local-correlation gate.

The [protocol](PROTOCOL.md), [matching code](prepare.py), [matching receipt](matching_receipt.json), [diagnostic code](diagnose.R), [summary](ldz_summary.tsv), [variant-level diagnostic values](ldz_details.tsv) and [verification](verification.json) document the check. Independent arithmetic reproduced the standardized residuals and flag counts; that is not a validation of the kriging model, true sample sizes, or error rates. Matched GWAS-derived data remain on the SSD. A defensible multi-signal analysis requires study-matched or well-validated LD, correct effective sample sizes, stronger allele and input-model checks and appropriate statistical review. The failed LAVA gates remain failed.
