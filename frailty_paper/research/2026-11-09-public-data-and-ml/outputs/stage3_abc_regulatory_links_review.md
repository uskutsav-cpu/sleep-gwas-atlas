# Predicted enhancer–gene links at exact candidate positions

I intersected the 22 exploratory markers plus rs323509 with the public [Nasser et al. 2021 Activity-by-Contact (ABC) predictions](https://doi.org/10.1038/s41586-021-03446-x) in 131 cell types and tissues. The source includes element–gene predictions with ABC score ≥0.015. Two candidate positions overlapped a source interval and returned two predicted links:

| Candidate | Predicted target label | Source context | ABC score | Element class |
|---|---|---|---:|---|
| rs7734045 | `LOC389332` | H1-derived cultured neuronal progenitor cells | 0.0471 | Intergenic |
| rs7731706 | `TRPC7-AS2` | Fetal large intestine | 0.0220 | Genic |

These are **predicted regulatory links**, not evidence that the candidate allele changes element activity, that the target mediates insomnia or frailty, or that the cultured/fetal context is the disease-relevant cell type. Target labels are retained exactly as supplied; stable IDs and current naming have not been independently resolved. The ABC score is a model score, not a probability that the gene causes disease. The H1 neuronal-progenitor model is a cultured-cell context, not direct adult-brain evidence.

The [complete matched rows](stage3_abc_candidate_regulatory_links.csv) include source interval, target label, tissue context and all source fields. The [receipt](stage3_abc_candidate_regulatory_receipt.json) records the 7,717,392 source rows scanned, full-download local SHA-256, coordinates and counts; the provider checksum was not independently verified. The source resource is the thresholded Nasser 2021 `AllPredictions.AvgHiC.ABC0.015.minus150` file on the [Engreitz Lab resource page](https://www.engreitzlab.org/resources/). Its intervals are treated as GRCh37 and 0-based half-open; the [AnnotSV source documentation](https://www.lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf) identifies this same file as GRCh37, while the source coordinates follow interval conventions. 1-based candidate positions from the existing project were converted for intersection. Candidate builds and representations still need independent authentication.

This is an **exact-position** screen; it did not examine LD proxies, credible-set variants, unlisted biosamples, predictions below the source threshold or other regulatory models. No overlap for the other 21 candidates is not proof that they lack regulatory activity or target genes. The source interval and gene links are model predictions, not an Ensembl enhancer-to-gene annotation or experimental validation. The network and candidate table retain these edges separately from QTL and colocalization evidence.
